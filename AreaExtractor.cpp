#include <osmtools/AreaExtractor.h>

namespace osmtools {
namespace detail {
namespace AreaExtractor {

bool MultiPolyResolver::closedPolysFromWays(const std::vector<RawWay> & ways, std::vector<RawWay> & closedWays) {
	using std::swap;
	typedef std::unordered_set<uint32_t> ActiveWaysContainer;
	ActiveWaysContainer activeWays; //contains only ways with at least two nodes
	bool allOk = true;
	std::vector<RawWay> resultWays;
	for(uint32_t i = 0, s = ways.size(); i < s; ++i) {
		if (ways[i].size() > 1)
			activeWays.insert(i);
	}
		
	while (activeWays.size()) {
		RawWay activeWay = ways[*activeWays.begin()];
		activeWays.erase(activeWays.begin());
		for(ActiveWaysContainer::iterator it(activeWays.begin()), end(activeWays.end()); activeWay.front() != activeWay.back() && it != end;) {
			bool matched = true;
			const RawWay & otherWay = ways[*it];
			//now check if we can append this way to our current active way
			if (activeWay.back() == otherWay.front()) {
				activeWay.insert(activeWay.end(), otherWay.begin()+1, otherWay.end());
			}
			else if (activeWay.back() == otherWay.back()) {
				activeWay.insert(activeWay.end(), otherWay.rbegin()+1, otherWay.rend());
			}
			else if (activeWay.front() == otherWay.front()) {
				RawWay tmp;
				tmp.reserve(activeWay.size()+otherWay.size());
				tmp.insert(tmp.end(), otherWay.rbegin(), otherWay.rend());
				tmp.insert(tmp.end(), activeWay.begin()+1, activeWay.end());
				swap(activeWay, tmp);
			}
			else if (activeWay.front() == otherWay.back()) {
				RawWay tmp;
				tmp.reserve(activeWay.size()+otherWay.size());
				tmp.insert(tmp.end(), otherWay.begin(), otherWay.end());
				tmp.insert(tmp.end(), activeWay.begin()+1, activeWay.end());
				swap(activeWay, tmp);
			}
			else {
				matched = false;
			}
			
			if (matched) {
				activeWays.erase(it);
				it = activeWays.begin();
			}
			else {
				++it;
			}
		}
		//check if way is an area and closed
		if (activeWay.size() > 2 && activeWay.front() == activeWay.back()) {
			resultWays.push_back(activeWay);
		}
		else {
// 			std::cout << "Unclosed way detected, discarding" << std::endl;
			allOk = false;
		}
	}
	
	closedWays.swap(resultWays);
	
	return allOk;
}
	
bool MultiPolyResolver::multiPolyFromWays(const std::vector<RawWay> & innerIn, const std::vector<RawWay> & outerIn, std::vector<RawWay> & innerOut, std::vector<RawWay> & outerOut) {
	using std::swap;
	if (!outerIn.size()) {
		std::cout << "MultiPolyResolver::multiPolyFromWays: outerIn is empty" << std::endl;
		return false;
	}
	
	std::vector<RawWay> outer;
	bool ok = closedPolysFromWays(outerIn, outer);
	swap(outer, outerOut);
	
	if (innerIn.size()) {
		std::vector<RawWay> inner;
		ok = closedPolysFromWays(innerIn, inner) && ok;
		swap(inner, innerOut);
	}
	
	return ok;
}

}}//end namespace detail::AreaExtractor

}//end namespace