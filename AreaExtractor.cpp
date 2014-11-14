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

generics::RCPtr<osmpbf::AbstractTagFilter> AreaExtractor::setUpMainFilter(ExtractionTypes extractionTypes, bool needsName) {
	generics::RCPtr<osmpbf::AbstractTagFilter> mainFilter;
	generics::RCPtr<osmpbf::OrTagFilter> areaFilter(new osmpbf::OrTagFilter());

	generics::RCPtr<osmpbf::AndTagFilter> wayFilter(new osmpbf::AndTagFilter());
	generics::RCPtr<osmpbf::AndTagFilter> relationFilter(new osmpbf::AndTagFilter());
	
	//setup primitive type filters
	wayFilter->addChild(new osmpbf::PrimitiveTypeFilter(osmpbf::WayPrimitive));
	relationFilter->addChild(new osmpbf::PrimitiveTypeFilter(osmpbf::RelationPrimitive));
	
	if (extractionTypes & ET_BOUNDARIES) {
		areaFilter->addChild(new osmpbf::KeyOnlyTagFilter("boundary"));
	}
	if (extractionTypes & ET_LANDUSE) {
		areaFilter->addChild(new osmpbf::KeyOnlyTagFilter("landuse"));
	}
	if (extractionTypes & ET_AREA) {
		osmpbf::AbstractTagFilter * areaTagFilter = 0;
		std::vector<osmpbf::AbstractTagFilter*> areaExclusions;
		if ((extractionTypes & ET_BUILDING) != ET_BUILDING) {
			areaExclusions.push_back(new osmpbf::BoolTagFilter("building", false));
		}
		if (! extractionTypes & ET_BOUNDARIES) {
			areaExclusions.push_back(new osmpbf::BoolTagFilter("boundary", false));
		}
		if (! extractionTypes & ET_LANDUSE) {
			areaExclusions.push_back(new osmpbf::BoolTagFilter("landuse", false));
		}
		if (areaExclusions.size()) {
			osmpbf::AndTagFilter * andTagFilter = new osmpbf::AndTagFilter();
			andTagFilter->addChild(new osmpbf::KeyOnlyTagFilter("area"));
			andTagFilter->addChildren(areaExclusions.begin(), areaExclusions.end());
			areaTagFilter = andTagFilter;
		}
		else {
			areaTagFilter = new osmpbf::KeyOnlyTagFilter("area");
		}
		areaFilter->addChild(areaTagFilter);
	}
	
	{//set the way filter
		wayFilter->addChild(areaFilter.get());
	}
	
	{//set the relation filter
		osmpbf::AndTagFilter * multiPolyFilter = osmpbf::newAnd(new osmpbf::StringTagFilter("type", "multipoly"), areaFilter.get());
		if (extractionTypes & ET_BOUNDARIES) {
			relationFilter->addChild( osmpbf::newOr(new osmpbf::StringTagFilter("type", "boundary"), multiPolyFilter) );
		}
		else {
			relationFilter->addChild(multiPolyFilter);
		}
	}
	
	if (needsName) {
		mainFilter.reset(osmpbf::newAnd(new osmpbf::KeyOnlyTagFilter("name"), osmpbf::newOr(wayFilter.get(), relationFilter.get())) );
	}
	else {
		mainFilter.reset( osmpbf::newOr(wayFilter.get(), relationFilter.get()) );
	}
	return mainFilter;
}

}//end namespace