#ifndef OSM_TOOLS_GRID_REGION_TREE_H
#define OSM_TOOLS_GRID_REGION_TREE_H
#include <sserialize/spatial/GridRegionTree.h>
#include <sserialize/spatial/DistanceCalculator.h>
#include <mutex>

namespace osmtools {

template<typename TValue>
class OsmGridRegionTree {
public:
	typedef sserialize::spatial::GeoPoint Point;
	typedef TValue value_type;
private:
	class FixedSizeDiagRefiner {
	private:
		sserialize::spatial::DistanceCalculator  m_dc;
		double m_minDiag;
		uint32_t m_latCount;
		uint32_t m_lonCout;
	public:
		FixedSizeDiagRefiner(double minDiag, uint32_t latCount, uint32_t lonCount) :
		m_dc( std::shared_ptr<sserialize::spatial::detail::DistanceCalculator>(new sserialize::spatial::detail::GeodesicDistanceCalculator()) ),
		m_minDiag(minDiag),
		m_latCount(latCount),
		m_lonCout(lonCount)
		{}
		~FixedSizeDiagRefiner() {}
		inline bool operator()(const sserialize::spatial::GeoRect& maxBounds, const std::vector<sserialize::spatial::GeoRegion*>& rId2Ptr, const std::vector<uint32_t> & sortedRegions, sserialize::spatial::GeoGrid & newGrid) const {
			if (sortedRegions.size() && m_dc.calc(maxBounds.minLat(), maxBounds.minLon(), maxBounds.maxLat(), maxBounds.maxLon()) > m_minDiag) {
				std::vector< uint32_t >::const_iterator it(sortedRegions.cbegin());
				sserialize::spatial::GeoRect myBounds(rId2Ptr[*it]->boundary());
				++it;
				for(std::vector<uint32_t>::const_iterator end(sortedRegions.cend()); it != end; ++it) {
					myBounds.enlarge(rId2Ptr[*it]->boundary());
				}
				newGrid = sserialize::spatial::GeoGrid( myBounds / maxBounds, m_latCount, m_lonCout);
				return true;
			}
			return false;
		}
	};
private:
	sserialize::spatial::GridRegionTree m_grt;
	std::vector<sserialize::spatial::GeoRegion*> m_regions;
	std::vector<value_type> m_values;
	std::mutex m_mtx;
	uint32_t m_latRefineCount;
	uint32_t m_lonRefineCount;
	double m_refineMinDiag;
public:
	OsmGridRegionTree() : m_latRefineCount(2), m_lonRefineCount(2), m_refineMinDiag(250) {}
	~OsmGridRegionTree() {
		for(std::vector<sserialize::spatial::GeoRegion*>::iterator it(m_regions.begin()), end(m_regions.end()); it != end; ++it) {
			delete *it;
		}
	}
	inline void printStats(std::ostream & out) {
		m_grt.printStats(out);
	}
	inline std::size_t size() const { return regions().size(); }
	inline const std::vector<sserialize::spatial::GeoRegion*> & regions() const { return m_regions; }
	inline const std::vector<value_type> & values() const { return m_values; }
	inline std::vector<value_type> & values() { return m_values; }
	///this is thread-safe
	inline uint32_t push_back(const sserialize::spatial::GeoRegion & p, const value_type & value) {
		std::unique_lock<std::mutex> (m_mtx);
		uint32_t pid = m_regions.size();
		m_regions.push_back(static_cast<sserialize::spatial::GeoRegion*>(p.copy()));
		m_values.push_back(value);
		return pid;
	}
	///this is thread-safe
	inline uint32_t push_back(const sserialize::spatial::GeoRegion & p, value_type && value) {
		std::unique_lock<std::mutex> lck(m_mtx);
		uint32_t pid = m_regions.size();
		m_regions.push_back(static_cast<sserialize::spatial::GeoRegion*>(p.copy()));
		m_values.emplace_back(value);
		return pid;
	}
	void setRefinerOptions(uint32_t latRefineCount, uint32_t lonRefineCount, double minDiag) {
		m_latRefineCount = latRefineCount;
		m_lonRefineCount = lonRefineCount;
		m_refineMinDiag = minDiag;
	}
	
	///you can still add more regions after calling this but they will not be part of the tree
	///if you want tehm in the tree aswell then you should call this again (which will rebuild the tree)
	void addPolygonsToRaster(unsigned int gridLatCount, unsigned int gridLonCount) {
		FixedSizeDiagRefiner refiner(m_refineMinDiag, m_latRefineCount, m_lonRefineCount);
		sserialize::spatial::GeoRect initRect( sserialize::spatial::GeoShape::bounds(m_regions.cbegin(), m_regions.cend()) );
		sserialize::spatial::GeoGrid initGrid(initRect, gridLatCount, gridLonCount);
		m_grt = sserialize::spatial::GridRegionTree(initGrid, m_regions.begin(), m_regions.end(), refiner);
		m_grt.shrink_to_fit();
	}
	
	///This is thread safe if you do not after calling addPolygonsToRaster()
	template<typename T_SET_TYPE = std::set<uint32_t> >
	inline void test(const Point & p, T_SET_TYPE & dest) const {
		std::insert_iterator<T_SET_TYPE> inserter(dest, dest.end());
		m_grt.find(p, inserter);
	}
	
	///This is thread safe if you do not after calling addPolygonsToRaster()
	template<typename T_SET_TYPE = std::set<uint32_t> >
	inline void test(double lat, double lon, T_SET_TYPE & dest) const {
		return test(Point(lat, lon), dest);
	}
	
	///This is thread safe if you do not after calling addPolygonsToRaster()
	template<typename T_SET_TYPE = std::set<uint32_t> >
	inline T_SET_TYPE test(const Point & p) const {
		T_SET_TYPE polys;
		test(p, polys);
		return polys;
	}
	
	///This is thread safe if you do not after calling addPolygonsToRaster()
	template<typename T_SET_TYPE = std::set<uint32_t> >
	inline T_SET_TYPE test(double lat, double lon) const {
		return test(Point(lat, lon));
	}
};

}//end namespace
#endif