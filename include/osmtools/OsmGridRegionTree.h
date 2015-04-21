#ifndef OSM_TOOLS_GRID_REGION_TREE_H
#define OSM_TOOLS_GRID_REGION_TREE_H
#include <sserialize/spatial/GridRegionTree.h>
#include <sserialize/spatial/DistanceCalculator.h>
#include <sserialize/utility/RangeGenerator.h>
#include <sserialize/templated/MMVector.h>
#include <sserialize/utility/printers.h>
#include <osmtools/types.h>
#include <mutex>

namespace osmtools {

///After adding regions to the grid-tree, they are converted to OsmGeoPolygon and OsmGeoMultiPolygon
template<typename TValue>
class OsmGridRegionTree {
public:
	typedef sserialize::spatial::GeoPoint Point;
	typedef TValue value_type;
	typedef std::vector<sserialize::spatial::GeoRegion*> RegionsContainer;
	typedef std::vector<TValue> ValuesContainer;
	typedef osmtools::OsmGeoPolygon GeoPolygon;
	typedef osmtools::OsmGeoMultiPolygon GeoMultiPolygon;
private:
	typedef sserialize::MMVector<sserialize::spatial::GeoPoint> PointDataContainer;
	typedef sserialize::MMVector<osmtools::OsmGeoPolygon> PolygonsContainer;

	template<typename T_GP_TYPE>
	struct ConvertGP {
		static sserialize::spatial::GeoRegion* conv(PointDataContainer & pointsDest, sserialize::spatial::GeoRegion * r) {
			T_GP_TYPE * gp = static_cast<T_GP_TYPE*>(r);
			sserialize::spatial::GeoPoint* pointsBegin = pointsDest.end();
			uint32_t pointsSize = gp->size();
			pointsDest.push_back(gp->begin(), gp->end());
			return new osmtools::OsmGeoPolygon(osmtools::OsmGeoPolygon::PointsContainer(pointsBegin, pointsSize));
		}
	};
		
	template<typename T_GMP_TYPE>
	struct ConvertGMP {
		static sserialize::spatial::GeoRegion* conv(PointDataContainer & pointsDest, PolygonsContainer & polygonsDest, sserialize::spatial::GeoRegion * r) {
			typedef typename T_GMP_TYPE::GeoPolygon MyGeoPolygon;
			T_GMP_TYPE * gmp = static_cast<T_GMP_TYPE*>(r);

			assert(gmp->outerPolygons().size() + gmp->innerPolygons().size() + polygonsDest.size() <= polygonsDest.capacity());


			OsmGeoPolygon* outerBegin = polygonsDest.end();
			OsmGeoPolygon* innerBegin = 0;
			
			for(MyGeoPolygon & gp : gmp->outerPolygons()) {
				sserialize::spatial::GeoPoint* pointsBegin = pointsDest.end();
				uint32_t pointsSize = gp.size();
				pointsDest.push_back(gp.begin(), gp.end());
				polygonsDest.emplace_back(gp.boundary(), osmtools::OsmGeoPolygon::PointsContainer(pointsBegin, pointsSize));
			}
			innerBegin = polygonsDest.end();
			for(MyGeoPolygon & gp : gmp->innerPolygons()) {
				sserialize::spatial::GeoPoint* pointsBegin = pointsDest.end();
				uint32_t pointsSize = gp.size();
				pointsDest.push_back(gp.begin(), gp.end());
				polygonsDest.emplace_back(gp.boundary(), osmtools::OsmGeoPolygon::PointsContainer(pointsBegin, pointsSize));
			}
			
			return new osmtools::OsmGeoMultiPolygon(gmp->size(),
													osmtools::OsmGeoMultiPolygon::PolygonList(outerBegin, gmp->outerPolygons().size()),
													osmtools::OsmGeoMultiPolygon::PolygonList(innerBegin, gmp->innerPolygons().size()),
													gmp->outerPolygonsBoundary(),
													gmp->innerPolygonsBoundary()
			);
		}
	};
	
	class FixedSizeDiagRefiner final {
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
	void flattenRegions() {
		PointDataContainer tmpPoints(sserialize::MM_SHARED_MEMORY);
		PolygonsContainer tmpPolygons(sserialize::MM_SHARED_MEMORY);
		
		PointDataContainer::size_type pointCount = 0;
		PolygonsContainer::size_type polygonsCount = 0;
		for(sserialize::spatial::GeoRegion* r : m_regions) {
			pointCount += r->size();
			if (r->type() == sserialize::spatial::GS_MULTI_POLYGON) {
				if (dynamic_cast<sserialize::spatial::GeoMultiPolygon*>(r)) {
					sserialize::spatial::GeoMultiPolygon* mp = static_cast<sserialize::spatial::GeoMultiPolygon*>(r);
					polygonsCount += mp->innerPolygons().size() + mp->outerPolygons().size();
				}
				else if (dynamic_cast<osmtools::OsmGeoMultiPolygon*>(r)) {
					osmtools::OsmGeoMultiPolygon* mp = static_cast<osmtools::OsmGeoMultiPolygon*>(r);
					polygonsCount += mp->innerPolygons().size() + mp->outerPolygons().size();
				}
				else {
					sserialize::TypeMissMatchException("OsmGridRegionTree");
				}
			}
		}
		tmpPoints.reserve(pointCount);
		tmpPolygons.reserve(polygonsCount);
// 		for(sserialize::spatial::GeoRegion* &r : m_regions) {
		for(uint32_t i(0), s(m_regions.size()); i < s; ++i) {
			sserialize::spatial::GeoRegion* &r = m_regions[i];
			sserialize::spatial::GeoRegion* tmp = r;
			PolygonPointsContainer::size_type prevCount = tmpPoints.size();
			if (i == 13501) {
				std::cout << "da" << std::endl;
			}
			switch (r->type()) {
			case sserialize::spatial::GS_POLYGON:
				if (dynamic_cast<sserialize::spatial::GeoPolygon*>(r)) {
					r = ConvertGP<sserialize::spatial::GeoPolygon>::conv(tmpPoints, r);
				}
				else if (dynamic_cast<osmtools::OsmGeoPolygon*>(r)) {
					r = ConvertGP<osmtools::OsmGeoPolygon>::conv(tmpPoints, r);
				}
				else {
					sserialize::TypeMissMatchException("OsmGridRegionTree");
				}
				break;
			case sserialize::spatial::GS_MULTI_POLYGON:
				if (dynamic_cast<sserialize::spatial::GeoMultiPolygon*>(r)) {
					r = ConvertGMP<sserialize::spatial::GeoMultiPolygon>::conv(tmpPoints, tmpPolygons, r);
				}
				else if (dynamic_cast<osmtools::OsmGeoMultiPolygon*>(r)) {
					r = ConvertGMP<osmtools::OsmGeoMultiPolygon>::conv(tmpPoints, tmpPolygons, r);
				}
				else {
					sserialize::TypeMissMatchException("OsmGridRegionTree");
				}
				break;
			default:
				throw sserialize::TypeMissMatchException("OsmGridRegionTree");
				break;
			}
			assert(tmpPolygons.size() <= polygonsCount);
			assert(r->size() == tmp->size());
			r->recalculateBoundary();
			assert(tmp->size() == tmpPoints.size()-prevCount);
			assert(r->size() == tmp->size());
			
			delete tmp;
		}
		assert(tmpPoints.size() == pointCount);
		assert(tmpPolygons.size() == polygonsCount);
		m_polygonPoints = std::move(tmpPoints);
		m_polygonsContainer = std::move(tmpPolygons);
	}
private:
	sserialize::spatial::GridRegionTree m_grt;
	PointDataContainer m_polygonPoints;
	PolygonsContainer m_polygonsContainer;
	RegionsContainer m_regions;
	ValuesContainer m_values;
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
	///You should only call this prior to calling addPolygonsToRaster
	///polygon ids are invalid afterwards
	template<typename T_COMPARE>
	void sort(T_COMPARE compfunc) {
		sserialize::RangeGenerator rg(0, size());
		std::vector<uint32_t> tmp(rg.cbegin(), rg.cend());
		std::vector<value_type> * valuesPtr = &m_values;
		std::sort(tmp.begin(), tmp.end(), [compfunc, valuesPtr](uint32_t a, uint32_t b) {
			return compfunc(valuesPtr->at(a), valuesPtr->at(b));
		});
		sserialize::reorder(m_values, tmp);
		sserialize::reorder(m_regions, tmp);
	}
	inline void printStats(std::ostream & out) {
		out << "OsmGridRegionTree::printStats--BEGIN\n";
		m_grt.printStats(out);
		std::cout << "#points: " << m_polygonPoints.size() << "=" << sserialize::prettyFormatSize(m_polygonPoints.size()*sizeof(PolygonPointsContainer::value_type)) << "\n";
		std::cout << "#GeoMultiPolygons: " << m_polygonsContainer.size() << "=" << sserialize::prettyFormatSize(m_polygonsContainer.size()*sizeof(PolygonsContainer::value_type)) << "\n";
		out << "OsmGridRegionTree::printStats--END\n";
	}
	inline std::size_t size() const { return regions().size(); }
	inline const RegionsContainer & regions() const { return m_regions; }
	inline const ValuesContainer & values() const { return m_values; }
	inline ValuesContainer & values() { return m_values; }
	///this is thread-safe
	uint32_t push_back(const sserialize::spatial::GeoRegion & p, const value_type & value) {
		std::unique_lock<std::mutex> (m_mtx);
		uint32_t pid = m_regions.size();
		m_regions.push_back(static_cast<sserialize::spatial::GeoRegion*>( p.copy() ));
		m_values.push_back(value);
		return pid;
	}
	///this is thread-safe
	uint32_t push_back(const sserialize::spatial::GeoRegion & p, value_type && value) {
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
	///if you want them in the tree aswell then you should call this again (which will rebuild the tree)
	void addPolygonsToRaster(unsigned int gridLatCount, unsigned int gridLonCount) {
		flattenRegions();
		FixedSizeDiagRefiner refiner(m_refineMinDiag, m_latRefineCount, m_lonRefineCount);
		sserialize::spatial::GeoRect initRect( sserialize::spatial::GeoShape::bounds(m_regions.cbegin(), m_regions.cend()) );
		sserialize::spatial::GeoGrid initGrid(initRect, gridLatCount, gridLonCount);
		typedef sserialize::spatial::GridRegionTree::TypeTraits<FixedSizeDiagRefiner, osmtools::OsmGeoPolygon, osmtools::OsmGeoMultiPolygon> MyTypeTraits;
		m_grt = sserialize::spatial::GridRegionTree(initGrid, m_regions.begin(), m_regions.end(), MyTypeTraits(), refiner);
		m_grt.shrink_to_fit();
	}
	
	///This is thread safe if you do not after calling addPolygonsToRaster()
	template<typename T_SET_TYPE = std::set<uint32_t> >
	void test(const Point & p, T_SET_TYPE & dest) const {
		std::insert_iterator<T_SET_TYPE> inserter(dest, dest.end());
		m_grt.find(p, inserter);
	}
	
	///This is thread safe if you do not after calling addPolygonsToRaster()
	template<typename T_SET_TYPE = std::set<uint32_t> >
	void test(double lat, double lon, T_SET_TYPE & dest) const {
		return test(Point(lat, lon), dest);
	}
	
	///This is thread safe if you do not after calling addPolygonsToRaster()
	template<typename T_SET_TYPE = std::set<uint32_t> >
	T_SET_TYPE test(const Point & p) const {
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