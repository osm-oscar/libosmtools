#ifndef OSM_TOOLS_GRID_REGION_TREE_H
#define OSM_TOOLS_GRID_REGION_TREE_H
#include <sserialize/spatial/GridRegionTree.h>
#include <sserialize/spatial/DistanceCalculator.h>
#include <sserialize/iterator/RangeGenerator.h>
#include <sserialize/utility/printers.h>
#include <osmtools/types.h>
#include <mutex>

namespace osmtools {

///After adding regions to the grid-tree, they are converted to OsmGeoPolygon and OsmGeoMultiPolygon
class OsmGridRegionTreeBase {
public:
	typedef sserialize::spatial::GeoPoint Point;
	typedef std::vector<sserialize::spatial::GeoRegion*> RegionsContainer;
	typedef osmtools::OsmGeoPolygon GeoPolygon;
	typedef osmtools::OsmGeoMultiPolygon GeoMultiPolygon;
private:
	typedef GeoPointStorageBackend PointDataContainer;
	typedef OsmGeoPolygonStorageBackend PolygonsContainer;

	template<typename T_GP_TYPE>
	struct ConvertGP {
		static sserialize::spatial::GeoRegion* conv(PointDataContainer & pointsDest, const sserialize::spatial::GeoRegion * r) {
			const T_GP_TYPE * gp = static_cast<const T_GP_TYPE*>(r);
			sserialize::OffsetType pointsBegin = pointsDest.size();
			uint32_t pointsSize = gp->size();
			pointsDest.push_back(gp->cbegin(), gp->cend());
			return new osmtools::OsmGeoPolygon(osmtools::OsmGeoPolygon::PointsContainer(&pointsDest, pointsBegin, pointsSize));
		}
	};

	template<typename T_GMP_TYPE>
	struct ConvertGMP {
		static sserialize::spatial::GeoRegion* conv(PointDataContainer & pointsDest, PolygonsContainer & polygonsDest, const sserialize::spatial::GeoRegion * r) {
			typedef typename T_GMP_TYPE::GeoPolygon MyGeoPolygon;
			const T_GMP_TYPE * gmp = static_cast<const T_GMP_TYPE*>(r);
			sserialize::OffsetType outerBegin = polygonsDest.size();
			sserialize::OffsetType innerBegin = 0;
			
			for(const MyGeoPolygon & gp : gmp->outerPolygons()) {
				sserialize::OffsetType pointsBegin = pointsDest.size();
				uint32_t pointsSize = gp.size();
				pointsDest.push_back(gp.cbegin(), gp.cend());
				polygonsDest.emplace_back(gp.boundary(), osmtools::OsmGeoPolygon::PointsContainer(&pointsDest, pointsBegin, pointsSize));
			}
			innerBegin = polygonsDest.size();
			for(const MyGeoPolygon & gp : gmp->innerPolygons()) {
				sserialize::OffsetType pointsBegin = pointsDest.size();
				uint32_t pointsSize = gp.size();
				pointsDest.push_back(gp.cbegin(), gp.cend());
				polygonsDest.emplace_back(gp.boundary(), osmtools::OsmGeoPolygon::PointsContainer(&pointsDest, pointsBegin, pointsSize));
			}
			
			return new osmtools::OsmGeoMultiPolygon(gmp->size(),
													osmtools::OsmGeoMultiPolygon::PolygonList(&polygonsDest, outerBegin, gmp->outerPolygons().size()),
													osmtools::OsmGeoMultiPolygon::PolygonList(&polygonsDest, innerBegin, gmp->innerPolygons().size()),
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
		FixedSizeDiagRefiner(double minDiag, uint32_t latCount, uint32_t lonCount);
		~FixedSizeDiagRefiner();
		bool operator()(const sserialize::spatial::GeoRect& maxBounds, const std::vector<sserialize::spatial::GeoRegion*>& rId2Ptr, const std::vector<uint32_t> & sortedRegions, sserialize::spatial::GeoGrid & newGrid) const;
	};
private:
	sserialize::spatial::GridRegionTree m_grt;
	PointDataContainer m_polygonPoints;
	PolygonsContainer m_polygonsContainer;
	RegionsContainer m_regions;
	uint32_t m_latRefineCount;
	uint32_t m_lonRefineCount;
	double m_refineMinDiag;
protected:
	template<typename T_REORDER_MAP>
	void reorderRegions(const T_REORDER_MAP & rm) {
		sserialize::reorder(m_regions, rm);
	}
public:
	OsmGridRegionTreeBase();
	virtual ~OsmGridRegionTreeBase();
	const sserialize::spatial::GridRegionTree & grt() const;
	void clearGRT();
	void clear();
	///@thread-safety no
	void push_back(const sserialize::spatial::GeoRegion * p);
	///"snap" points to the accuracy of sserialize::Static::spatial::GeoPoint
	///This will also recalculate the bbox of all regions
	void snapPoints();
	void printStats(std::ostream & out);
	std::size_t size() const;
	const PointDataContainer & points() const;
	const RegionsContainer & regions() const;
	void setRefinerOptions(uint32_t latRefineCount, uint32_t lonRefineCount, double minDiag);
	///you can still add more regions after calling this but they will not be part of the tree
	///if you want them in the tree aswell then you should call this again (which will rebuild the tree)
	void addPolygonsToRaster(unsigned int gridLatCount, unsigned int gridLonCount);
	
	template<typename T_OUTPUT_ITERATOR1, typename T_OUTPUT_ITERATOR2>
	void test(const Point & p, T_OUTPUT_ITERATOR1 definiteEnclosing, T_OUTPUT_ITERATOR2 candidateEnclosing) const {
		m_grt.find(p, definiteEnclosing, candidateEnclosing);
	}
	
	///This is thread safe if you do not after calling addPolygonsToRaster()
	///Inserts all hit regions into dest
	template<typename T_CONTAINER_TYPE>
	void test(const Point & p, T_CONTAINER_TYPE & dest) const {
		std::insert_iterator<T_CONTAINER_TYPE> inserter(dest, dest.end());
		m_grt.find(p, inserter);
	}
	
	///This is thread safe if you do not after calling addPolygonsToRaster()
	template<typename T_SET_TYPE = std::set<uint32_t> >
	void test(double lat, double lon, T_SET_TYPE & dest) const {
		test(Point(lat, lon), dest);
	}
	
	template<typename T_OUTPUT_ITERATOR>
	void find(const Point & p, T_OUTPUT_ITERATOR & dest) const {
		m_grt.find(p, dest);
	}
	
	///This is thread safe if you do not after calling addPolygonsToRaster()
	template<typename T_OUTPUT_ITERATOR>
	void find(double lat, double lon, T_OUTPUT_ITERATOR & dest) const {
		find(Point(lat, lon), dest);
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
	T_SET_TYPE test(double lat, double lon) const {
		return test(Point(lat, lon));
	}
};

template<typename TValue>
class OsmGridRegionTree:public OsmGridRegionTreeBase {
public:
	typedef TValue value_type;
	typedef std::vector<TValue> ValuesContainer;
public:
	OsmGridRegionTree() {}
	virtual ~OsmGridRegionTree() {}
	void clear() {
		OsmGridRegionTreeBase::clear();
		m_values = ValuesContainer();
	}
	const ValuesContainer & values() const {
		return m_values;
	}
	ValuesContainer & values() {
		return m_values;
	}
	///You should only call this prior to calling addPolygonsToRaster
	///polygon ids are invalid afterwards
	template<typename T_COMPARE>
	void sort(T_COMPARE compfunc) {
		sserialize::RangeGenerator<std::size_t> rg(0, size());
		std::vector<uint32_t> tmp(rg.cbegin(), rg.cend());
		std::vector<value_type> * valuesPtr = &m_values;
		std::sort(tmp.begin(), tmp.end(), [compfunc, valuesPtr](uint32_t a, uint32_t b) {
			return compfunc(valuesPtr->at(a), valuesPtr->at(b));
		});
		sserialize::reorder(m_values, tmp);
		reorderRegions(tmp);
	}
	///this is thread-safe
	uint32_t push_back(const sserialize::spatial::GeoRegion & p, const value_type & value) {
		std::lock_guard<std::mutex> lck(m_mtx);
		uint32_t pid = size();
		OsmGridRegionTreeBase::push_back(&p);
		m_values.push_back(value);
		return pid;
	}
	///this is thread-safe
	uint32_t push_back(const sserialize::spatial::GeoRegion & p, value_type && value) {
		std::lock_guard<std::mutex> lck(m_mtx);
		uint32_t pid = (uint32_t) size();
		OsmGridRegionTreeBase::push_back(&p);
		m_values.emplace_back(value);
		return pid;
	}
private:
	ValuesContainer m_values;
	std::mutex m_mtx;
};

}//end namespace
#endif
