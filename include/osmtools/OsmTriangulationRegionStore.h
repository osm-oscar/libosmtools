#ifndef LIBOSMTOOLS_OSM_TRIANGULATION_REGION_STORE_H
#include <unordered_map>

#include <sserialize/utility/hashspecializations.h>

#include <osmtools/OsmGridRegionTree.h>
#include <osmtools/TriangulationGridLocater.h>

//CGAL includes
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

#include <assert.h>
#include <thread>
#include <mutex>

namespace osmtools {

/** This class splits the arrangement of regions into cells to support fast point-in-polygon look-ups
  *
  *
  *
  *
  */
class OsmTriangulationRegionStore {
public:
	typedef CGAL::Exact_predicates_exact_constructions_kernel K;
	typedef CGAL::Exact_intersections_tag Itag;
	typedef CGAL::Triangulation_vertex_base_2<K> Vb;
	typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
	typedef CGAL::Triangulation_data_structure_2<Vb,Fb> TDS;
	typedef CGAL::Constrained_Delaunay_triangulation_2<K,TDS,Itag> CDT;
	typedef CGAL::Constrained_triangulation_plus_2<CDT> CDTP;
	typedef CDT Triangulation;
	typedef Triangulation::Point Point;
	typedef osmtools::GridLocator<Triangulation> GridLocator;
	typedef GridLocator::Face_handle Face_handle;
	typedef sserialize::MMVector<uint32_t> RegionListContainer;
	typedef sserialize::CFLArray<RegionListContainer> RegionList;
private:
	struct FaceHandleHash {
		std::hash<double> m_hasher;
		std::size_t operator()(const Face_handle & f) const {
			std::size_t seed = 0;
			for(int i(0); i < 3; ++i) {
				auto p = f->vertex(i)->point();
				hash_combine(seed, CGAL::to_double(p.x()), m_hasher);
				hash_combine(seed, CGAL::to_double(p.y()), m_hasher);
			}
			return seed;
		}
	};
// 	typedef std::unordered_map<Face_handle, uint32_t, FaceHandleHash> FaceCellIdMap;
	typedef CGAL::Unique_hash_map<Face_handle, uint32_t> FaceCellIdMap;
private:
	GridLocator m_grid;
	FaceCellIdMap m_faceToCellId;
	RegionListContainer m_cellLists;
	std::vector<RegionList> m_cellIdToCellList;
	std::vector<uint32_t> m_refinedCellIdToUnrefined;
public:
	OsmTriangulationRegionStore() {}
	~OsmTriangulationRegionStore() {}
	template<typename TDummy>
	void init(OsmGridRegionTree<TDummy> & grt);
	inline uint32_t cellId(double lat, double lon) const;
	inline const RegionList & regions(uint32_t cellId);
};


template<typename TDummy>
void OsmTriangulationRegionStore::init(OsmGridRegionTree<TDummy> & grt) {
	//we first need to find all relevant regions and extract their segments. This sould be possible by just using the extracted regions since
	//we don't do any calculations with our points so segments with the same endpoints should stay the same in different regions
	typedef std::pair<double, double> RawGeoPoint;
	typedef typename OsmGridRegionTree<TDummy>::GeoPolygon GeoPolygon;
	typedef typename OsmGridRegionTree<TDummy>::GeoMultiPolygon GeoMultiPolygon;
	std::unordered_map<RawGeoPoint, uint32_t> gpToId;
	std::vector<Point> pts;
	std::unordered_set< std::pair<uint32_t, uint32_t> > segments;
	
	
	auto handleSinglePolygon = [&gpToId,&pts, &segments](const GeoPolygon * gp) {
		typename GeoPolygon::const_iterator it(gp->cbegin()), prev(gp->cbegin()), end(gp->cend());
		for(++it; it != end; ++it, ++prev) {
			RawGeoPoint itGp = *it;
			RawGeoPoint prevGp = *prev;
			if (!gpToId.count(itGp)) {
				gpToId[itGp] = pts.size();
				pts.push_back(Point(itGp.first, itGp.second));
			}
			if (!gpToId.count(prevGp)) {
				gpToId[prevGp] = pts.size();
				pts.push_back(Point(prevGp.first, prevGp.second));
			}
			segments.insert( std::pair<uint32_t, uint32_t>(gpToId.at(itGp), gpToId.at(prevGp)) );
		}
	};
	std::cout << "OsmTriangulationRegionStore: extracting segments..." << std::flush;
	for(sserialize::spatial::GeoRegion* r : grt.regions()) {
		if (r->type() == sserialize::spatial::GS_POLYGON) {
			const GeoPolygon * gp = static_cast<const GeoPolygon*>(r);
			handleSinglePolygon(gp);
		}
		else if (r->type() == sserialize::spatial::GS_MULTI_POLYGON) {
			const GeoMultiPolygon * gmp = static_cast<const GeoMultiPolygon*>(r);
			for(const GeoPolygon & gp : gmp->outerPolygons()) {
				handleSinglePolygon(&gp);
			}
			for(const GeoPolygon & gp : gmp->innerPolygons()) {
				handleSinglePolygon(&gp);
			}
		}
	}
	std::cout << "done" << std::endl;
	for(const std::pair<uint32_t, uint32_t> & s : segments) {
		assert(s.first < pts.size());
		assert(s.second < pts.size());
	}
	std::cout << "OsmTriangulationRegionStore: creating triangulation..." << std::flush;
	m_grid.tds().insert_constraints(pts.cbegin(), pts.cend(), segments.cbegin(), segments.cend());
	std::cout << "done" << std::endl;
	
	//we now have to assign every face its cellid
	//a face lives in multiple regions, so every face has a unique list of region ids
	//faces that are connected and have the same region-list get the same cellid
	typedef sserialize::CFLArray< sserialize::MMVector<uint32_t> > MyCFLArray;
	

	std::cout << "Setting initial cellids" << std::flush;
	{
		struct Context {
			std::unordered_map<MyCFLArray, uint32_t> cellListToCellId;
			OsmGridRegionTree<TDummy> * grt;
			RegionListContainer * p_cellLists;
			FaceCellIdMap * p_faceToCellId;
			Triangulation::Finite_faces_iterator facesIt;
			Triangulation::Finite_faces_iterator facesEnd;
			std::mutex lock;
		} ctx;
		ctx.grt = &grt;
		ctx.p_cellLists = &m_cellLists;
		ctx.p_faceToCellId = &m_faceToCellId;
		ctx.facesIt = m_grid.tds().finite_faces_begin();
		ctx.facesEnd = m_grid.tds().finite_faces_end();
		
		struct WorkFunc {
			Context * ctx;
			std::unordered_set<uint32_t> tmpCellList;
			Point centroid;
			Face_handle fh;
			WorkFunc(Context * ctx) : ctx(ctx) {}
			WorkFunc(const WorkFunc & other) : ctx(other.ctx) {}
			void operator()() {
				while (true) {
					std::unique_lock<std::mutex> lck(ctx->lock);
					if (ctx->facesIt != ctx->facesEnd) {
						fh = ctx->facesIt;
						++(ctx->facesIt);
						centroid = CGAL::centroid(fh->vertex(0)->point(), fh->vertex(1)->point(), fh->vertex(2)->point());
					}
					else {
						return;
					}
					lck.unlock();
					double x = CGAL::to_double(centroid.x());
					double y = CGAL::to_double(centroid.y());
					tmpCellList.clear();
					ctx->grt->test(x, y, tmpCellList);
					MyCFLArray tmp(tmpCellList.begin(), tmpCellList.end());
					std::sort(tmp.begin(), tmp.end());
					uint32_t faceCellId = 0;
					lck.lock();
					if (!ctx->cellListToCellId.count(tmp)) {
						faceCellId = ctx->cellListToCellId.size();
						sserialize::MMVector<uint32_t>::size_type off = ctx->p_cellLists->size();
						ctx->p_cellLists->push_back(tmp.begin(), tmp.end());
						tmp = MyCFLArray(ctx->p_cellLists, off, tmp.size());
						ctx->cellListToCellId[tmp] = faceCellId;
					}
					else {
						faceCellId = ctx->cellListToCellId.at(tmp);
					}
					(*(ctx->p_faceToCellId))[fh] = faceCellId;
				}
			}
		};
		std::vector<std::thread> threads;
		for(uint32_t i(0), s(std::thread::hardware_concurrency()); i < s; ++i) {
			threads.push_back(std::thread(WorkFunc(&ctx)));
		}
		for(std::thread & x : threads) {
			x.join();
		}

		m_cellIdToCellList.resize(ctx.cellListToCellId.size());
		for(const auto & x : ctx.cellListToCellId) {
			m_cellIdToCellList.at(x.second) = x.first;
		}
	}
	std::cout << "done" << std::endl;
	
	std::cout << "Found " << m_cellIdToCellList.size() << " unrefined cells" << std::endl;
	//now every cell has an id but cells that are not connected may not have different cells
	//we now have to check for each id if the correspondig faces are all connected through cells with the same id
	//this essential is a graph traversel to get all connected components where each face is a node and there's an edge between nodes
	//if their correspondig faces are neighbours and share the same id
	std::cout << "Refining cells" << std::flush;
	{
		FaceCellIdMap tmp;
		std::vector<Face_handle> stack;
		uint32_t cellId = 0;
		for(CDT::Finite_faces_iterator it(m_grid.tds().finite_faces_begin()), end(m_grid.tds().finite_faces_end()); it != end; ++it) {
			Face_handle rfh = it;
			if (tmp.is_defined(rfh)) {
				continue;
			}
			//a new connected component is going to be created
			stack.push_back(rfh);
			while(stack.size()) {
				Face_handle fh = stack.back();
				stack.pop_back();
				if (tmp.is_defined(fh)) {
					continue;
				}
				tmp[fh] = cellId;
				assert(m_faceToCellId.is_defined(fh));
				uint32_t fhId = m_faceToCellId[fh];
				for(int i=0; i < 3; ++i) {
					Face_handle nfh = fh->neighbor(i);
					if (m_faceToCellId.is_defined(nfh) && m_faceToCellId[nfh] == fhId && !tmp.is_defined(nfh)) {
						stack.push_back(nfh);
					}
				}
			}
			
			stack.clear();
			m_refinedCellIdToUnrefined.push_back(m_faceToCellId[rfh]);
			++cellId;
		}
		assert(cellId == m_refinedCellIdToUnrefined.size());
		m_faceToCellId = std::move(tmp);
		std::cout << "done" << std::endl;
		std::cout << "Found " << cellId << " cells" << std::endl;
	}
}


///By definition: items that are not in any cell are in cell 0xFFFFFFFF
uint32_t OsmTriangulationRegionStore::cellId(double lat, double lon) const {
	Face_handle fh = m_grid.locate(lat, lon);
	if (fh->is_valid() && m_faceToCellId.is_defined(fh)) {
		return m_faceToCellId[fh];
	}
	else {
		return 0xFFFFFFFF;
	}
}

const OsmTriangulationRegionStore::RegionList& OsmTriangulationRegionStore::regions(uint32_t cellId) {
	return m_cellIdToCellList.at(m_refinedCellIdToUnrefined.at(cellId) );
}

}//end namespace

#endif