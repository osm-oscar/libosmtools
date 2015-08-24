#ifndef LIBOSMTOOLS_OSM_TRIANGULATION_REGION_STORE_H
#define LIBOSMTOOLS_OSM_TRIANGULATION_REGION_STORE_H
#include <unordered_map>

#include <sserialize/algorithm/hashspecializations.h>
#include <sserialize/stats/TimeMeasuerer.h>
#include <sserialize/utility/debug.h>
#include <sserialize/containers/ItemIndexFactory.h>
#include <sserialize/containers/SimpleBitVector.h>

#include <osmtools/OsmGridRegionTree.h>
#include <osmtools/TriangulationGridLocater.h>
#include <sserialize/Static/TriangulationGeoHierarchyArrangement.h>

//CGAL includes
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesher_2.h>

#include <CGAL/Delaunay_mesh_size_criteria_2.h>


#include <assert.h>
#include <thread>
#include <mutex>

namespace osmtools {
class OsmTriangulationRegionStore;

namespace detail {
namespace OsmTriangulationRegionStore {

template<typename T_TRIANGULATION>
typename T_TRIANGULATION::Point centroid(const typename T_TRIANGULATION::Face_handle & fh) {
	return CGAL::centroid(fh->vertex(0)->point(), fh->vertex(1)->point(), fh->vertex(2)->point());
}

//CellTriangulationGraph
class CTGraphBase {
private:
	friend class osmtools::OsmTriangulationRegionStore;
public:
	static constexpr uint32_t NullFace = 0xFFFFFFFF;
	struct FaceNode {
		static constexpr uint32_t NullNeighbor = CTGraphBase::NullFace;
		///by definition: if neighbours[i] == 0xFFFFFFFF => no neighbor
		uint32_t neighbours[3];
	};
	typedef std::vector<FaceNode> NodesContainer;
private:
	uint32_t m_cellId;
	NodesContainer m_nodes;
protected:
	NodesContainer & nodes() { return m_nodes; }
	const NodesContainer & nodes() const { return m_nodes; }
public:
	CTGraphBase() {}
	virtual ~CTGraphBase() {}
	inline uint32_t size() const { return m_nodes.size(); }
	inline uint32_t cellId() const { return m_cellId; }
	///@param bfsTree (nodeId, hopdistance from faceNodeId)
	void calcMaxHopDistance(uint32_t& maxHopDistRoot);
	uint32_t calcDiameter(uint32_t * startNode, uint32_t * endNode);
	const FaceNode & node(uint32_t pos) const { return m_nodes.at(pos); }
	FaceNode & node(uint32_t pos) { return m_nodes.at(pos); }
};

//Graph of the cells of the OsmTriangulationRegionStore

class CellGraph {
private:
	friend class osmtools::OsmTriangulationRegionStore;
public:
	typedef sserialize::MMVector<uint32_t> NodePointersContainer;
	
	class CellNode final {
	private:
		friend class CellGraph;
		friend class osmtools::OsmTriangulationRegionStore;
	private:
		typedef sserialize::CFLArray<CellGraph::NodePointersContainer> NodePointers;
	public:
		typedef NodePointers::const_iterator const_iterator;
		typedef NodePointers::iterator iterator;
	private:
		NodePointers m_d;
	private:
		void rebind(CellGraph::NodePointersContainer * d) { m_d.rebind(d); }
	public:
		CellNode(CellGraph::NodePointersContainer * d, CellGraph::NodePointersContainer::size_type offset, CellGraph::NodePointersContainer::size_type size) :
		m_d(d, offset, size) {}
		CellNode() {}
		~CellNode() {}
		inline uint32_t size() const { return m_d.size(); }
		inline iterator begin() { return m_d.begin(); }
		inline const_iterator begin() const { return m_d.begin(); }
		inline const_iterator cbegin() const { return m_d.cbegin(); }
		
		inline iterator end() { return m_d.end(); }
		inline const_iterator end() const { return m_d.end(); }
		inline const_iterator cend() const { return m_d.cend(); }
		
	};
	typedef std::vector<CellNode> NodesContainer;
private:
	NodePointersContainer * m_nodePtrs;
	NodesContainer m_nodes;
protected:
	NodesContainer & nodes() { return m_nodes; }
	const NodesContainer & nodes() const { return m_nodes; }
public:
	CellGraph() : m_nodePtrs(0) {}
	CellGraph(const CellGraph & other);
	CellGraph(CellGraph && other);
	virtual ~CellGraph();
	CellGraph & operator=(const CellGraph & other);
	CellGraph & operator=(CellGraph && other);

	inline uint32_t size() const { return m_nodes.size(); }
	inline const CellNode & node(uint32_t pos) const { return m_nodes.at(pos); }
	inline CellNode & node(uint32_t pos) { return m_nodes.at(pos); }
	
	sserialize::UByteArrayAdapter & append(sserialize::UByteArrayAdapter & dest, const std::unordered_map<uint32_t, uint32_t> & myIdsToGhCellIds) const;
};


template<typename TDS>
class CentroidDistanceBaseMeshCriteria {
public:
	typedef typename TDS::Face_handle Face_handle;
	typedef typename TDS::Point Point;
	struct Is_bad {
		sserialize::spatial::DistanceCalculator m_dc;
		Is_bad(const sserialize::spatial::DistanceCalculator & dc) : m_dc(dc) {}
		double maxCentroidDist(Face_handle fh) const {
			Point p(::osmtools::detail::OsmTriangulationRegionStore::centroid<TDS>(fh));
			double latp = CGAL::to_double(p.x());
			double lonp = CGAL::to_double(p.y());
			double q = 0.0;
			for(int j(0); j < 3; ++j) {
				Point fp( fh->vertex(j)->point() );
				double lat = CGAL::to_double(fp.x());
				double lon = CGAL::to_double(fp.y());
				double tmp = m_dc.calc(latp, lonp, lat, lon);
				q = std::max<double>(tmp, q);
			}
			return q;
		}
	};
protected:
	sserialize::spatial::DistanceCalculator m_dc;
protected:
	const sserialize::spatial::DistanceCalculator & dc() const { return m_dc; }
public:
	CentroidDistanceBaseMeshCriteria() : m_dc(sserialize::spatial::DistanceCalculator::DCT_GEODESIC_ACCURATE) {}
};

template<typename TDS>
class CentroidDistanceMeshCriteria: public CentroidDistanceBaseMeshCriteria<TDS> {
public:
	typedef CentroidDistanceBaseMeshCriteria<TDS> MyParentClass;
	typedef typename MyParentClass::Face_handle Face_handle;
	typedef typename MyParentClass::Point Point;
public:
	typedef double Quality;
	struct Is_bad: MyParentClass::Is_bad {
		double m_r;
		sserialize::spatial::DistanceCalculator m_dc;
		Is_bad(double maxDist, const sserialize::spatial::DistanceCalculator & dc) : MyParentClass::Is_bad(dc), m_r(maxDist), m_dc(dc) {}
		
		CGAL::Mesh_2::Face_badness operator()(Quality q) const {
			if (q > m_r) {
				return CGAL::Mesh_2::IMPERATIVELY_BAD;
			}
			else {
				return CGAL::Mesh_2::NOT_BAD;
			}
		}
		
		CGAL::Mesh_2::Face_badness operator()(Face_handle fh, Quality & q) const {
			q = MyParentClass::Is_bad::maxCentroidDist(fh);
			return (*this)(q);
		};
	};
private:
	double m_r;
public:
	CentroidDistanceMeshCriteria(double maxDist) : m_r(maxDist) {}
	Is_bad is_bad_object() const { return Is_bad(m_r, MyParentClass::dc()); }
};


template<typename TDS>
class LipschitzMeshCriteria: CentroidDistanceBaseMeshCriteria<TDS> {
public:
	typedef CentroidDistanceBaseMeshCriteria<TDS> MyParentClass;
	typedef typename MyParentClass::Face_handle Face_handle;
	typedef typename MyParentClass::Point Point;
	typedef double Quality;
	struct Is_bad: MyParentClass::Is_bad {
		double m_s;
		TDS * m_tds;
		Is_bad(double maxDist, const sserialize::spatial::DistanceCalculator & dc, TDS * tds) : MyParentClass::Is_bad(dc), m_s(maxDist), m_tds(tds) {}
		
		inline Quality quality(CGAL::Mesh_2::Face_badness fb) const {
			if (fb == CGAL::Mesh_2::NOT_BAD) {
				return m_s;
			}
			else {
				return std::numeric_limits<double>::max();
			}
		}
		
		CGAL::Mesh_2::Face_badness operator()(Quality q) const {
			if (q > m_s) {
				return CGAL::Mesh_2::IMPERATIVELY_BAD;
			}
			else {
				return CGAL::Mesh_2::NOT_BAD;
			}
		}
		
		CGAL::Mesh_2::Face_badness operator()(Face_handle fh, Quality & q) {
			double maxSlope = 0.0;
			double myCD = MyParentClass::Is_bad::maxCentroidDist(fh);
			//by definition: slope goes from nfh to fh
			//=> if this is the smallest triangle, then the slope is negative for all neighbors, if not it's positive
			for(int j(0); j < 3; ++j) {
				Face_handle nfh(fh->neighbor(j));
				if (!m_tds->is_infinite(nfh)) {
					double nfhD = MyParentClass::Is_bad::maxCentroidDist(nfh);
					if (myCD > nfhD) {
						maxSlope = std::max<double>(myCD / nfhD, maxSlope);
					}
				}
			}
			if (maxSlope > 2200) {
				for(int j(0); j < 3; ++j) {
					Face_handle nfh(fh->neighbor(j));
					if (!m_tds->is_infinite(nfh)) {
						double nfhD = MyParentClass::Is_bad::maxCentroidDist(nfh);
						if (myCD > nfhD) {
							maxSlope = std::max<double>(myCD / nfhD, maxSlope);
						}
					}
				}
			}
			q = maxSlope;
			return (*this)(q);
		};
	};
private:
	double m_s;
	TDS * m_tds;
public:
	LipschitzMeshCriteria(double maxSlope, TDS * tds) : m_s(maxSlope), m_tds(tds) {}
	Is_bad is_bad_object() const { return Is_bad(m_s, MyParentClass::dc(), m_tds); }
};

template<typename T_BASE_MESH_CRITERIA>
class RefineTrianglesWithCellIdMeshCriteria: public T_BASE_MESH_CRITERIA {
public:
	typedef T_BASE_MESH_CRITERIA MyParentClass;
	typedef typename MyParentClass::Face_handle Face_handle;
	typedef typename MyParentClass::Point Point;
	typedef typename MyParentClass::Quality Quality;
	struct Is_bad: MyParentClass::Is_bad {
		Is_bad(const typename MyParentClass::Is_bad & base) : MyParentClass::Is_bad(base) {}
		
		CGAL::Mesh_2::Face_badness operator()(Quality q) const {
			return MyParentClass::Is_bad::operator()(q);
		}
		
		CGAL::Mesh_2::Face_badness operator()(Face_handle fh, Quality & q) {
			if (fh->info().hasCellId()) {
				return MyParentClass::Is_bad::operator()(fh, q);
			}
			else {
				q = MyParentClass::Is_bad::quality(CGAL::Mesh_2::NOT_BAD);
				return CGAL::Mesh_2::NOT_BAD;
			}
		};
	};
public:
	RefineTrianglesWithCellIdMeshCriteria(const MyParentClass & base) : MyParentClass(base) {}
	Is_bad is_bad_object() const { return Is_bad(MyParentClass::is_bad_object()); }
};

}}//end namespace detail::OsmTriangulationRegionStore

/** This class splits the arrangement of regions into cells to support fast point-in-polygon look-ups
  *
  *
  *
  *
  */
class OsmTriangulationRegionStore {
public:

	class FaceInfo {
	private:
		//default initialized to UnsetFacesCellId
		uint32_t m_cellId;
	public:
		FaceInfo();
		void clear();
		inline void setCellId(uint32_t cellId) { m_cellId = cellId; }
		inline uint32_t cellId() const { return m_cellId;}
		inline bool hasCellId() const { return m_cellId != OsmTriangulationRegionStore::UnsetFacesCellId; }
	};

	typedef CGAL::Exact_predicates_exact_constructions_kernel K;
	typedef CGAL::Exact_intersections_tag Itag;
// 	typedef CGAL::No_intersection_tag Itag;
// 	typedef CGAL::Triangulation_vertex_base_2<K> Vb;
// // 	typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
// 	typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
// 	typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
// 	typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;
// // 	typedef CGAL::Constrained_triangulation_plus_2<CDT> CDTP;


// 	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef CGAL::Triangulation_vertex_base_2<K> Vb;
	typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, K> FbWithInfo;
	typedef CGAL::Constrained_triangulation_face_base_2<K, FbWithInfo> CFb;
	typedef CGAL::Delaunay_mesh_face_base_2<K, CFb> Fb;
	typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
	typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds, Itag> CDT;

	typedef CDT Triangulation;
	typedef Triangulation::Point Point;
	typedef osmtools::GridLocator<Triangulation> GridLocator;
	typedef GridLocator::Face_handle Face_handle;
	typedef sserialize::MMVector<uint32_t> RegionListContainer;
	typedef sserialize::CFLArray<RegionListContainer> RegionList;
	typedef GridLocator::TriangulationDataStructure::Finite_faces_iterator Finite_faces_iterator;
	typedef GridLocator::TriangulationDataStructure::All_faces_iterator All_faces_iterator;
	typedef sserialize::SimpleBitVector SimpleBitVector;
public:
	class CTGraph: public detail::OsmTriangulationRegionStore::CTGraphBase {
	private:
		friend class OsmTriangulationRegionStore;
	private:
		std::vector<Face_handle> m_faces;
		CGAL::Unique_hash_map<Face_handle, uint32_t> m_faceToNodeId;
	public:
		CTGraph() {}
		virtual ~CTGraph() {}
		Face_handle face(uint32_t faceNodeId);
		uint32_t node(const Face_handle & fh);
		using CTGraphBase::node;
	};
	typedef detail::OsmTriangulationRegionStore::CellGraph CellGraph;
	
	typedef detail::OsmTriangulationRegionStore::CentroidDistanceMeshCriteria<CDT> CentroidDistanceMeshCriteria;
	typedef detail::OsmTriangulationRegionStore::LipschitzMeshCriteria<CDT> LipschitzMeshCriteria;
	//If this tat is selected, then MyRefineTag is mandatory
	typedef detail::OsmTriangulationRegionStore::RefineTrianglesWithCellIdMeshCriteria<LipschitzMeshCriteria> RegionOnlyLipschitzMeshCriteria;
	
	typedef enum { CGALRefineTag, MyRefineTag } RefinementAlgoTags;
	
	//Refines the triangulation if the distance between the centroid of the triangle and any of its defining points is larger than maxDist
public:
	static constexpr uint32_t InfiniteFacesCellId = 0xFFFFFFFF;
	static constexpr uint32_t UnsetFacesCellId = 0xFFFFFFFE;
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
	static Point centroid(const Face_handle & fh);
	//calculate hop-distances from rfh to all other faces of the cell of rfh and store their hop-distance in cellTriangMap and the cells triangs in cellTriangs
	void hopDistances(const Face_handle & rfh, std::vector<Face_handle> & cellTriangs, CGAL::Unique_hash_map<Face_handle, uint32_t> & cellTriangMap, uint32_t & maxHopDist);

	template<typename T_REFINER>
	void cgalRefineMesh(T_REFINER & refiner);

	template<typename T_REFINER, typename T_Dummy>
	void myRefineMesh(T_REFINER & refiner, OsmGridRegionTree<T_Dummy> & grt, uint32_t threadCount);

	void setInfinteFacesCellIds();
	template<typename T_DUMMY>
	void assignCellIds(OsmGridRegionTree<T_DUMMY> & grt, uint32_t threadCount);
private:
	GridLocator m_grid;
	RegionListContainer m_cellLists;
	std::vector<RegionList> m_cellIdToCellList;
	std::vector<uint32_t> m_refinedCellIdToUnrefined;
	bool m_isConnected;
	std::mutex m_lock;
public:
	OsmTriangulationRegionStore();
	OsmTriangulationRegionStore(const OsmTriangulationRegionStore & other) = delete;
	~OsmTriangulationRegionStore() {}
	void clear();
	const Triangulation & tds() const { return m_grid.tds(); }
	Triangulation & tds() { return m_grid.tds(); }
	inline uint32_t cellCount() const { return m_refinedCellIdToUnrefined.size(); }
	inline uint32_t unrefinedCellCount() const { return m_cellIdToCellList.size(); }
	///@param threadCount pass 0 for automatic deduction (uses std::thread::hardware_concurrency())
	///@param meshCriteria must be a modell of CGAL::MeshingCriteria_2
	template<typename TDummy, typename T_TRIANG_REFINER = OsmTriangulationRegionStore::RegionOnlyLipschitzMeshCriteria>
	void init(OsmGridRegionTree<TDummy> & grt, uint32_t threadCount, T_TRIANG_REFINER * meshCriteria = 0, RefinementAlgoTags refineAlgo = MyRefineTag);
	void initGrid(uint32_t gridLatCount, uint32_t gridLonCount);

	void makeConnected();
	///Splits cells into connected cells into smaller cells if they are larger than cellSizeTh
	///This is done in multiple runs where each cell is split into up to numVoronoiSplitRuns smaller cells until each cell is smaller than cellSizeTh
	///cells are not split into equally sized cells but rather by their voronoi diagram
	///refine cells by connectedness so that all cells form a connected polygon (with holes)
	void refineBySize(uint32_t cellSizeTh, uint32_t runs, uint32_t splitPerRun, uint32_t threadCount);
	
	void clearRefinement();
	inline uint32_t unrefinedCellId(uint32_t cellId) { return m_refinedCellIdToUnrefined.at(cellId); }
	uint32_t cellId(double lat, double lon);
	inline uint32_t cellId(const sserialize::spatial::GeoPoint & gp) { return cellId(gp.lat(), gp.lon()); }
	///@thread-safety no
	uint32_t cellId(const Face_handle & fh);
	inline const RegionListContainer & regionLists() const { return m_cellLists; }
	inline RegionListContainer & regionLists() { return m_cellLists; }
	const RegionList & regions(uint32_t cellId);
	Finite_faces_iterator finite_faces_begin() { return m_grid.tds().finite_faces_begin(); }
	Finite_faces_iterator finite_faces_end() { return m_grid.tds().finite_faces_end(); }
	
	CellGraph cellGraph();
	void ctGraph(const Face_handle& rfh, CTGraph& cg);
	void cellInfo(std::vector<Face_handle> & cellRepresentatives, std::vector<uint32_t> & cellSizes);
	
	template<typename T_OUTPUT_ITERATOR>
	void regionCells(uint32_t regionId, T_OUTPUT_ITERATOR out);
	
	void printStats(std::ostream & out);
	
	bool selfTest();
	///serializes to sserialize::Static::spatial::TriangulationRegionArrangement
	sserialize::UByteArrayAdapter & append(sserialize::UByteArrayAdapter & dest, sserialize::ItemIndexFactory & idxFactory);
	
	///serializes to sserialize::Static::spatial::TriangulationGeoHierarchyArrangement
	sserialize::UByteArrayAdapter & append(sserialize::UByteArrayAdapter & dest, const std::unordered_map<uint32_t, uint32_t> & myIdsToGhCellIds);
	
	bool equal(const sserialize::Static::spatial::TriangulationGeoHierarchyArrangement & ra, const std::unordered_map<uint32_t, uint32_t> & myIdsToGhCellIds);
};

template<typename T_REFINER>
void OsmTriangulationRegionStore::cgalRefineMesh(T_REFINER & /*refiner*/) {
// 		CGAL::refine_Delaunay_mesh_2(m_grid.tds(), CGAL::Delaunay_mesh_size_criteria_2<CDT>(0.125, 0.5));
	//this only workd with an epic kernel, but the normal default triangulation needs an epec kernel
// 		CGAL::refine_Delaunay_mesh_2(m_grid.tds(), refiner);
}

template<typename T_REFINER, typename T_Dummy>
void OsmTriangulationRegionStore::myRefineMesh(T_REFINER & refiner, OsmGridRegionTree<T_Dummy> & grt, uint32_t threadCount) {
	uint32_t refineCount = 0;
	std::vector<Point> refinePoints;
	sserialize::MinMax<typename T_REFINER::Quality> qs;
	bool trWasRefined = true;
	typename T_REFINER::Is_bad bo(refiner.is_bad_object());
	typename T_REFINER::Quality q = 0;
	typename T_REFINER::Face_handle rfh;
	for(uint32_t round(0); round < 10000 && trWasRefined; ++round) {
		std::cout << "Trianglerefinement round " << round << std::flush;
		trWasRefined = false;
		qs.reset();
		refinePoints.clear();
		for(Finite_faces_iterator fit(finite_faces_begin()), fend(finite_faces_end()); fit != fend; ++fit) {
			rfh = fit;
			if (CGAL::Mesh_2::NOT_BAD != bo.operator()(rfh, q)) {
				refinePoints.push_back(centroid(fit));
			}
			qs.update(q);
		}
		std::cout << " adds up to " << refinePoints.size() << " points. " << std::flush;
		uint32_t tmp = m_grid.tds().insert(refinePoints.begin(), refinePoints.end());
		refineCount += tmp;
		trWasRefined = refinePoints.size();
		std::cout << "Added " << tmp << " extra points. Quality: min=" << qs.min() << "max=" << qs.max() << std::endl;
		assignCellIds(grt, threadCount);
	}
	std::cout << "Refined triangulation with a total of " << refineCount << " extra points" << std::endl;
}

template<typename T_DUMMY>
void OsmTriangulationRegionStore::assignCellIds(OsmGridRegionTree<T_DUMMY> & grt, uint32_t threadCount) {
	struct Context {
		std::unordered_map<RegionList, uint32_t> cellListToCellId;
		OsmGridRegionTree<T_DUMMY> * grt;
		RegionListContainer * p_cellLists;
		sserialize::ProgressInfo pinfo;
		uint32_t finishedFaces;
		Triangulation::Finite_faces_iterator facesIt;
		Triangulation::Finite_faces_iterator facesEnd;
		std::mutex lock;
	} ctx;
	ctx.grt = &grt;
	ctx.p_cellLists = &m_cellLists;
	ctx.finishedFaces = 0;
	ctx.facesIt = m_grid.tds().finite_faces_begin();
	ctx.facesEnd = m_grid.tds().finite_faces_end();
	
	setInfinteFacesCellIds();
	//cells that are not in any region get cellid 0
	ctx.cellListToCellId[RegionList(ctx.p_cellLists, 0, 0)] = 0;
	for(uint32_t i(0), s(m_cellIdToCellList.size()); i < s; ++i) {
		ctx.cellListToCellId[m_cellIdToCellList.at(i)] = i;
	}

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
					if (fh->info().hasCellId()) {
						continue;
					}
					centroid = OsmTriangulationRegionStore::centroid(fh);
				}
				else {
					return;
				}
				lck.unlock();
				double x = CGAL::to_double(centroid.x());
				double y = CGAL::to_double(centroid.y());
				tmpCellList.clear();
				ctx->grt->test(x, y, tmpCellList);
				RegionList tmp(tmpCellList.begin(), tmpCellList.end());
				std::sort(tmp.begin(), tmp.end());
				uint32_t faceCellId = 0;
				lck.lock();
				if (!ctx->cellListToCellId.count(tmp)) {
					faceCellId = ctx->cellListToCellId.size();
					sserialize::MMVector<uint32_t>::size_type off = ctx->p_cellLists->size();
					ctx->p_cellLists->push_back(tmp.begin(), tmp.end());
					tmp = RegionList(ctx->p_cellLists, off, tmp.size());
					ctx->cellListToCellId[tmp] = faceCellId;
				}
				else {
					faceCellId = ctx->cellListToCellId.at(tmp);
				}
				assert((tmpCellList.size() || faceCellId == 0) && (faceCellId != 0 || !tmpCellList.size()));
				fh->info().setCellId(faceCellId);
				ctx->finishedFaces += 1;
				ctx->pinfo(ctx->finishedFaces);
			}
		}
	};
	ctx.pinfo.begin(m_grid.tds().number_of_faces(), "Setting initial cellids");
	std::vector<std::thread> threads;
	for(uint32_t i(0); i < threadCount; ++i) {
		threads.push_back(std::thread(WorkFunc(&ctx)));
	}
	for(std::thread & x : threads) {
		x.join();
	}
	ctx.pinfo.end();

	m_cellIdToCellList.resize(ctx.cellListToCellId.size());
	for(const auto & x : ctx.cellListToCellId) {
		m_cellIdToCellList.at(x.second) = x.first;
	}
}

template<typename TDummy, typename T_TRIANG_REFINER>
void OsmTriangulationRegionStore::init(OsmGridRegionTree<TDummy> & grt, uint32_t threadCount, T_TRIANG_REFINER * meshCriteria, RefinementAlgoTags refineAlgo) {
	if (!threadCount) {
		threadCount = std::thread::hardware_concurrency();
	}
	this->clear();
	{
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
				if ((itGp.first < -170.0 && prevGp.first > 170.0) || (itGp.first > 170.0 && prevGp.first < -170)) {
					std::cout << "Skipped edge crossing latitude boundary(-180->180)" << std::endl;
					continue;
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
#ifndef NDEBUG
		for(const std::pair<uint32_t, uint32_t> & s : segments) {
			assert(s.first < pts.size());
			assert(s.second < pts.size());
		}
#endif
		std::cout << "OsmTriangulationRegionStore: creating triangulation..." << std::flush;
		m_grid.tds().insert_constraints(pts.cbegin(), pts.cend(), segments.cbegin(), segments.cend());
		std::cout << "done" << std::endl;
	}
	
	//we now have to assign every face its cellid
	//a face lives in multiple regions, so every face has a unique list of region ids
	//faces that are connected and have the same region-list get the same cellid
	
	assignCellIds(grt, threadCount);
	
	//refine the triangulation
	if (meshCriteria) {
		switch (refineAlgo) {
		case MyRefineTag:
			myRefineMesh(*meshCriteria, grt, threadCount);
			break;
		case CGALRefineTag:
			cgalRefineMesh(*meshCriteria);
			break;
		default:
			break;
		}
	}
	
	m_refinedCellIdToUnrefined.reserve(m_cellIdToCellList.size());
	for(uint32_t i(0), s(m_cellIdToCellList.size()); i < s; ++i) {
		m_refinedCellIdToUnrefined.push_back(i);
	}
	std::cout << "Found " << m_cellIdToCellList.size() << " unrefined cells" << std::endl;
	assert(selfTest());
}


template<typename T_OUTPUT_ITERATOR>
void OsmTriangulationRegionStore::regionCells(uint32_t regionId, T_OUTPUT_ITERATOR out) {
	std::unordered_set<uint32_t> unrefinedMatching;
	for(uint32_t uCellId(0), s(m_cellIdToCellList.size()); uCellId < s; ++uCellId) {
		const RegionList & rl = m_cellIdToCellList[uCellId];
		if (std::find(rl.begin(), rl.end(), regionId) != rl.end()) {
			unrefinedMatching.insert(uCellId);
		}
	}
	for(uint32_t cellId(0), s(m_refinedCellIdToUnrefined.size()); cellId < s; ++cellId) {
		if (unrefinedMatching.count(m_refinedCellIdToUnrefined[cellId])) {
			*out = cellId;
			++out;
		}
	}
}

}//end namespace

#endif