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
#include <CGAL/Triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

//for ExtenedInt64q kernel
#include <libratss/CGAL/ExtendedInt64Cartesian.h>

#include <assert.h>
#include <thread>
#include <mutex>

namespace osmtools {
class OsmTriangulationRegionStore;

namespace detail {
namespace OsmTriangulationRegionStore {

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
	inline uint32_t size() const { return (uint32_t) m_nodes.size(); }
	inline uint32_t cellId() const { return m_cellId; }
	///@param maxHopDistRoot to the root node from where the bfs-tree is deepest.
	///       This is only true if calcMaxHopDistance returned true,
	///       otherwise it is an approximation
	bool calcMaxHopDistance(uint32_t& maxHopDistRoot);
	uint32_t calcDiameter(uint32_t * startNode, uint32_t * endNode);
	const FaceNode & node(uint32_t pos) const { return m_nodes.at(pos); }
	FaceNode & node(uint32_t pos) { return m_nodes.at(pos); }
};

template<typename TDS>
class CTGraph: public detail::OsmTriangulationRegionStore::CTGraphBase {
public:
	typedef typename TDS::Face_handle Face_handle;
private:
	friend class osmtools::OsmTriangulationRegionStore;
private:
	std::vector<Face_handle> m_faces;
	CGAL::Unique_hash_map<Face_handle, uint32_t> m_faceToNodeId;
public:
	CTGraph() {}
	virtual ~CTGraph() {}
	Face_handle face(uint32_t faceNodeId) const { return m_faces.at(faceNodeId); }
	uint32_t node(const Face_handle & fh) {
		if (m_faceToNodeId.is_defined(fh)) {
			return m_faceToNodeId[fh];
		}
		throw std::out_of_range("OsmTriangulationRegionStore::CellGraph::node");
		return CTGraph::NullFace;
	}
	using CTGraphBase::node;
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

	inline uint32_t size() const { return (uint32_t) m_nodes.size(); }
	inline const CellNode & node(uint32_t pos) const { return m_nodes.at(pos); }
	inline CellNode & node(uint32_t pos) { return m_nodes.at(pos); }
	
	sserialize::UByteArrayAdapter & append(sserialize::UByteArrayAdapter & dest, const std::unordered_map<uint32_t, uint32_t> & myIdsToGhCellIds) const;
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
		void setCellId(uint32_t cellId);
		uint32_t cellId() const;
		bool hasCellId() const;
	};

	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
// 	typedef CGAL::Exact_predicates_exact_constructions_kernel K;
// 	typedef CGAL::Filtered_simple_cartesian_extended_integer_kernel K;
	static constexpr bool KernelHasThreadSafeNumberType = std::is_same<K, CGAL::Filtered_simple_cartesian_extended_integer_kernel>::value;
// 	typedef CGAL::Simple_cartesian_extended_integer_kernel K;
	typedef CGAL::Exact_intersections_tag Itag;
	typedef CGAL::Triangulation_vertex_base_2<K> Vb;
	typedef CGAL::Delaunay_mesh_vertex_base_2<K, Vb> MVb;
	typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo, K> FbWithInfo;
	typedef CGAL::Constrained_triangulation_face_base_2<K, FbWithInfo> CFb;
	typedef CGAL::Delaunay_mesh_face_base_2<K, CFb> MFb;
	typedef CGAL::Triangulation_data_structure_2<MVb, MFb> Tds;
	typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds, Itag> CDTBase;
	typedef CGAL::Constrained_triangulation_plus_2<CDTBase> CDTP;

	//choose here to either use the normal CDT or the CDTPlus
	//WARNING: CDTPlus results in a double free or corruption on planet!
	typedef CDTBase CDT;
	typedef CDT Triangulation;
	
	typedef Triangulation::Point Point;
	typedef osmtools::GridLocator<Triangulation, KernelHasThreadSafeNumberType> GridLocator;
	typedef GridLocator::Vertex_handle Vertex_handle;
	typedef GridLocator::Face_handle Face_handle;
	typedef sserialize::MMVector<uint32_t> RegionListContainer;
	typedef sserialize::CFLArray<RegionListContainer> RegionList;
	typedef GridLocator::TriangulationDataStructure::Finite_faces_iterator Finite_faces_iterator;
	typedef GridLocator::TriangulationDataStructure::All_faces_iterator All_faces_iterator;
	typedef sserialize::SimpleBitVector SimpleBitVector;
public:
	typedef detail::OsmTriangulationRegionStore::CTGraph<Tds> CTGraph;
	typedef detail::OsmTriangulationRegionStore::CellGraph CellGraph;

	typedef enum {
		TRAS_NoRefinement=0x0,
		TRAS_Osmtools=0x1,
		TRAS_DelaunayMesher=0x2,
		TRAS_ConformingTriangulation=0x4,
		TRAS_GabrielTriangulation=0x8
	} TriangulationRefinementAlgorithmSelector;
	
	class CellRefinerInterface {
	public:
		struct State {
			///Triangle count of all cells (includes new cells)
			std::vector<uint32_t> cellSizes;
			///A face representative all cells (includes new cells)
			std::vector<Face_handle> cellRep;
			///the cell triangle graph of the current source cell that is about to be split
			CTGraph cg;
			///the list of new cell Ids, maps from CTGraph::NodeId to CellId
			std::vector<uint32_t> newFaceCellIds;
			//the list of new cell representatives, maps to CTGraph::NodeId
			std::vector<uint32_t> newCellReps;
			///the set of new cell ids
			std::unordered_set<uint32_t> currentCells;
		};
	public:
		virtual ~CellRefinerInterface() {}
		///@return indicate if refinement is necessary at all
		virtual bool init(const ::osmtools::OsmTriangulationRegionStore &) { return false; }
		//tells the refiner that refinement begins
		virtual void begin() {}
		//tells the refiner that refinement ends
		virtual void end() {}
		///@return true if source cell needs further refinement
		virtual bool refine(const State &);
		///@return true if cell needs further refinement
		virtual bool refine(uint32_t /*cellId*/, const State &) { return false; }
		virtual CellRefinerInterface * copy() = 0;
	};
	
public:
	static uint32_t InfiniteFacesCellId;
	static uint32_t UnsetFacesCellId;
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
	
	typedef enum {
		CS_EMPTY=0x0,
		CS_HAVE_TRIANGULATION=0x1, //triangulation is there, but no cells are assigned
		CS_HAVE_CELLS=0x2, //faces have cell ids assigned
		CS_HAVE_REFINED_TRIANGULATION=0x4,
		CS_HAVE_SNAPPED_TRIANGULATION=0x8, //triangulation has snapped geometry
		CS_HAVE_GRID=0x10, //grid for fast look-up is available
		CS_HAVE_REFINED_CELLS=0x20 // cells are refined
	} ConstructionState;
	
private:
	static Point centroid(const Face_handle & fh);
	//calculate hop-distances from rfh to all other faces of the cell of rfh and store their hop-distance in cellTriangMap and the cells triangs in cellTriangs
	void hopDistances(const Face_handle & rfh, std::vector<Face_handle> & cellTriangs, CGAL::Unique_hash_map<Face_handle, uint32_t> & cellTriangMap, uint32_t & maxHopDist);

	template<typename T_TRIAN_REFINER>
	void myRefineMesh(T_TRIAN_REFINER refiner, uint32_t threadCount);

	void setInfinteFacesCellIds();
	
	void refineTriangulationFinalize();
	
	std::shared_ptr<OsmGridRegionTreeBase> & grt() { return m_grt; }
private:
	std::shared_ptr<OsmGridRegionTreeBase> m_grt;
	GridLocator m_grid;
	RegionListContainer m_cellLists;
	std::vector<RegionList> m_cellIdToCellList;
	std::vector<uint32_t> m_refinedCellIdToUnrefined;
	bool m_isConnected;
	std::mutex m_lock;
	int m_cs; //construction state
public:
	OsmTriangulationRegionStore();
	OsmTriangulationRegionStore(const OsmTriangulationRegionStore & other) = delete;
	~OsmTriangulationRegionStore() {}
	void clear();
	void clearCells();
	void clearRefinement();
	const std::shared_ptr<OsmGridRegionTreeBase> & grt() const { return m_grt; }
	const GridLocator & grid() const { return m_grid; }
	const Triangulation & tds() const { return m_grid.tds(); }
	Triangulation & tds() { return m_grid.tds(); }
	inline uint32_t cellCount() const { return (uint32_t) m_refinedCellIdToUnrefined.size(); }
	inline uint32_t unrefinedCellCount() const { return (uint32_t) m_cellIdToCellList.size(); }
	///@param threadCount pass 0 for automatic deduction (uses std::thread::hardware_concurrency())
	void init(std::shared_ptr<OsmGridRegionTreeBase> grt, uint32_t threadCount);

	///@param meshCriteria must be a modell of CGAL::MeshingCriteria_2 and provide function bool usesCellIds() if refineAlgo == MyRefineTag
	///@param refineAlgo selects the refinement algo
	template<typename T_TRIANG_REFINER>
	void refineTriangulation(const T_TRIANG_REFINER & meshCriteria, TriangulationRefinementAlgorithmSelector refineAlgo, uint32_t threadCount);
	///used for CGALConformingTriangulationTag or CGALGabrielTriangulationTag
	void refineTriangulation(TriangulationRefinementAlgorithmSelector refineAlgo);
	
	template<typename T_REMOVED_EDGES = sserialize::Static::spatial::detail::Triangulation::PrintRemovedEdges>
	void snapTriangulation(sserialize::Static::spatial::Triangulation::GeometryCleanType geoCleanType, T_REMOVED_EDGES re = T_REMOVED_EDGES());
	
	void initGrid(uint32_t gridLatCount, uint32_t gridLonCount);
	
	void assignCellIds(uint32_t threadCount);
	
	///Splits cells into connected cells
	void makeConnected();
	///First calls makeConnected(). Then splits cells into smaller cells as long as refiner returns true
	///This is done in multiple runs where each cell is split into up to numVoronoiSplitRuns smaller cells.
	///Cells are not split into equally sized cells but rather by their voronoi diagram
	///refine cells by connectedness so that all cells form a connected polygon (with holes)
	void refineCells(std::shared_ptr<CellRefinerInterface> refiner, uint32_t runs, uint32_t splitPerRun, uint32_t threadCount);
	
	void refineBySize(uint32_t cellSizeTh, uint32_t runs, uint32_t splitPerRun, uint32_t threadCount);
	
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
	std::vector<sserialize::spatial::GeoPoint> cellCenterOfMass(const std::unordered_map<uint32_t, uint32_t> & myIdsToGhCellIds);

	template<typename T_OUTPUT_ITERATOR>
	void regionCells(uint32_t regionId, T_OUTPUT_ITERATOR out);
	
	void printStats(std::ostream & out);
	
	bool selfTest();
	///serializes to sserialize::Static::spatial::TriangulationRegionArrangement
	sserialize::UByteArrayAdapter & append(sserialize::UByteArrayAdapter& dest,
											sserialize::ItemIndexFactory& idxFactory,
											sserialize::Static::spatial::Triangulation::GeometryCleanType gct);
	
	///serializes to sserialize::Static::spatial::TriangulationGeoHierarchyArrangement
	sserialize::UByteArrayAdapter & append(sserialize::UByteArrayAdapter& dest,
											const std::unordered_map< uint32_t, uint32_t >& myIdsToGhCellIds,
											sserialize::Static::spatial::Triangulation::GeometryCleanType gct);
	
	bool equal(const sserialize::Static::spatial::TriangulationGeoHierarchyArrangement & ra, const std::unordered_map<uint32_t, uint32_t> & myIdsToGhCellIds);
};

namespace detail {
namespace OsmTriangulationRegionStore {
	
class RefineByTriangleCount: public ::osmtools::OsmTriangulationRegionStore::CellRefinerInterface {
public:
	RefineByTriangleCount(uint32_t cellSizeTh);
	virtual ~RefineByTriangleCount() {}
public:
	virtual bool init(const ::osmtools::OsmTriangulationRegionStore & store) override;
	virtual void begin() override;
	virtual void end() override;
	virtual bool refine(uint32_t cellId, const State & state) override;
	virtual CellRefinerInterface * copy() override;
private:
	uint32_t m_cellSizeTh;
};

class RefineBySize: public ::osmtools::OsmTriangulationRegionStore::CellRefinerInterface {
public:
	RefineBySize(double maxCellDiameter);
	virtual ~RefineBySize() {}
public:
	virtual bool init(const ::osmtools::OsmTriangulationRegionStore & store) override;
	virtual void begin() override;
	virtual void end() override;
	virtual bool refine(const State & state) override;
	virtual bool refine(uint32_t cellId, const State & state) override;
	virtual CellRefinerInterface * copy() override;
private:
	double m_maxCellDiameter;
	sserialize::spatial::DistanceCalculator m_dc;
};

}}//end namespace detail::OsmTriangulationRegionStore



template<typename T_TRIANG_REFINER>
void OsmTriangulationRegionStore::refineTriangulation(const T_TRIANG_REFINER & meshCriteria, TriangulationRefinementAlgorithmSelector refineAlgo, uint32_t threadCount) {
	if (refineAlgo == TRAS_NoRefinement || !(m_cs & CS_HAVE_TRIANGULATION)) {
		return;
	}
	
	//refine the triangulation
	switch (refineAlgo) {
	case TRAS_Osmtools:
		myRefineMesh(meshCriteria, threadCount);
		m_cs |= CS_HAVE_REFINED_TRIANGULATION;
		refineTriangulationFinalize();
		break;
	case TRAS_DelaunayMesher:
	{
		CGAL::Delaunay_mesher_2<Triangulation, T_TRIANG_REFINER> mesher(m_grid.tds(), meshCriteria);
		mesher.refine_mesh();
		m_cs |= CS_HAVE_REFINED_TRIANGULATION;
		refineTriangulationFinalize();
	}
	case TRAS_ConformingTriangulation:
	case TRAS_GabrielTriangulation:
	default:
		refineTriangulation(refineAlgo);
		break;
	}
}

template<typename T_REMOVED_EDGES>
void
OsmTriangulationRegionStore::snapTriangulation(sserialize::Static::spatial::Triangulation::GeometryCleanType geoCleanType, T_REMOVED_EDGES re)
{
	if (geoCleanType != sserialize::Static::spatial::Triangulation::GCT_NONE) {
		sserialize::Static::spatial::Triangulation::prepare(tds(), re,  geoCleanType, 0.01);
		if (m_cs & CS_HAVE_REFINED_CELLS) {
			std::cerr << "WARNING: OsmTriangulationRegionStore::snapTriangulation: removing cell refinement" << std::endl;
			clearRefinement();
		}
		if (m_cs & CS_HAVE_CELLS) {
			std::cerr << "WARNING: OsmTriangulationRegionStore::snapTriangulation: removing cells" << std::endl;
			m_cs &= ~CS_HAVE_CELLS;
		}
		m_cs |= CS_HAVE_SNAPPED_TRIANGULATION;
	}
	SSERIALIZE_VERY_EXPENSIVE_ASSERT(selfTest());
}

template<typename T_REFINER>
void OsmTriangulationRegionStore::myRefineMesh(T_REFINER refiner, uint32_t threadCount) {
	struct RefinePoint {
		Vertex_handle vh;
		Point p;
		RefinePoint(const Vertex_handle & vh, const Point & p) : vh(vh), p(p) {}
	};
	struct MyBackInserter {
		MyBackInserter(std::vector<RefinePoint> * dest) : dest(dest) {}
		MyBackInserter & operator=(const Point & p) {
			dest->emplace_back(vh, p);
			return *this;
		}
		MyBackInserter & operator*() { return *this; }
		MyBackInserter & operator++() { return *this; }
		std::vector<RefinePoint> * dest;
		Vertex_handle vh;
	};
	
	if (T_REFINER::usesCellIds() && !(m_cs & CS_HAVE_CELLS)) {
		assignCellIds(threadCount);
	}
	
	uint32_t refineCount = 0;
	std::vector<RefinePoint> refinePoints;
	sserialize::MinMax<typename T_REFINER::Quality> qs;
	bool trWasRefined = true;
	typename T_REFINER::Is_bad bo(refiner.is_bad_object());
	typename T_REFINER::Construct_refine_points crp(refiner.construct_refine_points_object());
	typename T_REFINER::Quality q = 0;
	typename T_REFINER::Face_handle rfh;
	MyBackInserter refinePointsInserter(&refinePoints);
	for(uint32_t round(0); round < 10000 && trWasRefined; ++round) {
		std::cout << "Trianglerefinement round " << round << std::flush;
		trWasRefined = false;
		qs.reset();
		refinePoints.clear();
		for(Finite_faces_iterator fit(finite_faces_begin()), fend(finite_faces_end()); fit != fend; ++fit) {
			rfh = fit;
			if (CGAL::Mesh_2::NOT_BAD != bo.operator()(rfh, q)) {
				crp.template calc<Triangulation, MyBackInserter>(rfh, refinePointsInserter);
			}
			qs.update(q);
		}
		std::cout << " adds up to " << refinePoints.size() << " points. " << std::flush;
		for(const RefinePoint & rp : refinePoints) {
			m_grid.tds().insert(rp.p, rp.vh->face());
			++refineCount;
		}
		trWasRefined = refinePoints.size();
		std::cout << "Added " << refinePoints.size() << " extra points. Quality: min=" << qs.min() << ", max=" << qs.max() << std::endl;
	}
	std::cout << "Refined triangulation with a total of " << refineCount << " extra points" << std::endl;
}

template<typename T_OUTPUT_ITERATOR>
void OsmTriangulationRegionStore::regionCells(uint32_t regionId, T_OUTPUT_ITERATOR out) {
	std::unordered_set<uint32_t> unrefinedMatching;
	for(uint32_t uCellId(0), s((uint32_t) m_cellIdToCellList.size()); uCellId < s; ++uCellId) {
		const RegionList & rl = m_cellIdToCellList[uCellId];
		if (std::find(rl.begin(), rl.end(), regionId) != rl.end()) {
			unrefinedMatching.insert(uCellId);
		}
	}
	for(uint32_t cellId(0), s((uint32_t) m_refinedCellIdToUnrefined.size()); cellId < s; ++cellId) {
		if (unrefinedMatching.count(m_refinedCellIdToUnrefined[cellId])) {
			*out = cellId;
			++out;
		}
	}
}

}//end namespace

#endif
