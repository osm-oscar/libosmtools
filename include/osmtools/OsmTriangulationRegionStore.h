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

struct Construct_refine_points {
public:
	typedef enum {T_CENTROID=1, T_LONGEST_EDGE=2, T_ON_EDGES=3} Type;
public:
	Construct_refine_points(Type t, const sserialize::spatial::DistanceCalculator & dc) : m_t(t), m_dc(dc) {}
	template<typename T_TDS, typename T_OUTPUT_ITERATOR>
	void calc(const typename T_TDS::Face_handle fh, T_OUTPUT_ITERATOR out) {
		using ::osmtools::detail::OsmTriangulationRegionStore::centroid;
		using TDS = T_TDS;
		using Point = typename TDS::Point;
		switch (m_t) {
		case T_CENTROID:
			*out = CGAL::centroid(fh->vertex(0)->point(), fh->vertex(1)->point(), fh->vertex(2)->point());
			++out;
			break;
		case T_LONGEST_EDGE:
		{
			std::array<Point, 3> pts = {{
				fh->vertex(0)->point(),
				fh->vertex(1)->point(),
				fh->vertex(2)->point()
			}};
			std::array<double, 3> el;
			double longest = 0;
			for(int i(0); i < 3; ++i) {
				double p1x = CGAL::to_double(pts[i].x());
				double p1y = CGAL::to_double(pts[i].y());
				double p2x = CGAL::to_double(pts[TDS::cw(i)].x());
				double p2y = CGAL::to_double(pts[TDS::cw(i)].y());
				el[i] = m_dc.calc(p1x, p1y, p2x, p2y);
				longest = std::max(el[i], longest);
			}
			for(int i(0); i < 3; ++i) {
				if (el[i] == longest) {
					*out = CGAL::midpoint(fh->vertex(i)->point(), fh->vertex(TDS::cw(i))->point());
					++out;
				}
			}
			break;
		}
		case T_ON_EDGES:
			for(int i(0); i < 3; ++i) {
				*out = CGAL::midpoint(fh->vertex(i)->point(), fh->vertex(TDS::cw(i))->point());
				++out;
			}
			break;
		default:
			throw sserialize::InvalidEnumValueException("Construct_refine_points with type="  + std::to_string(m_t));
			break;
		};
	}
private:
	int m_t;
	sserialize::spatial::DistanceCalculator m_dc;
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
			//clip centroid distance to 1.0 (which is pretty small)
			q = std::max<double>(q, 1.0); 
			return q;
		}
	};
	
	using Construct_refine_points = ::osmtools::detail::OsmTriangulationRegionStore:: Construct_refine_points;
	
public:
	CentroidDistanceBaseMeshCriteria(Construct_refine_points::Type crpt = Construct_refine_points::T_CENTROID) : m_dc(sserialize::spatial::DistanceCalculator::DCT_GEODESIC_ACCURATE),
	m_crpt(crpt)
	{}
public:
	void crpt(Construct_refine_points::Type type) {
		m_crpt = type;
	}
	Construct_refine_points::Type crpt() const { return m_crpt; }
public:
	Construct_refine_points construct_refine_points_object() {
		return Construct_refine_points(m_crpt, dc());
	}
public:
	static bool usesCellIds() { return false; }
protected:
	const sserialize::spatial::DistanceCalculator & dc() const { return m_dc; }
protected:
	sserialize::spatial::DistanceCalculator m_dc;
	Construct_refine_points::Type m_crpt;
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
		
		inline Quality quality(CGAL::Mesh_2::Face_badness fb) const {
			if (fb == CGAL::Mesh_2::NOT_BAD) {
				return m_r;
			}
			else {
				return std::numeric_limits<double>::max();
			}
		}
		
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
class EdgeLengthRatioMeshCriteria: public CentroidDistanceBaseMeshCriteria<TDS> {
public:
	typedef CentroidDistanceBaseMeshCriteria<TDS> MyParentClass;
	typedef typename MyParentClass::Face_handle Face_handle;
	typedef typename MyParentClass::Point Point;
public:
	typedef double Quality;
	struct Is_bad {
		double m_r;
		sserialize::spatial::DistanceCalculator m_dc;
		Is_bad(double maxRatio, const sserialize::spatial::DistanceCalculator & dc) : m_r(maxRatio), m_dc(dc) {}
		
		CGAL::Mesh_2::Face_badness operator()(Quality q) const {
			if (q > m_r) {
				return CGAL::Mesh_2::IMPERATIVELY_BAD;
			}
			else {
				return CGAL::Mesh_2::NOT_BAD;
			}
		}
		
		Quality quality(CGAL::Mesh_2::Face_badness fb) const {
			if (fb == CGAL::Mesh_2::NOT_BAD) {
				return m_r;
			}
			else {
				return std::numeric_limits<double>::max();
			}
		}
		
		CGAL::Mesh_2::Face_badness operator()(Face_handle fh, Quality & q) const {
			std::array<Point, 3> pts = {{
				fh->vertex(0)->point(),
				fh->vertex(1)->point(),
				fh->vertex(2)->point()
			}};
			double longest = 0;
			double shortest = std::numeric_limits<double>::max();
			for(int i(0); i < 3; ++i) {
				double p1x = CGAL::to_double(pts[i].x());
				double p1y = CGAL::to_double(pts[i].y());
				double p2x = CGAL::to_double(pts[TDS::cw(i)].x());
				double p2y = CGAL::to_double(pts[TDS::cw(i)].y());
				double dist = m_dc.calc(p1x, p1y, p2x, p2y);
				longest = std::max(dist, longest);
				shortest = std::min(dist, shortest);
			}
			shortest = std::max<double>(shortest, std::numeric_limits<double>::epsilon());
			q = longest/shortest;
			return (*this)(q);
		};
	};
private:
	double m_r;
public:
	///@param maxRatio the maxium ratio between the shortest and longest edge
	EdgeLengthRatioMeshCriteria(double maxRatio) : m_r(maxRatio) {}
	Is_bad is_bad_object() const { return Is_bad(m_r, MyParentClass::dc()); }
};


template<typename TDS>
class LipschitzMeshCriteria: public CentroidDistanceBaseMeshCriteria<TDS> {
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
	static bool usesCellIds() { return true; }
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
	
	//Refines the triangulation if the distance between the centroid of the triangle and any of its defining points is larger than maxDist
	typedef detail::OsmTriangulationRegionStore::CentroidDistanceMeshCriteria<CDT> CentroidDistanceMeshCriteria;
	typedef detail::OsmTriangulationRegionStore::LipschitzMeshCriteria<CDT> LipschitzMeshCriteria;
	//If this type is selected, then MyRefineTag is mandatory
	typedef detail::OsmTriangulationRegionStore::RefineTrianglesWithCellIdMeshCriteria<LipschitzMeshCriteria> RegionOnlyLipschitzMeshCriteria;
	
	typedef enum {
		__CGALPredefinedRefineTag=0x10000,
		NoRefineTag=0x0,
		MyRefineTag=0x1,
		CGALRefineTag=0x2,
		CGALConformingTriangulationTag=0x4|__CGALPredefinedRefineTag,
		CGALGabrielTriangulationTag=0x8|__CGALPredefinedRefineTag,
		
	} RefinementAlgoTags;
	
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
		CS_CLEAN_GEOMETRY=0x4 //triangulation has clean geometry
		CS_HAVE_GRID=0x8, //grid for fast look-up is available
		CS_HAVE_REFINED_CELLS=0x10 // cells are refined
	} ConstructionState;
	
private:
	static Point centroid(const Face_handle & fh);
	//calculate hop-distances from rfh to all other faces of the cell of rfh and store their hop-distance in cellTriangMap and the cells triangs in cellTriangs
	void hopDistances(const Face_handle & rfh, std::vector<Face_handle> & cellTriangs, CGAL::Unique_hash_map<Face_handle, uint32_t> & cellTriangMap, uint32_t & maxHopDist);

	template<typename T_REFINER, typename T_Dummy>
	void myRefineMesh(T_REFINER & refiner, OsmGridRegionTree<T_Dummy> & grt, uint32_t threadCount);

	void setInfinteFacesCellIds();
	template<typename T_DUMMY>
	void assignCellIds(OsmGridRegionTree<T_DUMMY> & grt, uint32_t threadCount, bool reUseOld = true);
private:
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
	const GridLocator & grid() const { return m_grid; }
	const Triangulation & tds() const { return m_grid.tds(); }
	Triangulation & tds() { return m_grid.tds(); }
	inline uint32_t cellCount() const { return (uint32_t) m_refinedCellIdToUnrefined.size(); }
	inline uint32_t unrefinedCellCount() const { return (uint32_t) m_cellIdToCellList.size(); }
	///@param threadCount pass 0 for automatic deduction (uses std::thread::hardware_concurrency())
	///@param meshCriteria must be a modell of CGAL::MeshingCriteria_2 and provide function bool usesCellIds() if refineAlgo == MyRefineTag
	///@param refineAlgo selects the refinement algo, if CGALConformingTriangulationTag or CGALGabrielTriangulationTag, then meshCriteria is ignored
	template<typename TDummy,
		typename T_TRIANG_REFINER = OsmTriangulationRegionStore::RegionOnlyLipschitzMeshCriteria,
		typename T_REMOVED_EDGES = sserialize::Static::spatial::detail::Triangulation::PrintRemovedEdges>
	void init(OsmGridRegionTree< TDummy >& grt,
				uint32_t threadCount,
				T_TRIANG_REFINER* meshCriteria = 0,
				osmtools::OsmTriangulationRegionStore::RefinementAlgoTags refineAlgo = MyRefineTag,
				sserialize::Static::spatial::Triangulation::GeometryCleanType geoCleanType = sserialize::Static::spatial::Triangulation::GCT_NONE,
				T_REMOVED_EDGES re = T_REMOVED_EDGES());
	
	void initGrid(uint32_t gridLatCount, uint32_t gridLonCount);
	///Splits cells into connected cells
	void makeConnected();
	///First calls makeConnected(). Then splits cells into smaller cells as long as refiner returns true
	///This is done in multiple runs where each cell is split into up to numVoronoiSplitRuns smaller cells.
	///Cells are not split into equally sized cells but rather by their voronoi diagram
	///refine cells by connectedness so that all cells form a connected polygon (with holes)
	void refineCells(std::shared_ptr<CellRefinerInterface> refiner, uint32_t runs, uint32_t splitPerRun, uint32_t threadCount);
	
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

template<typename T_REFINER, typename T_Dummy>
void OsmTriangulationRegionStore::myRefineMesh(T_REFINER & refiner, OsmGridRegionTree<T_Dummy> & grt, uint32_t threadCount) {
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
	
	if (T_REFINER::usesCellIds()) {
		assignCellIds(grt, threadCount);
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

///assign cellIds, while keeping the old cellIds if reUseOld is true
template<typename T_DUMMY>
void OsmTriangulationRegionStore::assignCellIds(OsmGridRegionTree<T_DUMMY> & grt, uint32_t threadCount, bool reUseOld) {
	if (!reUseOld) {
		m_cellIdToCellList.clear();
		m_cellLists.clear();
		m_refinedCellIdToUnrefined.clear();
		m_isConnected = false;
		for(All_faces_iterator it(m_grid.tds().all_faces_begin()), end(m_grid.tds().all_faces_end()); it != end; ++it) {
			it->info().clear();
		}
	}
	struct CellListKey {
		uint64_t hash;
		RegionList list;
		CellListKey() : hash(0) {}
		explicit CellListKey(uint64_t h, const RegionList & l) : hash(h), list(l) {}
		explicit CellListKey(uint64_t h, RegionListContainer* c, RegionList::size_type off, RegionList::size_type size) : hash(h), list(c, off, size) {}
		CellListKey(const CellListKey &) = default;
		CellListKey(CellListKey &&) = default;
		CellListKey & operator=(const CellListKey&) = default;
		CellListKey & operator=(CellListKey&&) = default;
		inline bool operator!=(const CellListKey & other) const { return hash != other.hash || list != other.list; }
		inline bool operator==(const CellListKey & other) const  { return hash == other.hash && list == other.list; }
	};
	
	struct CellListKeyHasher {
		CellListKeyHasher() = default;
		CellListKeyHasher(const CellListKeyHasher &) = default;
		CellListKeyHasher & operator=(const CellListKeyHasher &) = default;
		inline std::size_t operator()(const CellListKey & v) const { return v.hash; }
	};
	
	struct Context {
		std::unordered_map<CellListKey, uint32_t, CellListKeyHasher> cellListToCellId;
		OsmGridRegionTree<T_DUMMY> * grt;
		RegionListContainer * p_cellLists;
		sserialize::ProgressInfo pinfo;
		uint32_t finishedFaces;
		Triangulation::Finite_faces_iterator facesIt;
		Triangulation::Finite_faces_iterator facesEnd;
		std::mutex iteratorLock;
		std::mutex cellListLock;
	} ctx;
	ctx.grt = &grt;
	ctx.p_cellLists = &m_cellLists;
	ctx.finishedFaces = 0;
	ctx.facesIt = m_grid.tds().finite_faces_begin();
	ctx.facesEnd = m_grid.tds().finite_faces_end();
	
	setInfinteFacesCellIds();
	//cells that are not in any region get cellid 0
	{
		std::hash<RegionList> hasher;
		{
			RegionList tmp(ctx.p_cellLists, 0, 0);
			ctx.cellListToCellId[CellListKey(hasher(tmp), tmp)] = 0;
		}
		for(uint32_t i(0), s((uint32_t) m_cellIdToCellList.size()); i < s; ++i) {
			const auto & tmp = m_cellIdToCellList.at(i);
			ctx.cellListToCellId[CellListKey(hasher(tmp), tmp)] = i;
		}
	}

	struct WorkFunc {
		Context * ctx;
		RegionList::container_type tmpCellListContainer;
		std::back_insert_iterator<RegionList::container_type> tmpCellListInserter;
		std::hash<RegionList> cellListHasher;
		CellListKey tmpCellListKey;
		Point centroid;
		Face_handle fh;
		WorkFunc(Context * ctx) : ctx(ctx), tmpCellListInserter(tmpCellListContainer) {
			tmpCellListContainer.reserve(64);
		}
		WorkFunc(const WorkFunc & other) : WorkFunc(other.ctx) {}
		void operator()() {
			while (true) {
				{
					std::lock_guard<std::mutex> lck(ctx->iteratorLock);
					while (true) {
						if (ctx->facesIt == ctx->facesEnd) {
							return;
						}
						if (!ctx->facesIt->info().hasCellId()) {
							break;
						}
						++(ctx->facesIt);
					}
					fh = ctx->facesIt;
					++(ctx->facesIt);
					centroid = OsmTriangulationRegionStore::centroid(fh);
				}
				
				uint32_t faceCellId = 0;
				{
					double x = CGAL::to_double(centroid.x());
					double y = CGAL::to_double(centroid.y());
					tmpCellListContainer.clear();
					ctx->grt->find(x, y, tmpCellListInserter);
					std::sort(tmpCellListContainer.begin(), tmpCellListContainer.end());
					SSERIALIZE_NORMAL_ASSERT(sserialize::is_strong_monotone_ascending(tmpCellListContainer.begin(), tmpCellListContainer.end()));
					tmpCellListKey.list = RegionList(&tmpCellListContainer);
					tmpCellListKey.hash = cellListHasher(tmpCellListKey.list);
				}
				
				{
					std::lock_guard<std::mutex> lck(ctx->cellListLock);
					auto cellListToCellIdIt = ctx->cellListToCellId.find(tmpCellListKey);
					if (cellListToCellIdIt == ctx->cellListToCellId.end()) {
						faceCellId = (uint32_t) ctx->cellListToCellId.size();
						auto off = ctx->p_cellLists->size();
						ctx->p_cellLists->push_back(tmpCellListKey.list.begin(), tmpCellListKey.list.end());
						ctx->cellListToCellId[CellListKey(tmpCellListKey.hash, ctx->p_cellLists, off, tmpCellListKey.list.size())] = faceCellId;
						SSERIALIZE_CHEAP_ASSERT_EQUAL(ctx->cellListToCellId.size(), faceCellId+1);
					}
					else {
						faceCellId = cellListToCellIdIt->second;
					}
				}
				
				SSERIALIZE_CHEAP_ASSERT((tmpCellListKey.list.size() || faceCellId == 0) && (faceCellId != 0 || !tmpCellListKey.list.size()));
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
		m_cellIdToCellList.at(x.second) = x.first.list;
	}
	
	SSERIALIZE_EXPENSIVE_ASSERT(selfTest());
}

template<typename TDummy, typename T_TRIANG_REFINER, typename T_REMOVED_EDGES>
void
OsmTriangulationRegionStore::init(
	OsmGridRegionTree<TDummy> & grt,
	uint32_t threadCount,
	T_TRIANG_REFINER * meshCriteria, RefinementAlgoTags refineAlgo,
	sserialize::Static::spatial::Triangulation::GeometryCleanType geoCleanType, T_REMOVED_EDGES re
	)
{
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
		
		auto handlePolygonPoints = [&gpToId](const GeoPolygon * gp) {
			typename GeoPolygon::const_iterator it(gp->cbegin()), end(gp->cend());
			for(; it != end; ++it) {
				RawGeoPoint itGp = *it;
				if (!gpToId.count(itGp)) {
					uint32_t gpId = (uint32_t) gpToId.size();
					gpToId[itGp] = gpId;
				}
			}
		};
		
		std::cout << "OsmTriangulationRegionStore: extracting points..." << std::flush;
		for(sserialize::spatial::GeoRegion* r : grt.regions()) {
			if (r->type() == sserialize::spatial::GS_POLYGON) {
				const GeoPolygon * gp = static_cast<const GeoPolygon*>(r);
				handlePolygonPoints(gp);
			}
			else if (r->type() == sserialize::spatial::GS_MULTI_POLYGON) {
				const GeoMultiPolygon * gmp = static_cast<const GeoMultiPolygon*>(r);
				for(const GeoPolygon & gp : gmp->outerPolygons()) {
					handlePolygonPoints(&gp);
				}
				for(const GeoPolygon & gp : gmp->innerPolygons()) {
					handlePolygonPoints(&gp);
				}
			}
		}
		std::cout << "done" << std::endl;
		
		auto handlePolygonSegments = [&gpToId, &segments](const GeoPolygon * gp) {
			typename GeoPolygon::const_iterator it(gp->cbegin()), prev(gp->cbegin()), end(gp->cend());
			for(++it; it != end; ++it, ++prev) {
				RawGeoPoint itGp = *it;
				RawGeoPoint prevGp = *prev;
				if ((itGp.first < -170.0 && prevGp.first > 170.0) || (itGp.first > 170.0 && prevGp.first < -170)) {
					std::cout << "Skipped edge crossing latitude boundary(-180->180)\n";
					continue;
				}
				segments.insert( std::pair<uint32_t, uint32_t>(gpToId.at(itGp), gpToId.at(prevGp)) );
			}
		};
		std::cout << "OsmTriangulationRegionStore: extracting segments..." << std::flush;
		for(sserialize::spatial::GeoRegion* r : grt.regions()) {
			if (r->type() == sserialize::spatial::GS_POLYGON) {
				const GeoPolygon * gp = static_cast<const GeoPolygon*>(r);
				handlePolygonSegments(gp);
			}
			else if (r->type() == sserialize::spatial::GS_MULTI_POLYGON) {
				const GeoMultiPolygon * gmp = static_cast<const GeoMultiPolygon*>(r);
				for(const GeoPolygon & gp : gmp->outerPolygons()) {
					handlePolygonSegments(&gp);
				}
				for(const GeoPolygon & gp : gmp->innerPolygons()) {
					handlePolygonSegments(&gp);
				}
			}
		}
		std::cout << "done" << std::endl;
		
		std::cout << "Found " << gpToId.size() << " different points creating " << segments.size() << " different segments" << std::endl;
		
		std::cout << "Converting points to CGAL points..." << std::flush;
		pts.resize(gpToId.size());
		for(const auto & x : gpToId) {
			pts.at(x.second) = Point(x.first.first, x.first.second);
		}
		gpToId = decltype(gpToId)();
		std::cout << "done" << std::endl;
		
#ifdef SSERIALIZE_EXPENSIVE_ASSERT_ENABLED
		for(const std::pair<uint32_t, uint32_t> & s : segments) {
			SSERIALIZE_EXPENSIVE_ASSERT(s.first < pts.size());
			SSERIALIZE_EXPENSIVE_ASSERT(s.second < pts.size());
		}
#endif
		std::cout << "OsmTriangulationRegionStore: creating triangulation..." << std::flush;
		sserialize::TimeMeasurer tm;
		tm.begin();
		m_grid.tds().insert_constraints(pts.cbegin(), pts.cend(), segments.cbegin(), segments.cend());
		tm.end();
		std::cout << "took " << tm << std::endl;
	}
	
	
	//refine the triangulation
	if (refineAlgo == MyRefineTag && meshCriteria) {
		myRefineMesh(*meshCriteria, grt, threadCount);
	}
	else if (refineAlgo == CGALRefineTag) {
		CGAL::Delaunay_mesher_2<Triangulation, T_TRIANG_REFINER> mesher(m_grid.tds(), *meshCriteria);
		mesher.refine_mesh();
	}
	else if (refineAlgo & __CGALPredefinedRefineTag) {
		CGAL::Triangulation_conformer_2<Triangulation> conform(m_grid.tds());
		switch (refineAlgo) {
		case CGALConformingTriangulationTag:
			conform.make_conforming_Delaunay();
			break;
		case CGALGabrielTriangulationTag:
			conform.make_conforming_Delaunay();
			conform.make_conforming_Gabriel();
			break;
		default:
			throw sserialize::UnsupportedFeatureException("Unsupported refinement algorithm: " + std::to_string(refineAlgo));
			break;
		}
	}

	//we now have to assign every face its cellid
	//a face lives in multiple regions, so every face has a unique list of region ids
	//faces that are connected and have the same region-list get the same cellid
	assignCellIds(grt, threadCount);

	
	if (geoCleanType != sserialize::Static::spatial::Triangulation::GCT_NONE) {
		sserialize::Static::spatial::Triangulation::prepare(tds(), re,  geoCleanType, 0.01);
		assignCellIds(grt, threadCount, false);
	}
	
	m_refinedCellIdToUnrefined.reserve(m_cellIdToCellList.size());
	for(uint32_t i(0), s((uint32_t) m_cellIdToCellList.size()); i < s; ++i) {
		m_refinedCellIdToUnrefined.push_back(i);
	}
	std::cout << "Found " << m_cellIdToCellList.size() << " unrefined cells" << std::endl;
	SSERIALIZE_VERY_EXPENSIVE_ASSERT(selfTest());
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
