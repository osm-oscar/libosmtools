#ifndef LIBOSMTOOLS_OSM_TRIANGULATION_REGION_STORE_H
#define LIBOSMTOOLS_OSM_TRIANGULATION_REGION_STORE_H
#include <unordered_map>

#include <sserialize/utility/hashspecializations.h>

#include <osmtools/OsmGridRegionTree.h>
#include <osmtools/TriangulationGridLocater.h>

#include <sserialize/utility/TimeMeasuerer.h>

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
class OsmTriangulationRegionStore;

namespace detail {
namespace OsmTriangulationRegionStore {

class CellGraphBase {
private:
	friend class osmtools::OsmTriangulationRegionStore;
public:
	static const uint32_t NullFace = 0xFFFFFFFF;
	struct FaceNode {
		static const uint32_t NullNeighbor = CellGraphBase::NullFace;
		///by definition: if neighbours[i] == 0xFFFFFFFF => no neighbor
		uint32_t neighbours[3];
	};
	struct SimpleBitVector {
		std::vector<uint64_t> m_d;
		inline void resize(uint32_t count) { m_d.resize(count/64+1, 0); }
		inline void set(uint32_t pos) { m_d.at(pos/64) |= (static_cast<uint64_t>(1) << (pos%64)); }
		inline bool isSet(uint32_t pos) { return m_d.at(pos/64) & (static_cast<uint64_t>(1) << (pos%64)); }
		inline void reset() { m_d.assign(m_d.size(), 0); }
	};
	typedef std::vector<FaceNode> NodesContainer;
private:
	uint32_t m_cellId;
	NodesContainer m_nodes;
protected:
	NodesContainer & nodes() { return m_nodes; }
	const NodesContainer & nodes() const { return m_nodes; }
public:
	CellGraphBase() {}
	virtual ~CellGraphBase() {}
	inline uint32_t size() const { return m_nodes.size(); }
	inline uint32_t cellId() const { return m_cellId; }
	///@param bfsTree (nodeId, hopdistance from faceNodeId)
	void calcMaxHopDistance(std::vector< std::pair< uint32_t, uint32_t > >& bfsTree);
	const FaceNode & node(uint32_t pos) const { return m_nodes.at(pos); }
	FaceNode & node(uint32_t pos) { return m_nodes.at(pos); }
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
	typedef GridLocator::TriangulationDataStructure::Finite_faces_iterator Finite_faces_iterator;
public:
	class CellGraph: public detail::OsmTriangulationRegionStore::CellGraphBase {
	private:
		friend class OsmTriangulationRegionStore;
	private:
		std::vector<Face_handle> m_faces;
		CGAL::Unique_hash_map<Face_handle, uint32_t> m_faceToNodeId;
	public:
		CellGraph() {}
		virtual ~CellGraph() {}
		Face_handle face(uint32_t faceNodeId);
		uint32_t node(const Face_handle & fh);
		using CellGraphBase::node;
	};
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
private:
	GridLocator m_grid;
	FaceCellIdMap m_faceToCellId;
	RegionListContainer m_cellLists;
	std::vector<RegionList> m_cellIdToCellList;
	std::vector<uint32_t> m_refinedCellIdToUnrefined;
	std::mutex m_lock;
public:
	OsmTriangulationRegionStore() {}
	OsmTriangulationRegionStore(const OsmTriangulationRegionStore & other) = delete;
	~OsmTriangulationRegionStore() {}
	void clear();
	uint32_t cellCount() const { return m_refinedCellIdToUnrefined.size(); }
	///@param cellSizeTh threshold over which a cell is split into smaller cells, pass std::numeric_limits<uint32_t>::max() to disable
	///@param threadCount pass 0 for automatic deduction (uses std::thread::hardware_concurrency())
	template<typename TDummy>
	void init(OsmGridRegionTree<TDummy> & grt, uint32_t gridLatCount, uint32_t gridLonCount, uint32_t cellSizeTh, uint32_t threadCount);
	inline uint32_t cellId(double lat, double lon);
	inline uint32_t cellId(const sserialize::spatial::GeoPoint & gp) { return cellId(gp.lat(), gp.lon()); }
	///@thread-safety no
	uint32_t cellId(const Face_handle & fh);
	inline const RegionListContainer & regionLists() const { return m_cellLists; }
	inline RegionListContainer & regionLists() { return m_cellLists; }
	inline const RegionList & regions(uint32_t cellId);
	Finite_faces_iterator finite_faces_begin() { return m_grid.tds().finite_faces_begin(); }
	Finite_faces_iterator finite_faces_end() { return m_grid.tds().finite_faces_end(); }
	
	void cellGraph(const Face_handle& rfh, CellGraph& cg);
	std::vector<Face_handle> cellRepresentatives();
	
	void printStats(std::ostream & out);
};


template<typename TDummy>
void OsmTriangulationRegionStore::init(OsmGridRegionTree<TDummy> & grt, uint32_t gridLatCount, uint32_t gridLonCount, uint32_t cellSizeTh, uint32_t threadCount) {
	if (!threadCount) {
		threadCount = std::thread::hardware_concurrency();
	}

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
	}
	
	
	//we now have to assign every face its cellid
	//a face lives in multiple regions, so every face has a unique list of region ids
	//faces that are connected and have the same region-list get the same cellid

	{
		struct Context {
			std::unordered_map<RegionList, uint32_t> cellListToCellId;
			OsmGridRegionTree<TDummy> * grt;
			RegionListContainer * p_cellLists;
			FaceCellIdMap * p_faceToCellId;
			sserialize::ProgressInfo pinfo;
			uint32_t finishedFaces;
			Triangulation::Finite_faces_iterator facesIt;
			Triangulation::Finite_faces_iterator facesEnd;
			std::mutex lock;
		} ctx;
		ctx.grt = &grt;
		ctx.p_cellLists = &m_cellLists;
		ctx.p_faceToCellId = &m_faceToCellId;
		ctx.finishedFaces = 0;
		ctx.facesIt = m_grid.tds().finite_faces_begin();
		ctx.facesEnd = m_grid.tds().finite_faces_end();
		
		//set the infinite_face to cellId=0xFFFFFFFF
		m_faceToCellId[m_grid.tds().infinite_face()] = 0xFFFFFFFF;
		//cells that are not in any region get cellid 0
		ctx.cellListToCellId[RegionList(ctx.p_cellLists, 0, 0)] = 0;
		
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
						faceCellId = ctx->cellListToCellId.size(); //+1 account for the infinte_face
						sserialize::MMVector<uint32_t>::size_type off = ctx->p_cellLists->size();
						ctx->p_cellLists->push_back(tmp.begin(), tmp.end());
						tmp = RegionList(ctx->p_cellLists, off, tmp.size());
						ctx->cellListToCellId[tmp] = faceCellId;
					}
					else {
						faceCellId = ctx->cellListToCellId.at(tmp);
					}
					assert((tmpCellList.size() || faceCellId == 0) && (faceCellId != 0 || !tmpCellList.size()));
					(*(ctx->p_faceToCellId))[fh] = faceCellId;
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
	
	std::cout << "Found " << m_cellIdToCellList.size() << " unrefined cells" << std::endl;
	//now every cell has an id but cells that are not connected may not have different cells
	//we now have to check for each id if the correspondig faces are all connected through cells with the same id
	//this essential is a graph traversel to get all connected components where each face is a node and there's an edge between nodes
	//if their correspondig faces are neighbours and share the same id
	std::cout << "Refining cells..." << std::flush;
	{
		FaceCellIdMap tmp;
		std::vector<Face_handle> stack;
		uint32_t cellId = 1; //faces that are not in any region are in cellid 0, no need to refine them
		m_refinedCellIdToUnrefined.push_back(0);
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
	//check if there are any cells that are too large
	if (cellSizeTh < std::numeric_limits<uint32_t>::max()) {
		sserialize::TimeMeasurer tm;
		tm.begin();
		std::cout << "Splitting cells larger than " << cellSizeTh << " triangles" << std::endl;
		std::vector<uint32_t> cellSizes(m_refinedCellIdToUnrefined.size(), 0);
		std::vector<Face_handle> cellRep(m_refinedCellIdToUnrefined.size());
		for(CDT::Finite_faces_iterator it(m_grid.tds().finite_faces_begin()), end(m_grid.tds().finite_faces_end()); it != end; ++it) {
			uint32_t faceCellId = m_faceToCellId[it];
			if (!cellSizes.at(faceCellId)) {
				cellRep.at(faceCellId) = it;
			}
			cellSizes.at(faceCellId) += 1;
		}

		//Stuff needed to handle the explicit dual-graph
		CellGraph cg;
		std::vector< std::pair<uint32_t, uint32_t> > bfsTree;
		CellGraph::SimpleBitVector processedBfsTreeNodes;
		std::vector<uint32_t> stack;
		
		for(uint32_t round(0); true; ++round) {
			std::cout << "Round " << round << std::endl;
			uint32_t prevCellIdCount = cellRep.size();
			//skipt cellId=0 since that is the infinite_face
			sserialize::ProgressInfo pinfo;
			pinfo.begin(cellRep.size(), "Splitting");
			for(uint32_t cellId(1); cellId < cellRep.size(); ++cellId) {
				if (cellSizes.at(cellId) < cellSizeTh) {
					continue;
				}
				uint32_t processedNodesCount = 0;
				
				cellGraph(cellRep.at(cellId), cg);
				assert(cg.size() == cellSizes.at(cellId));
				
				processedBfsTreeNodes.resize(cg.size());
				processedBfsTreeNodes.reset();
				
				//get the face with maximum hop-distance, this is currently O(n^2)
				cg.calcMaxHopDistance(bfsTree);
				//set new cellRep
				cellRep.at(cellId) = cg.face(bfsTree.front().first);
				
				//now split the cell by adding all faces with a maximum distance of maxDist/2 into the current cell
				//and all remaining faces into newly connected cells
				//the faces in bfsTree are already sorted according to the hop-distance
				uint32_t maxHopDistance = bfsTree.back().second / 2;
				auto bfsIt (bfsTree.begin()), bfsEnd(bfsTree.end());
				while(bfsIt->second <= maxHopDistance) {
					assert(!processedBfsTreeNodes.isSet(bfsIt->first));
					processedBfsTreeNodes.set(bfsIt->first);++processedNodesCount;
					++bfsIt;
				}
				//update the size of the old cell
				cellSizes.at(cellId) = bfsIt-bfsTree.begin();
				//bfsIt points to the first node with hopdistance > maxHopDistance
				//reset the value of the map to assign new ids to the faces
				
				//now assign new ids to all faces that are not part of the old cell anymore
				for(; bfsIt != bfsEnd; ++bfsIt) {
					if (processedBfsTreeNodes.isSet(bfsIt->first)) {
						continue;
					}
					//create the connected components for this cell
					uint32_t newCellId = m_refinedCellIdToUnrefined.size();
					uint32_t newCellSize = 0;
					m_refinedCellIdToUnrefined.push_back(m_refinedCellIdToUnrefined.at(cellId));
					cellRep.push_back(cg.face(bfsIt->first));
					stack.clear();
					stack.push_back(bfsIt->first);
					processedBfsTreeNodes.set(bfsIt->first);
					while(stack.size()) {
						uint32_t nodeId = stack.back();
						stack.pop_back();
						const CellGraph::FaceNode & fn = cg.node(nodeId);
						
						m_faceToCellId[ cg.face(nodeId) ] = newCellId;
						++newCellSize;
						
						for(int i(0); i < 3; ++i) {
							uint32_t neighborNodeId = fn.neighbours[i];
							if (neighborNodeId != CellGraph::FaceNode::NullNeighbor && !processedBfsTreeNodes.isSet(neighborNodeId)) {
								stack.push_back(neighborNodeId);
								assert(!processedBfsTreeNodes.isSet(neighborNodeId));
								processedBfsTreeNodes.set(neighborNodeId);++processedNodesCount;
							}
						}
					}
					cellSizes.push_back(newCellSize);
				}
				pinfo(cellId);
				assert(processedNodesCount == cg.size());
				
#if defined(DEBUG_CHECK_ALL) || !defined(NDEBUG)
		{
			assert(cellSizes.size() == cellCount());
			std::vector<uint32_t> triangCountOfCells(cellCount(), 0);
			for(Finite_faces_iterator it(finite_faces_begin()), end(finite_faces_end()); it != end; ++it) {
				uint32_t fid = this->cellId(it);
				triangCountOfCells.at(fid) += 1;
			}
			for(uint32_t cellId(0), s(cellCount()); cellId < s; ++cellId) {
				assert(cellSizes.at(cellId) == triangCountOfCells.at(cellId));
			}
		}
#endif
				
				
				
			}
			pinfo.end();
			assert(cellRep.size() == cellSizes.size());
			assert(cellRep.size() == m_refinedCellIdToUnrefined.size());
			//if no new cell was created then all cells are smaller than cellSizeTh
			if (prevCellIdCount == cellRep.size()) {
				break;
			}
		}
		tm.end();
		std::cout << "Took " << tm << " to split the cells" << std::endl;
		std::cout << "Found " << m_refinedCellIdToUnrefined.size() << " cells" << std::endl;
#if defined(DEBUG_CHECK_ALL) || !defined(NDEBUG)
		for(uint32_t x : cellSizes) {
			assert(x <= cellSizeTh);
		}
		{
			assert(cellSizes.size() == cellCount());
			std::vector<uint32_t> triangCountOfCells(cellCount(), 0);
			for(Finite_faces_iterator it(finite_faces_begin()), end(finite_faces_end()); it != end; ++it) {
				uint32_t fid = cellId(it);
				triangCountOfCells.at(fid) += 1;
			}
			for(uint32_t cellId(0), s(cellCount()); cellId < s; ++cellId) {
				assert(cellSizes.at(cellId) == triangCountOfCells.at(cellId));
				assert(triangCountOfCells.at(cellId) <= cellSizeTh && (cellId == 0 || triangCountOfCells.at(cellId)));
			}
		}
#endif
	}
	
	m_grid.initGrid(gridLatCount, gridLonCount);
}

///By definition: items that are not in any cell are in cell 0
uint32_t OsmTriangulationRegionStore::cellId(double lat, double lon) {
	std::unique_lock<std::mutex> lck(m_lock);
	Face_handle fh = m_grid.locate(lat, lon);
	if (fh->is_valid() && m_faceToCellId.is_defined(fh)) {
		return m_faceToCellId[fh];
	}
	else {
		return 0;
	}
}

const OsmTriangulationRegionStore::RegionList& OsmTriangulationRegionStore::regions(uint32_t cellId) {
	return m_cellIdToCellList.at(m_refinedCellIdToUnrefined.at(cellId) );
}

}//end namespace

#endif