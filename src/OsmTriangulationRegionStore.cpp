#include <osmtools/OsmTriangulationRegionStore.h>
#include <sserialize/utility/ThreadPool.h>

namespace osmtools {
namespace detail {
namespace OsmTriangulationRegionStore {

//definitions
constexpr uint32_t CTGraphBase::NullFace;
constexpr uint32_t CTGraphBase::FaceNode::NullNeighbor;

void CTGraphBase::calcMaxHopDistance(uint32_t & maxHopDistRoot) {
	struct WorkContext {
		std::mutex lock;
		uint32_t processNode;
		uint32_t maxHopDist;
		uint32_t maxHopDistRoot;
		NodesContainer * nodes;
	};
	
	struct Worker {
		WorkContext * ctx;
		SimpleBitVector processedNodes;
		std::vector< std::pair<uint32_t, uint32_t> > bfsTree; //holds (nodeId, hopDist)
		Worker(WorkContext * wctx) : ctx(wctx) {
			uint32_t nodeCount = ctx->nodes->size();
			processedNodes.resize(nodeCount);
			bfsTree.reserve(nodeCount);
		}
		Worker(const Worker & other) : Worker(other.ctx) {}
		void calc(uint32_t rootNodeId) {
			processedNodes.reset();
			bfsTree.clear();
			bfsTree.emplace_back(rootNodeId, 0);
			processedNodes.set(rootNodeId);
			for(uint32_t i(0); i < bfsTree.size(); ++i) {
				std::pair<uint32_t, uint32_t> & nodeInfo = bfsTree[i];
				const FaceNode & fn = (*(ctx->nodes))[nodeInfo.first];
				for(int j(0); j < 3; ++j) {
					uint32_t nid = fn.neighbours[j];
					if (nid != FaceNode::NullNeighbor && !processedNodes.isSet(nid)) {
						bfsTree.emplace_back(nid, nodeInfo.second+1);
						processedNodes.set(nid);
					}
				}
			}
		}
		void operator()() {
			while(true) {
				std::unique_lock<std::mutex> lck(ctx->lock);
				if (ctx->processNode >= ctx->nodes->size()) {
					return;
				}
				uint32_t myRootNode = ctx->processNode;
				ctx->processNode += 1;
				lck.unlock();
				calc(myRootNode);
				lck.lock();
				uint32_t maxHopDist = bfsTree.back().second;
				if (maxHopDist > ctx->maxHopDist) {
					ctx->maxHopDist = maxHopDist;
					ctx->maxHopDistRoot = myRootNode;
				}
			}
		}
	};
	WorkContext wctx;
	wctx.processNode = 0;
	wctx.maxHopDist = 0;
	wctx.maxHopDistRoot = 0;
	wctx.nodes = &m_nodes;

	Worker w(&wctx);
	sserialize::ThreadPool::execute(w, (wctx.nodes->size() > 1000 ? 0 : 1));
	
	maxHopDistRoot = wctx.maxHopDistRoot;
}

CellGraph::CellGraph(CellGraph && other) :
m_nodePtrs(other.m_nodePtrs),
m_nodes(std::forward<NodesContainer>(other.m_nodes))
{
	other.m_nodePtrs = 0;
}

CellGraph::CellGraph(const detail::OsmTriangulationRegionStore::CellGraph & other) :
m_nodePtrs(new NodePointersContainer(*other.m_nodePtrs)),
m_nodes(other.m_nodes)
{
	for(CellNode & cn : m_nodes) {
		cn.rebind(m_nodePtrs);
	}
}

CellGraph::~CellGraph() {
	delete m_nodePtrs;
}

CellGraph& CellGraph::operator=(CellGraph && other) {
	m_nodePtrs = other.m_nodePtrs;
	other.m_nodePtrs = 0;
	
	m_nodes = std::move(other.m_nodes);
	return *this;
}

CellGraph& CellGraph::operator=(const CellGraph & other) {
	m_nodePtrs = new NodePointersContainer(*other.m_nodePtrs);
	m_nodes = other.m_nodes;
	for(CellNode & cn : m_nodes) {
		cn.rebind(m_nodePtrs);
	}
	return *this;
}

}}//end namespace detail::OsmTriangulationRegionStore

OsmTriangulationRegionStore::Face_handle OsmTriangulationRegionStore::CTGraph::face(uint32_t faceNodeId) {
	return m_faces.at(faceNodeId);
}

uint32_t OsmTriangulationRegionStore::CTGraph::node(const OsmTriangulationRegionStore::Face_handle& fh) {
	if (m_faceToNodeId.is_defined(fh)) {
		return m_faceToNodeId[fh];
	}
	throw std::out_of_range("OsmTriangulationRegionStore::CellGraph::node");
	return CTGraph::NullFace;
}


OsmTriangulationRegionStore::Point OsmTriangulationRegionStore::centroid(const OsmTriangulationRegionStore::Face_handle& fh) {
	return CGAL::centroid(fh->vertex(0)->point(), fh->vertex(1)->point(), fh->vertex(2)->point());
}

void OsmTriangulationRegionStore::cellInfo(std::vector<Face_handle> & cellRep, std::vector<uint32_t> & cellSizes) {
	cellRep.clear();
	cellSizes.clear();
	cellSizes.resize(m_refinedCellIdToUnrefined.size(), 0);
	cellRep.resize(m_refinedCellIdToUnrefined.size());
	for(CDT::Finite_faces_iterator it(m_grid.tds().finite_faces_begin()), end(m_grid.tds().finite_faces_end()); it != end; ++it) {
		uint32_t faceCellId = m_faceToCellId[it];
		if (!cellSizes.at(faceCellId)) {
			cellRep.at(faceCellId) = it;
		}
		cellSizes.at(faceCellId) += 1;
	}
}

void OsmTriangulationRegionStore::ctGraph(const Face_handle & rfh, CTGraph& cg) {
	std::vector<Face_handle> & cgFaces = cg.m_faces;
	CGAL::Unique_hash_map<Face_handle, uint32_t> & faceToNodeId = cg.m_faceToNodeId;
	
	faceToNodeId.clear();
	cgFaces.clear();
	cg.m_nodes.clear();

	if (!m_faceToCellId.is_defined(rfh)) {
		throw std::out_of_range("OsmTriangulationRegionStore::ctGraph called with invalid cell representative");
	}
	
	uint32_t myCellId = m_faceToCellId[rfh];
	cg.m_cellId = myCellId;
	
	if (m_isConnected) {
		cgFaces.push_back(rfh);
		faceToNodeId[rfh] = 0;
		for(uint32_t i(0); i < cgFaces.size(); ++i) {
			Face_handle fh = cgFaces[i];
			for(int j(0); j < 3; ++j) {
				Face_handle nfh = fh->neighbor(j);
				assert(m_faceToCellId.is_defined(nfh));
				if (m_faceToCellId[nfh] == myCellId && !faceToNodeId.is_defined(nfh)) {
					faceToNodeId[nfh] = cgFaces.size();
					cgFaces.emplace_back(nfh);
				}
			}
		}
	}
	else {//unrefined graph may have unconnected cells
		for(Finite_faces_iterator it(finite_faces_begin()), end(finite_faces_end()); it != end; ++it) {
			Face_handle fh = it;
			for(int j(0); j < 3; ++j) {
				Face_handle nfh = fh->neighbor(j);
				assert(m_faceToCellId.is_defined(nfh));
				if (m_faceToCellId[nfh] == myCellId && !faceToNodeId.is_defined(nfh)) {
					faceToNodeId[nfh] = cgFaces.size();
					cgFaces.emplace_back(nfh);
				}
			}
		}
	}
	
	for(const Face_handle & fh : cgFaces) {
		CTGraph::FaceNode fn;
		for(int i(0); i < 3; ++i) {
			Face_handle nfh = fh->neighbor(i);
			if (faceToNodeId.is_defined(nfh)) {
				fn.neighbours[i] = faceToNodeId[nfh];
			}
			else {
				fn.neighbours[i] = CTGraph::FaceNode::NullNeighbor;
			}
		}
		cg.m_nodes.push_back(fn);
	}
}

OsmTriangulationRegionStore::CellGraph OsmTriangulationRegionStore::cellGraph() {
	if (!cellCount()) {
		return CellGraph();
	}
	
	std::vector< std::pair<uint32_t, uint32_t> > edges;
	{
		std::unordered_set< std::pair<uint32_t, uint32_t> > tmp;
		for(Finite_faces_iterator it(finite_faces_begin()), end(finite_faces_end()); it != end; ++it) {
			Face_handle fh = it;
			uint32_t myCellId = m_faceToCellId[fh];
			for(int j(0); j < 3; ++j) {
				Face_handle nfh = fh->neighbor(j);
				if (m_faceToCellId.is_defined(nfh)) {
					uint32_t oCellId = m_faceToCellId[nfh];
					if (oCellId != myCellId && oCellId != InfiniteFacesCellId) {
						tmp.emplace(myCellId, oCellId);
					}
				}
			}
		}
		edges.insert(edges.end(), tmp.begin(), tmp.end());
		std::sort(edges.begin(), edges.end());
	}
	//now build the graph
	CellGraph cg;
	assert(!cg.m_nodePtrs);
	cg.m_nodePtrs = new CellGraph::NodePointersContainer(sserialize::MM_PROGRAM_MEMORY);
	
	uint32_t prevNodeId = 0;
	CellGraph::NodePointersContainer::size_type prevOffset = 0;
	for(const std::pair<uint32_t, uint32_t> & e : edges) {
		if (e.first != prevNodeId) {
			assert(prevNodeId == cg.m_nodes.size());
			cg.m_nodes.emplace_back(cg.m_nodePtrs, prevOffset, cg.m_nodePtrs->size()-prevOffset);
			prevOffset = cg.m_nodePtrs->size();
			prevNodeId = e.first;
		}
		cg.m_nodePtrs->push_back(e.second);
	}
	assert(prevNodeId = cg.m_nodes.size());
	cg.m_nodes.emplace_back(cg.m_nodePtrs, prevOffset, cg.m_nodePtrs->size()-prevOffset);
	prevOffset = cg.m_nodePtrs->size();
	
	return cg;
}

void OsmTriangulationRegionStore::
hopDistances(const Face_handle & rfh, std::vector<Face_handle> & cellTriangs, CGAL::Unique_hash_map<Face_handle, uint32_t> & cellTriangMap, uint32_t & maxHopDist) {
	cellTriangMap.clear();
	cellTriangs.clear();
	maxHopDist = 0;
	cellTriangs.emplace_back(rfh);
	cellTriangMap[rfh] = 0;
	uint32_t fhId = m_faceToCellId[rfh];
	for(uint32_t i(0); i < cellTriangs.size(); ++i) {
		Face_handle fh = cellTriangs.at(i);
		assert(fh->is_valid());
		uint32_t fhHopDist = cellTriangMap[fh];
		for(int j=0; j < 3; ++j) {
			Face_handle nfh = fh->neighbor(j);
			if (!cellTriangMap.is_defined(nfh) && (m_faceToCellId.is_defined(nfh) && m_faceToCellId[nfh] == fhId)) {
				cellTriangs.emplace_back(nfh);
				cellTriangMap[nfh] = fhHopDist+1;
				maxHopDist = std::max<uint32_t>(maxHopDist, fhHopDist+1);
			}
		}
	}
}

void OsmTriangulationRegionStore::setInfinteFacesCellIds() {
	//set all infinte faces to cellId=0xFFFFFFFF
	for(All_faces_iterator it(m_grid.tds().all_faces_begin()), end(m_grid.tds().all_faces_end()); it != end; ++it) {
		if (m_grid.tds().is_infinite(it)) {
			m_faceToCellId[it] = InfiniteFacesCellId;
		}
	}
	//the following should be equivalent but isn't. Why?
// 	Triangulation::Face_circulator ifcBegin(m_grid.tds().incident_faces(m_grid.tds().infinite_vertex()));
// 	Triangulation::Face_circulator ifcEnd(ifcBegin);
// 	for(; ifcBegin != ifcEnd; ++ifcBegin) {
// 		m_faceToCellId[ifcBegin] = InfiniteFacesCellId;
// 	}
}


uint32_t OsmTriangulationRegionStore::cellId(const OsmTriangulationRegionStore::Face_handle& fh) {
	if (m_faceToCellId.is_defined(fh)) {
		return m_faceToCellId[fh];
	}
	throw std::out_of_range("OsmTriangulationRegionStore::cellId");
}

void OsmTriangulationRegionStore::clear() {

	m_grid = GridLocator();
	m_faceToCellId = FaceCellIdMap();
	m_cellLists = RegionListContainer();
	m_cellIdToCellList = decltype(m_cellIdToCellList)();
	m_refinedCellIdToUnrefined = decltype(m_refinedCellIdToUnrefined)();
	m_isConnected = false;
}

void OsmTriangulationRegionStore::clearRefinement() {
	for(Finite_faces_iterator it(finite_faces_begin()), end(finite_faces_end()); it != end; ++it) {
		assert(m_faceToCellId.is_defined(it));
		m_faceToCellId[it] = m_refinedCellIdToUnrefined.at(m_faceToCellId[it]);
	}
	m_refinedCellIdToUnrefined.clear();
	for(uint32_t i(0), s(m_cellIdToCellList.size()); i < s; ++i) {
		m_refinedCellIdToUnrefined.push_back(i);
	}
	m_isConnected = false;
	assert(selfTest());
}

void OsmTriangulationRegionStore::initGrid(uint32_t gridLatCount, uint32_t gridLonCount) {
	m_grid.initGrid(gridLatCount, gridLonCount);
	assert(selfTest());
}

void OsmTriangulationRegionStore::makeConnected() {
	if (m_isConnected) {
		return;
	}
	//Every cell has an id but cells that are not connected may not have different cells
	//we now have to check for each id if the correspondig faces are all connected through cells with the same id
	//this essential is a graph traversel to get all connected components where each face is a node and there's an edge between nodes
	//if their correspondig faces are neighbours and share the same id

	m_refinedCellIdToUnrefined.clear();
	
	std::cout << "Refining cells..." << std::flush;
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
	setInfinteFacesCellIds();
	
	std::cout << "done" << std::endl;
	std::cout << "Found " << cellId << " cells" << std::endl;
	
	m_isConnected = true;
	assert(selfTest());
}


void OsmTriangulationRegionStore::refineBySize(uint32_t cellSizeTh, uint32_t runs, uint32_t splitPerRun, uint32_t /*threadCount*/) {
	if (cellSizeTh > m_grid.tds().number_of_faces()) {
		return;
	}
	
	splitPerRun = std::max<uint32_t>(splitPerRun, 2);
	
	makeConnected();
	//all cells are connected now

	//check if there are any cells that are too large
	sserialize::TimeMeasurer tm;
	tm.begin();
	std::cout << "Splitting cells larger than " << cellSizeTh << " triangles" << std::endl;
	std::vector<uint32_t> cellSizes;
	std::vector<Face_handle> cellRep;
	cellInfo(cellRep, cellSizes);

	//Stuff needed to handle the explicit dual-graph
	CTGraph cg;
	std::vector<uint32_t> hopDists;
	std::vector<uint32_t> newFaceCellIds;
	std::unordered_set<uint32_t> currentCells;
	std::vector< std::pair<uint32_t, uint32_t> > stack;
	
	for(uint32_t round(0); round < runs; ++round) {
		std::cout << "Round " << round << std::endl;
		uint32_t prevCellIdCount = cellRep.size();
		//skipt cellId=0 since that is the infinite_face
		sserialize::ProgressInfo pinfo;
		pinfo.begin(cellRep.size(), "Splitting");
		for(uint32_t cellId(1), cellIdInitialSize(cellRep.size()); cellId < cellIdInitialSize; ++cellId) {
			if (cellSizes.at(cellId) < cellSizeTh) {
				continue;
			}
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
			currentCells.clear();
			newFaceCellIds.clear();
			stack.clear();
			hopDists.clear();
			
			newFaceCellIds.resize(cellSizes.at(cellId), std::numeric_limits<uint32_t>::max());
			hopDists.resize(cellSizes.at(cellId), std::numeric_limits<uint32_t>::max());
			
			assert(m_faceToCellId.is_defined(cellRep.at(cellId)) && m_faceToCellId[cellRep.at(cellId)] == cellId);
// 			uint32_t mpt = m_faceToCellId[cellRep.at(cellId)];
			assert(m_faceToCellId[cellRep.at(cellId)] == cellId);
			ctGraph(cellRep.at(cellId), cg);
			assert(cg.size() == cellSizes.at(cellId));
			
			uint32_t currentGenerator;
			cg.calcMaxHopDistance(currentGenerator);
			
			uint32_t currentCellId = cellId;
			cellRep.at(cellId) = cg.face(currentGenerator);
			currentCells.insert(currentCellId);
			
			bool cellsTooLarge = true;
			for(uint32_t voronoiSplitRun(0); voronoiSplitRun < splitPerRun && cellsTooLarge; ++voronoiSplitRun) {
				//do a depth-first search and mark all nodes with larger hopDists as our own node
				stack.clear();
				stack.emplace_back(currentGenerator, 0); //nodeId, next-neighbor to inspect
				hopDists.at(currentGenerator) = 0;
				newFaceCellIds.at(currentGenerator) = currentCellId;
				while (stack.size()) {
					while(stack.size() && stack.back().second == 3) {
						stack.pop_back();
					}
					if (!stack.size()) {
						break;
					}
					uint32_t nextHopDist = stack.size();
					
					std::pair<uint32_t, uint32_t> & cn = stack.back();
					const CTGraph::FaceNode & fn = cg.node(cn.first);
					uint32_t nid = fn.neighbours[cn.second];
					
					cn.second += 1;
					if (nid != CTGraph::FaceNode::NullNeighbor && hopDists.at(nid) > nextHopDist) {
						hopDists.at(nid) = nextHopDist;
						newFaceCellIds.at(nid) = currentCellId;
						stack.emplace_back(nid, 0);
					}
					assert(m_faceToCellId.is_defined(cellRep.at(cellId)) && m_faceToCellId[cellRep.at(cellId)] == cellId);
				}
				
				for(uint32_t x : currentCells) {
					cellSizes.at(x) = 0;
				}
				for(uint32_t & x : newFaceCellIds) {
					cellSizes.at(x) += 1;
				}
				cellsTooLarge = false;
				for(uint32_t x : currentCells) {
					assert(cellSizes.at(x));
					if (cellSizes.at(x) > cellSizeTh) {
						cellsTooLarge = true;
						break;
					}
				}
				if (cellsTooLarge && voronoiSplitRun+1 < splitPerRun) {//find a new generator
					std::vector<uint32_t>::const_iterator maxElem = std::max_element(hopDists.begin(), hopDists.end());
					currentGenerator = maxElem - hopDists.begin();
					currentCellId = m_refinedCellIdToUnrefined.size();
					m_refinedCellIdToUnrefined.push_back(m_refinedCellIdToUnrefined.at(cellId));
					cellSizes.push_back(0);
					currentCells.insert(currentCellId);
					cellRep.push_back(cg.face(currentGenerator));
				}
			}//end for-loop voronoi-split run
			assert(m_faceToCellId.is_defined(cellRep.at(cellId)) && m_faceToCellId[cellRep.at(cellId)] == cellId);
			assert(currentCells.size() <= splitPerRun);
			//cellSizes are correctly set, assign faces the new ids
			for(uint32_t nodeId(0), s(cg.size()); nodeId < s; ++nodeId) {
				assert(newFaceCellIds.at(nodeId) != std::numeric_limits<uint32_t>::max());
				m_faceToCellId[cg.face(nodeId)] = newFaceCellIds.at(nodeId);
			}
			assert(m_faceToCellId.is_defined(cellRep.at(cellId)) && m_faceToCellId[cellRep.at(cellId)] == cellId);
			
			pinfo(cellId);
			
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
		}//end for-loop cell-loop
		pinfo.end();
		assert(cellRep.size() == cellSizes.size());
		assert(cellRep.size() == m_refinedCellIdToUnrefined.size());
		//if no new cell was created then all cells are smaller than cellSizeTh
		if (prevCellIdCount == cellRep.size()) {
			break;
		}
	}//end for-loop split rounds
	tm.end();
	std::cout << "Took " << tm << " to split the cells" << std::endl;
	std::cout << "Found " << m_refinedCellIdToUnrefined.size() << " cells" << std::endl;
#if defined(DEBUG_CHECK_ALL) || !defined(NDEBUG)
	if (runs == 0xFFFFFFFF) {
		for(uint32_t x : cellSizes) {
			assert(x <= cellSizeTh);
		}
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
			assert((runs != 0xFFFFFFFF || triangCountOfCells.at(cellId) <= cellSizeTh) && (cellId == 0 || triangCountOfCells.at(cellId)));
		}
	}
#endif
	assert(selfTest());
}

OsmTriangulationRegionStore::OsmTriangulationRegionStore() :
m_isConnected(false)
{}

void OsmTriangulationRegionStore::printStats(std::ostream& out) {
	if (cellCount() <= 1)
		return;
	std::vector<uint32_t> triangCountOfCells(cellCount(), 0);
	for(Finite_faces_iterator it(finite_faces_begin()), end(finite_faces_end()); it != end; ++it) {
		uint32_t fid = cellId(it);
		triangCountOfCells.at(fid) += 1;
	}
	//skip cell 0 since that one is not created by any region ans thus should not contain many or any items
	
	std::vector<uint32_t>::const_iterator maxElem = std::max_element(triangCountOfCells.begin()+1, triangCountOfCells.end());
	std::vector<uint32_t>::const_iterator minElem = std::min_element(triangCountOfCells.begin()+1, triangCountOfCells.end());
	
	out << "Cell Triangle stats: \n";
	out << "\tmin: " << *minElem << " at " << minElem - triangCountOfCells.begin() << "\n";
	out << "\tmax: " << *maxElem << " at " << maxElem - triangCountOfCells.begin() << "\n";
	out << "\tmedian: " << sserialize::statistics::median(triangCountOfCells.begin()+1, triangCountOfCells.end(), 0) << "\n";
	out << "\tmean: " << sserialize::statistics::mean(triangCountOfCells.begin()+1, triangCountOfCells.end(), 0) << "\n";
	out << std::flush;
}

///By definition: items that are not in any cell are in cell 0
uint32_t OsmTriangulationRegionStore::cellId(double lat, double lon) {
	uint32_t cellId;
	if (m_grid.contains(lat, lon)) {
		std::unique_lock<std::mutex> lck(m_lock);
		Face_handle fh = m_grid.locate(lat, lon);
		assert(m_faceToCellId.is_defined(fh));
		cellId = m_faceToCellId[fh];
		lck.unlock();
		if (cellId == InfiniteFacesCellId) {
			cellId = 0;
		}
	}
	else {
		cellId = 0;
	}
	return cellId;
}

const OsmTriangulationRegionStore::RegionList& OsmTriangulationRegionStore::regions(uint32_t cellId) {
	return m_cellIdToCellList.at(m_refinedCellIdToUnrefined.at(cellId) );
}

bool OsmTriangulationRegionStore::selfTest() {
	for(All_faces_iterator it(m_grid.tds().all_faces_begin()), end(m_grid.tds().all_faces_end()); it != end; ++it) {
		if (!m_faceToCellId.is_defined(it)) {
			if (m_grid.tds().is_infinite(it)) {
				return false;
			}
			return false;
		}
	}
	
	for(const Face_handle & fh : m_grid.grid().storage()) {
		if (!m_faceToCellId.is_defined(fh)) {
			return false;
		}
	}
	if (m_isConnected) {
		std::vector<Face_handle> cellRepresentatives;
		std::vector<uint32_t> cellSizes;
		cellInfo(cellRepresentatives, cellSizes);
		assert(cellRepresentatives.size() == cellSizes.size() && cellSizes.size() == cellCount());
		CTGraph cg;
		for(uint32_t cellId(1), s(cellRepresentatives.size()); cellId < s; ++cellId) {
			ctGraph(cellRepresentatives.at(cellId), cg);
			if (cg.size() != cellSizes.at(cellId)) {
				return false;
			}
		}
	}
	return true;
}


}//end namespace