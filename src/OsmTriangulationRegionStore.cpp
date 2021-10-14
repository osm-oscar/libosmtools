#include <osmtools/OsmTriangulationRegionStore.h>
#include <sserialize/mt/ThreadPool.h>
#include <sserialize/Static/TracGraph.h>
#include <random>

namespace osmtools {
namespace detail {
namespace OsmTriangulationRegionStore {

bool CTGraphBase::calcMaxHopDistance(FaceId & maxHopDistRoot) {

	struct WorkContext {
		std::mutex lock;
		std::size_t processNode;
		cellsize_type maxHopDist;
		FaceId maxHopDistRoot;
		NodesContainer * nodes;
	};
	
	struct Worker {
		WorkContext * ctx;
		sserialize::SimpleBitVector processedNodes;
		std::vector< std::pair<FaceId, cellsize_type> > bfsTree; //holds (nodeId, hopDist)
		Worker(WorkContext * wctx) : ctx(wctx) {
			std::size_t nodeCount = ctx->nodes->size();
			processedNodes.resize(nodeCount);
			bfsTree.reserve(nodeCount);
		}
		Worker(const Worker & other) : Worker(other.ctx) {}
		void calc(FaceId rootNodeId) {
			processedNodes.reset();
			bfsTree.clear();
			bfsTree.emplace_back(rootNodeId, 0);
			processedNodes.set(rootNodeId.ut());
			for(std::size_t i(0); i < bfsTree.size(); ++i) {
				auto & nodeInfo = bfsTree[i];
				const FaceNode & fn = (*(ctx->nodes))[nodeInfo.first.ut()];
				for(int j(0); j < 3; ++j) {
					auto nid = fn.neighbours[j];
					if (nid != FaceNode::NullNeighbor && !processedNodes.isSet(nid.ut())) {
						bfsTree.emplace_back(nid, nodeInfo.second+1);
						processedNodes.set(nid.ut());
					}
				}
			}
		}
		//check if root is valid. Root is valid iff it has less than 3 neighbors
		//Inner nodes can never be the root of the bfs of the diameter
		//Thus rootNode musst be on the border which is the case iff it has less than 3 neighbors (since this is a triangulation!)
		bool validRoot(FaceId rootNode) const {
			const FaceNode & rn = (*(ctx->nodes))[rootNode.ut()];
			bool ok = false;
			for(int j(0); j < 3; ++j) {
				auto nid = rn.neighbours[j];
				if (nid == FaceNode::NullNeighbor) {
					ok = true;
					break;
				}
			}
			return ok;
		}
		//returns true if there are more nodes to be processed
		void process(FaceId myRootNode) {
			if (!validRoot(myRootNode)) {
				return;
			}
			calc(myRootNode);
			std::lock_guard<std::mutex> lck(ctx->lock);
			uint32_t maxHopDist = bfsTree.back().second;
			if (maxHopDist > ctx->maxHopDist) {
				ctx->maxHopDist = maxHopDist;
				ctx->maxHopDistRoot = myRootNode;
			}
		}
		
		void operator()() {
			FaceId::underlying_type myRootNode = 0;
			while(true) {
				{
					std::lock_guard<std::mutex> lck(ctx->lock);
					if (ctx->processNode >= ctx->nodes->size()) {
						return;
					}
					myRootNode = ctx->processNode;
					ctx->processNode += 1;
				}
				process( FaceId{myRootNode} );
			}
		}
	};
	WorkContext wctx;
	wctx.processNode = 0;
	wctx.maxHopDist = 0;
	wctx.maxHopDistRoot = FaceId{0};
	wctx.nodes = &m_nodes;
	Worker w(&wctx);

	bool hopDistIsExact = true;
	if (size() > 50*1000) {
		//approximate hop distance. We do this by selecting a random start node
		//from there we calculate a bfs-tree.
		//Use the node farthest away as a new root-node and calculate another bfs-tree
		//We do this until the maximum hop-distance doesn't change anymore
		//We use multiple runs with random start points to try to avoid local minima
		std::uniform_int_distribution<std::size_t> ud(0, size()-1);
		std::default_random_engine gen( std::chrono::system_clock::now().time_since_epoch().count() );
		
		cellsize_type prevMaxHopDist = wctx.maxHopDist;
		for(std::size_t run(0); run < 10; ++run) {
			FaceId rootNode{ ud(gen) };
			//we first have to get a root node at the border
			w.calc(rootNode);
			rootNode = w.bfsTree.back().first;
			while (true) {
				SSERIALIZE_CHEAP_ASSERT(w.validRoot(rootNode));
				w.process(rootNode);
				if (prevMaxHopDist < wctx.maxHopDist) {
					prevMaxHopDist = wctx.maxHopDist;
					rootNode = w.bfsTree.back().first;
				}
				else {
					break;
				}
			}
		}
		hopDistIsExact = false;
	}
	else {
		//This is an O(n^2) algorithm, hence if n > 10000 we can only approximate the real maximum hop distance
		sserialize::ThreadPool::execute(w, (wctx.nodes->size() > 1000 ? 0 : 1), sserialize::ThreadPool::CopyTaskTag());
	}
	maxHopDistRoot = wctx.maxHopDistRoot;
	return hopDistIsExact;
}

CellGraph::CellGraph(CellGraph && other) :
m_nodePtrs(other.m_nodePtrs),
m_nodes(std::move(other.m_nodes))
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

sserialize::UByteArrayAdapter& CellGraph::append(sserialize::UByteArrayAdapter& dest, const std::unordered_map<cellid_type, cellid_type> & myIdsToGhCellIds) const {
	#ifdef SSERIALIZE_NORMAL_ASSERT_ENABLED
	sserialize::UByteArrayAdapter::OffsetType debugBegin = dest.tellPutPtr();
	#endif
	dest.putUint8(1); //Version
	std::size_t edgeCount = 0;
	std::vector<cellid_type> ghCellIdsToMyIds(myIdsToGhCellIds.size(), 0xFFFFFFFF);
	for(const std::pair<const cellid_type, cellid_type> & x : myIdsToGhCellIds) {
		ghCellIdsToMyIds.at(x.second) = x.first;
	}
	
	std::size_t maxNeighborCount = 0;
	for(cellid_type cellId(0), s(sserialize::narrow_check<cellid_type>(size())); cellId < s; ++cellId) {
		if (!myIdsToGhCellIds.count(cellId)) {
			continue;
		}
		CellNode n(node(cellId));
		std::size_t tmp = 0;
		for(cellid_type nId : n) {
			if (myIdsToGhCellIds.count(nId)) {
				++tmp;
			}
		}
		maxNeighborCount = std::max(maxNeighborCount, tmp);
		edgeCount += tmp;
	}
	std::vector<uint8_t> bitConfig(2, 0);
	std::vector<cellid_type> edges; //offset array for node edges
	bitConfig[0] = sserialize::CompactUintArray::minStorageBits(maxNeighborCount);
	bitConfig[1] = sserialize::CompactUintArray::minStorageBits(edgeCount);
	sserialize::MultiVarBitArrayCreator mvac(bitConfig, dest);
	mvac.reserve(myIdsToGhCellIds.size());
	for(cellid_type ghCellId(0), s(sserialize::narrow_check<cellid_type>(ghCellIdsToMyIds.size())); ghCellId < s; ++ghCellId) {
		cellid_type myCellId = ghCellIdsToMyIds.at(ghCellId);
		CellNode n(node(myCellId));
		cellsize_type neighborCount = 0;
		std::size_t edgesBegin = edges.size();
		for(cellid_type nId : n) {
			if (myIdsToGhCellIds.count(nId)) {
				++neighborCount;
				edges.push_back(myIdsToGhCellIds.at(nId));
			}
		}
		mvac.set(ghCellId, 0, neighborCount);
		mvac.set(ghCellId, 1, edgesBegin);
	}
	#ifdef SSERIALIZE_NORMAL_ASSERT_ENABLED
	{
		sserialize::MultiVarBitArray mvat(mvac.flush());
		SSERIALIZE_NORMAL_ASSERT(mvat.size() == myIdsToGhCellIds.size());
	}
	#else
	mvac.flush();
	#endif
	SSERIALIZE_CHEAP_ASSERT(edgeCount == edges.size());
	sserialize::BoundedCompactUintArray::create(edges, dest);

	#ifdef SSERIALIZE_NORMAL_ASSERT_ENABLED
	sserialize::UByteArrayAdapter tmp(dest);
	tmp.setPutPtr(debugBegin);
	tmp.shrinkToPutPtr();
	sserialize::Static::spatial::TracGraph sg(tmp);
	SSERIALIZE_NORMAL_ASSERT(sg.size() == myIdsToGhCellIds.size());
	#endif
	return dest;
}



}}//end namespace detail::OsmTriangulationRegionStore


bool OsmTriangulationRegionStore::CellCriteriaInterface::refine(const State & state) {
	for(cellid_type x : state.currentCells) {
		SSERIALIZE_CHEAP_ASSERT(state.cellSizes.at(x));
		if (this->refine(x, state)) {
			return true;
		}
	}
	return false;
}

OsmTriangulationRegionStore::FaceInfo::FaceInfo() :
m_cellId(OsmTriangulationRegionStore::UnsetFacesCellId)
{}

void OsmTriangulationRegionStore::FaceInfo::clear() {
	m_cellId = OsmTriangulationRegionStore::UnsetFacesCellId;
}

void OsmTriangulationRegionStore::FaceInfo::setCellId(cellid_type cellId) {
	m_cellId = cellId;
}

uint32_t OsmTriangulationRegionStore::FaceInfo::cellId() const {
	return m_cellId;
}

bool OsmTriangulationRegionStore::FaceInfo::hasCellId() const {
	return m_cellId != OsmTriangulationRegionStore::UnsetFacesCellId;
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
		uint32_t faceCellId = it->info().cellId();
		if (!cellSizes.at(faceCellId)) {
			cellRep.at(faceCellId) = it;
		}
		cellSizes.at(faceCellId) += 1;
	}
}

std::vector< sserialize::spatial::GeoPoint > OsmTriangulationRegionStore::cellCenterOfMass(const std::unordered_map<uint32_t, uint32_t> & myIdsToGhCellIds) {
	std::vector< sserialize::spatial::GeoPoint> centerOfMass(myIdsToGhCellIds.size(), sserialize::spatial::GeoPoint(0.0, 0.0));
	std::vector<uint32_t> faceCount(cellCount(), 0);
	for(auto it(finite_faces_begin()), end(finite_faces_end()); it != end; ++it) {
		Point fp = centroid(it);
		double lat = CGAL::to_double(fp.x());
		double lon = CGAL::to_double(fp.y());
		uint32_t myCellId = cellId(it);
		if (myIdsToGhCellIds.count(myCellId)) {
			uint32_t ghCellId = myIdsToGhCellIds.at(myCellId);
			auto & p = centerOfMass.at(ghCellId);
			p.lat() += lat;
			p.lon() += lon;
			faceCount.at(ghCellId) += 1;
		}
	}
	//reweight
	for(std::size_t i(0), s(centerOfMass.size()); i < s; ++i) {
		uint32_t fc = faceCount.at(i);
		if (fc) {
			auto & x = centerOfMass.at(i);
			x.lat() /= fc;
			x.lon() /= fc;
		}
	}
	return centerOfMass;
}

void OsmTriangulationRegionStore::ctGraph(const Face_handle & rfh, CTGraph& cg) {
	std::vector<Face_handle> & cgFaces = cg.m_faces;
	CGAL::Unique_hash_map<Face_handle, FaceId> & faceToNodeId = cg.m_faceToNodeId;
	
	faceToNodeId.clear();
	cgFaces.clear();
	cg.m_nodes.clear();
	cg.m_cellId = 0xFFFFFFFF;
	
	if (!rfh->info().hasCellId()) {
		throw std::out_of_range("OsmTriangulationRegionStore::ctGraph called with invalid cell representative");
	}
	
	uint32_t myCellId = rfh->info().cellId();
	cg.m_cellId = myCellId;
	
	if (m_isConnected) {
		cgFaces.push_back(rfh);
		faceToNodeId[rfh] = FaceId(0);
		for(uint32_t i(0); i < cgFaces.size(); ++i) {
			Face_handle fh = cgFaces[i];
			for(int j(0); j < 3; ++j) {
				Face_handle nfh = fh->neighbor(j);
				SSERIALIZE_CHEAP_ASSERT(nfh->info().hasCellId());
				if (nfh->info().cellId() == myCellId && !faceToNodeId.is_defined(nfh)) {
					faceToNodeId[nfh] = FaceId(cgFaces.size());
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
				SSERIALIZE_CHEAP_ASSERT(nfh->info().hasCellId());
				if (nfh->info().cellId() == myCellId && !faceToNodeId.is_defined(nfh)) {
					faceToNodeId[nfh] = FaceId(cgFaces.size());
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
	
	SSERIALIZE_EXPENSIVE_ASSERT(selfTest());
	
	std::vector< std::pair<cellid_type, cellid_type> > edges;
	{
		std::unordered_set< decltype(edges)::value_type > tmp;
		for(Finite_faces_iterator it(finite_faces_begin()), end(finite_faces_end()); it != end; ++it) {
			Face_handle fh = it;
			uint32_t myCellId = fh->info().cellId();
			for(int j(0); j < 3; ++j) {
				Face_handle nfh = fh->neighbor(j);
				if (nfh->info().hasCellId()) {
					auto oCellId = nfh->info().cellId();
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
	SSERIALIZE_CHEAP_ASSERT(!cg.m_nodePtrs);
	cg.m_nodePtrs = new CellGraph::NodePointersContainer(sserialize::MM_PROGRAM_MEMORY);
	
	cellid_type prevNodeId = 0;
	CellGraph::NodePointersContainer::size_type prevOffset = 0;
	for(const auto & e : edges) {
		if (e.first != prevNodeId) {
			SSERIALIZE_CHEAP_ASSERT(prevNodeId == cg.m_nodes.size());
			cg.m_nodes.emplace_back(cg.m_nodePtrs, prevOffset, cg.m_nodePtrs->size()-prevOffset);
			prevOffset = cg.m_nodePtrs->size();
			prevNodeId = e.first;
		}
		cg.m_nodePtrs->push_back(e.second);
	}
	SSERIALIZE_CHEAP_ASSERT(prevNodeId == cg.m_nodes.size());
	cg.m_nodes.emplace_back(cg.m_nodePtrs, prevOffset, cg.m_nodePtrs->size()-prevOffset);
	prevOffset = cg.m_nodePtrs->size();
	
	return cg;
}

void OsmTriangulationRegionStore::
hopDistances(const Face_handle & rfh, std::vector<Face_handle> & cellTriangs, CGAL::Unique_hash_map<Face_handle, cellsize_type> & cellTriangMap, cellsize_type & maxHopDist) {
	cellTriangMap.clear();
	cellTriangs.clear();
	maxHopDist = 0;
	cellTriangs.emplace_back(rfh);
	cellTriangMap[rfh] = 0;
	cellid_type fhId = rfh->info().cellId();
	for(std::size_t i(0); i < cellTriangs.size(); ++i) {
		Face_handle fh = cellTriangs.at(i);
		SSERIALIZE_NORMAL_ASSERT(fh->is_valid());
		cellsize_type fhHopDist = cellTriangMap[fh];
		for(int j=0; j < 3; ++j) {
			Face_handle nfh = fh->neighbor(j);
			if (!cellTriangMap.is_defined(nfh) && (nfh->info().hasCellId() && nfh->info().cellId() == fhId)) {
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
			it->info().setCellId(InfiniteFacesCellId);
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
	if (fh->info().hasCellId()) {
		return fh->info().cellId();
	}
	throw std::out_of_range("OsmTriangulationRegionStore::cellId");
}

void OsmTriangulationRegionStore::clear() {
	assert( m_grid.tds().is_valid() );
	m_grid = GridLocator();
	m_cellLists = RegionListContainer();
	m_cellIdToCellList = decltype(m_cellIdToCellList)();
	m_refinedCellIdToUnrefined = decltype(m_refinedCellIdToUnrefined)();
	m_isConnected = false;
}


void OsmTriangulationRegionStore::clearCells() {
	m_cellIdToCellList.clear();
	m_cellLists.clear();
	m_refinedCellIdToUnrefined.clear();
	m_isConnected = false;
	for(All_faces_iterator it(m_grid.tds().all_faces_begin()), end(m_grid.tds().all_faces_end()); it != end; ++it) {
		it->info().clear();
	}
}

void OsmTriangulationRegionStore::clearRefinement() {
	for(Finite_faces_iterator it(finite_faces_begin()), end(finite_faces_end()); it != end; ++it) {
		if (it->info().hasCellId()) {
			it->info().setCellId( m_refinedCellIdToUnrefined.at(it->info().cellId()) );
		}
	}
	m_refinedCellIdToUnrefined.clear();
	for(uint32_t i(0), s((uint32_t) m_cellIdToCellList.size()); i < s; ++i) {
		m_refinedCellIdToUnrefined.push_back(i);
	}
	m_isConnected = false;
	m_cs &= ~CS_HAVE_REFINED_CELLS;
	SSERIALIZE_EXPENSIVE_ASSERT(selfTest());
}

void OsmTriangulationRegionStore::initGrid(std::size_t gridLatCount, std::size_t gridLonCount) {
	m_grid.initGrid(gridLatCount, gridLonCount);
	SSERIALIZE_EXPENSIVE_ASSERT(selfTest());
}

void OsmTriangulationRegionStore::makeConnected() {
	if (m_isConnected) {
		return;
	}
	
	if (!(m_cs & CS_HAVE_CELLS)) {
		throw sserialize::PreconditionViolationException("OsmTriangulationRegionStore::makeConnected: no cells set");
	}
	
	SSERIALIZE_EXPENSIVE_ASSERT(selfTest());
	//Every cell has an id but cells that are not connected may not have different cells
	//we now have to check for each id if the correspondig faces are all connected through cells with the same id
	//this essential is a graph traversel to get all connected components where each face is a node and there's an edge between nodes
	//if their correspondig faces are neighbours and share the same id

	m_refinedCellIdToUnrefined.clear();
	
	std::cout << "Refining cells..." << std::flush;
	FaceCellIdMap tmp;
	std::vector<Face_handle> stack;
	cellid_type cellId = 0;
	auto makeCC = [&tmp, &cellId, &stack, this](Face_handle rfh) {
		if (tmp.is_defined(rfh)) {
			return;
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
			SSERIALIZE_CHEAP_ASSERT(fh->info().hasCellId());
			auto fhId = fh->info().cellId();
			for(int i=0; i < 3; ++i) {
				Face_handle nfh = fh->neighbor(i);
				if (nfh->info().hasCellId() && nfh->info().cellId() == fhId && !tmp.is_defined(nfh)) {
					stack.push_back(nfh);
				}
			}
		}
		stack.clear();
		m_refinedCellIdToUnrefined.push_back(rfh->info().cellId());
		++cellId;
	};
	
	//first take care of all faces with cellId=0
	for(CDT::Finite_faces_iterator it(m_grid.tds().finite_faces_begin()), end(m_grid.tds().finite_faces_end()); it != end; ++it) {
		if (it->info().cellId() == 0) {
			makeCC(it);
		}
	}
	if (!m_refinedCellIdToUnrefined.size()) {
		m_refinedCellIdToUnrefined.push_back(0);
		++cellId;
	}
	//now take care of the rest
	for(CDT::Finite_faces_iterator it(m_grid.tds().finite_faces_begin()), end(m_grid.tds().finite_faces_end()); it != end; ++it) {
		makeCC(it);
	}
	SSERIALIZE_CHEAP_ASSERT(cellId == m_refinedCellIdToUnrefined.size());
	for(Finite_faces_iterator it(finite_faces_begin()), end(finite_faces_end()); it != end; ++it) {
		SSERIALIZE_CHEAP_ASSERT(tmp.is_defined(it));
		it->info().setCellId(tmp[it]);
	}
	setInfinteFacesCellIds();
	
	std::cout << "done" << std::endl;
	std::cout << "Found " << cellId << " cells" << std::endl;
	
	m_isConnected = true;
	SSERIALIZE_EXPENSIVE_ASSERT(selfTest());
}

void OsmTriangulationRegionStore::refineCells(std::shared_ptr<CellCriteriaInterface> refiner, std::size_t runs, std::size_t splitPerRun, std::size_t /*threadCount*/) {
	makeConnected();
	//all cells are connected now

	if (!refiner->init(*this)) {
		return;
	}
	
	int dd = refiner->dataDependence();
	CellCriteriaInterface::State state;
	
	splitPerRun = std::max<std::size_t>(splitPerRun, 2);
	
	//check if there are any cells that are too large
	sserialize::TimeMeasurer tm;
	tm.begin();
	refiner->begin();
	cellInfo(state.cellRep, state.cellSizes);

	//Stuff needed to handle the explicit dual-graph
	std::vector<cellsize_type> hopDists;
	std::vector< std::pair<FaceId, cellsize_type> > stack;
	
	for(std::size_t round(0); round < runs; ++round) {
		std::cout << "Round " << round << std::endl;
		std::size_t prevCellIdCount = state.cellRep.size();
		//skipt cellId=0 since that is the infinite_face
		sserialize::ProgressInfo pinfo;
		pinfo.begin(state.cellRep.size(), "Splitting");
		for(cellid_type cellId(1), cellIdInitialSize(cellid_type(state.cellRep.size())); cellId < cellIdInitialSize; ++cellId) {
			SSERIALIZE_CHEAP_ASSERT_EQUAL(cellid_type(state.cellRep.size()), state.cellRep.size());
			SSERIALIZE_CHEAP_ASSERT(state.cellRep.at(cellId)->info().hasCellId())
			SSERIALIZE_CHEAP_ASSERT_EQUAL(state.cellRep.at(cellId)->info().cellId(), uint32_t(cellId));
			SSERIALIZE_CHEAP_ASSERT_SMALLER(uint32_t(0), state.cellSizes.at(cellId));

			state.currentCells.clear();
			state.newFaceCellIds.clear();
			state.newCellReps.clear();
			stack.clear();
			hopDists.clear();

			if (dd & CellCriteriaInterface::DD_CELL_GRAPH) {
				ctGraph(state.cellRep.at(cellId), state.cg);
				SSERIALIZE_CHEAP_ASSERT_EQUAL(state.cg.size(), state.cellSizes.at(cellId));
			}
			
			if (dd & CellCriteriaInterface::DD_NEW_FACE_CELL_IDS) {
				state.newFaceCellIds.assign(state.cellSizes.at(cellId), cellId);
			}
			
			if (dd & (CellCriteriaInterface::DD_NEW_CELL_REPS | CellCriteriaInterface::DD_CURRENT_CELLS)) {
				throw sserialize::UnimplementedFunctionException("OsmTriangulationRegionStore::refineCells with a refiner with data dependence on DD_NEW_CELL_REPS | DD_CURRENT_CELLS");
			}
			
			if (!refiner->refine(cellId, state) ) {
				continue;
			}
			
			state.newFaceCellIds.assign(state.cellSizes.at(cellId), OsmTriangulationRegionStore::UnsetFacesCellId);
			hopDists.assign(state.cellSizes.at(cellId), std::numeric_limits<uint32_t>::max());
			
			if (! (dd & CellCriteriaInterface::DD_CELL_GRAPH) ) {
				ctGraph(state.cellRep.at(cellId), state.cg);
			}
			SSERIALIZE_CHEAP_ASSERT(state.cg.size() == state.cellSizes.at(cellId));
			
			FaceId currentGenerator;
			state.cg.calcMaxHopDistance(currentGenerator);
			
			cellid_type currentCellId = cellId;
			state.cellRep.at(cellId) = state.cg.face(currentGenerator);
			state.newCellReps.push_back(currentGenerator);
			state.currentCells.insert(currentCellId);
			
			bool cellsTooLarge = true;
			for(std::size_t voronoiSplitRun(0); voronoiSplitRun < splitPerRun && cellsTooLarge; ++voronoiSplitRun) {
				//do a depth-first search and mark all nodes with larger hopDists as our own node
				stack.clear();
				stack.emplace_back(currentGenerator, 0); //nodeId, next-neighbor to inspect
				hopDists.at(currentGenerator.ut()) = 0;
				state.newFaceCellIds.at(currentGenerator.ut()) = currentCellId;
				while (stack.size()) {
					while(stack.size() && stack.back().second == 3) {
						stack.pop_back();
					}
					if (!stack.size()) {
						break;
					}
					cellsize_type nextHopDist = cellsize_type(stack.size());
					
					auto & cn = stack.back();
					const CTGraph::FaceNode & fn = state.cg.node(cn.first.ut());
					auto nid = fn.neighbours[cn.second];
					
					cn.second += 1;
					if (nid != CTGraph::FaceNode::NullNeighbor && hopDists.at(nid.ut()) > nextHopDist) {
						hopDists.at(nid.ut()) = nextHopDist;
						state.newFaceCellIds.at(nid.ut()) = currentCellId;
						stack.emplace_back(nid, 0);
					}
					SSERIALIZE_CHEAP_ASSERT(state.cellRep.at(cellId)->info().hasCellId());
					SSERIALIZE_CHEAP_ASSERT_EQUAL(state.cellRep.at(cellId)->info().cellId(), cellId);
				}
				
				for(auto x : state.currentCells) {
					state.cellSizes.at(x) = 0;
				}
				for(auto & x : state.newFaceCellIds) {
					state.cellSizes.at(x) += 1;
				}
				cellsTooLarge = refiner->refine(state);
				if (cellsTooLarge && voronoiSplitRun+1 < splitPerRun) {//find a new generator
					auto maxElem = std::max_element(hopDists.begin(), hopDists.end());
					currentGenerator = FaceId(maxElem - hopDists.begin());
					//check if there are any triangles left, this is the case if there are triangles at least one hop-dist apart
					if (state.cellSizes.at( state.newFaceCellIds.at(currentGenerator.ut()) ) <= 1) {
						break;
					}
					currentCellId = (uint32_t) m_refinedCellIdToUnrefined.size();
					m_refinedCellIdToUnrefined.push_back(m_refinedCellIdToUnrefined.at(cellId));
					state.cellSizes.push_back(1);
					state.currentCells.insert(currentCellId);
					state.cellRep.push_back(state.cg.face(currentGenerator));
					state.newCellReps.push_back(currentGenerator);
				}
			}//end for-loop voronoi-split run
			SSERIALIZE_CHEAP_ASSERT(state.cellRep.at(cellId)->info().hasCellId());
			SSERIALIZE_CHEAP_ASSERT_EQUAL(state.cellRep.at(cellId)->info().cellId(), cellId);
			SSERIALIZE_CHEAP_ASSERT(state.currentCells.size() <= splitPerRun);
			//cellSizes are correctly set, assign faces the new ids
			for(std::size_t nodeId(0), s(state.cg.size()); nodeId < s; ++nodeId) {
				SSERIALIZE_CHEAP_ASSERT_NOT_EQUAL(std::numeric_limits<uint32_t>::max(), state.newFaceCellIds.at(nodeId));
				state.cg.face(FaceId{nodeId})->info().setCellId(state.newFaceCellIds.at(nodeId));
			}
			SSERIALIZE_CHEAP_ASSERT(state.cellRep.at(cellId)->info().hasCellId());
			SSERIALIZE_CHEAP_ASSERT_EQUAL(cellId, state.cellRep.at(cellId)->info().cellId());
			
			pinfo(cellId);
		}//end for-loop cell-loop
		pinfo.end();
		SSERIALIZE_CHEAP_ASSERT(state.cellRep.size() == state.cellSizes.size());
		SSERIALIZE_CHEAP_ASSERT(state.cellRep.size() == m_refinedCellIdToUnrefined.size());
		//if no new cell was created then all cells are smaller than cellSizeTh
		if (prevCellIdCount == state.cellRep.size()) {
			break;
		}
	}//end for-loop split rounds
	tm.end();
	std::cout << "Took " << tm << " to split the cells" << std::endl;
	std::cout << "Found " << m_refinedCellIdToUnrefined.size() << " cells" << std::endl;
#ifdef SSERIALIZE_EXPENSIVE_ASSERT_ENABLED
	{
		SSERIALIZE_EXPENSIVE_ASSERT(state.cellSizes.size() == cellCount());
		std::vector<uint32_t> triangCountOfCells(cellCount(), 0);
		for(Finite_faces_iterator it(finite_faces_begin()), end(finite_faces_end()); it != end; ++it) {
			uint32_t fid = cellId(it);
			triangCountOfCells.at(fid) += 1;
		}
		for(uint32_t cellId(0), s(cellCount()); cellId < s; ++cellId) {
			SSERIALIZE_EXPENSIVE_ASSERT(state.cellSizes.at(cellId) == triangCountOfCells.at(cellId));
			SSERIALIZE_EXPENSIVE_ASSERT(runs != 0xFFFFFFFF || !refiner->refine(cellId, state));
			SSERIALIZE_EXPENSIVE_ASSERT(cellId == 0 || triangCountOfCells.at(cellId));
		}
	}
#endif
	SSERIALIZE_EXPENSIVE_ASSERT(selfTest());
}

uint32_t OsmTriangulationRegionStore::InfiniteFacesCellId = 0xFFFFFFFF;
uint32_t OsmTriangulationRegionStore::UnsetFacesCellId = 0xFFFFFFFE;

OsmTriangulationRegionStore::OsmTriangulationRegionStore() :
m_isConnected(false)
{}

void
OsmTriangulationRegionStore::init(std::shared_ptr<OsmGridRegionTreeBase> grt, std::size_t threadCount) {
	if (!threadCount) {
		threadCount = std::thread::hardware_concurrency();
	}
	this->clear();
	m_grt = grt;
	{
		//we first need to find all relevant regions and extract their segments. This should be possible by just using the extracted regions since
		//we don't do any calculations with our points so segments with the same endpoints should stay the same in different regions
		typedef std::pair<double, double> RawGeoPoint;
		typedef typename OsmGridRegionTreeBase::GeoPolygon GeoPolygon;
		typedef typename OsmGridRegionTreeBase::GeoMultiPolygon GeoMultiPolygon;
		std::unordered_map<RawGeoPoint, uint32_t> gpToId;
		std::vector<Point> pts;
		std::vector< std::pair<uint32_t, uint32_t> > segments;
		
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
		for(sserialize::spatial::GeoRegion* r : m_grt->regions()) {
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
				if ((itGp.second < -179.0 && prevGp.second > 179.0) || (itGp.second > 179.0 && prevGp.second < -179)) {
					std::cout << "Skipped edge crossing longitude boundary(-180->180)\n";
					continue;
				}
				segments.emplace_back(gpToId.at(itGp), gpToId.at(prevGp));
			}
		};
		std::cout << "OsmTriangulationRegionStore: extracting segments..." << std::flush;
		for(sserialize::spatial::GeoRegion* r : m_grt->regions()) {
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
		std::sort(segments.begin(), segments.end());
		segments.resize(std::unique(segments.begin(), segments.end()) - segments.begin());
		
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
	m_cs = CS_HAVE_TRIANGULATION;
}

void OsmTriangulationRegionStore::refineTriangulation(TriangulationRefinementAlgorithmSelector refineAlgo) {
	if (! (m_cs & CS_HAVE_TRIANGULATION)) {
		std::cout << "No triangulation available" << std::endl;
		return;
	}
// 	CGAL::Triangulation_conformer_2<Triangulation> conform(m_grid.tds());
	switch (refineAlgo) {
	case TRAS_ConformingTriangulation:
		#if defined(LIBOSMTOOLS_OSMTRS_USE_EXACT_KERNEL)
		throw sserialize::UnsupportedFeatureException("OsmTriangulationRegionStore is built with exact kernel. Conforming triangulation needs in-exact kernel");
		#else
		conform.make_conforming_Delaunay();
		#endif
		m_cs |= CS_HAVE_REFINED_TRIANGULATION;
		refineTriangulationFinalize();
		break;
	case TRAS_GabrielTriangulation:
		#if defined(LIBOSMTOOLS_OSMTRS_USE_EXACT_KERNEL)
		throw sserialize::UnsupportedFeatureException("OsmTriangulationRegionStore is built with exact kernel. Gabriel triangulation needs in-exact kernel");
		#else
		conform.make_conforming_Delaunay();
		conform.make_conforming_Gabriel();
		#endif
		m_cs |= CS_HAVE_REFINED_TRIANGULATION;
		refineTriangulationFinalize();
		break;
	default:
		throw sserialize::UnsupportedFeatureException("Unsupported refinement algorithm: " + std::to_string(refineAlgo));
		break;
	}
}

void OsmTriangulationRegionStore::simplify() {
	#if defined(LIBOSMTOOLS_OSMTRS_USE_CONSTRAINED_TRIANGULATION_PLUS)
	throw sserialize::UnimplementedFunctionException("OsmTriangulationRegionStore::simplify() is not implemented yet.");
	#else
	throw sserialize::UnsupportedFeatureException("OsmTriangulationRegionStore::simplify() needs constrained triangulation plus.");
	#endif
}

void OsmTriangulationRegionStore::refineTriangulationFinalize() {
	if (m_cs & CS_HAVE_SNAPPED_TRIANGULATION) {
		std::cerr << "WARNING: OsmTriangulationRegionStore::refineTriangulation: triangulation was snapped before" << std::endl;
		m_cs &= ~CS_HAVE_SNAPPED_TRIANGULATION;
	}
	
	if (m_cs & CS_HAVE_REFINED_CELLS) {
		std::cerr << "WARNING: OsmTriangulationRegionStore::refineTriangulation: Removing cell refinement" << std::endl;
		clearRefinement();
	}
	
	if (m_cs & CS_HAVE_CELLS) {
		std::cerr << "WARNING: OsmTriangulationRegionStore::refineTriangulation: Removing cells" << std::endl;
		m_cs &= ~CS_HAVE_CELLS;
	}
}

///assign cellIds, while keeping the old cellIds if reUseOld is true
void OsmTriangulationRegionStore::assignCellIds(std::size_t threadCount) {
	if (! (m_cs & CS_HAVE_TRIANGULATION)) {
		return;
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
		std::shared_ptr<OsmGridRegionTreeBase> grt;
		RegionListContainer * p_cellLists;
		sserialize::ProgressInfo pinfo;
		uint32_t finishedFaces;
		Triangulation::Finite_faces_iterator facesIt;
		Triangulation::Finite_faces_iterator facesEnd;
		std::mutex iteratorLock;
		std::mutex cellListLock;
	} ctx;
	ctx.grt = m_grt;
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
	ctx.pinfo.begin(m_grid.tds().number_of_faces(), "Setting cellids");
	sserialize::ThreadPool::execute(WorkFunc(&ctx), threadCount, sserialize::ThreadPool::CopyTaskTag());
	ctx.pinfo.end();

	m_cellIdToCellList.resize(ctx.cellListToCellId.size());
	for(const auto & x : ctx.cellListToCellId) {
		m_cellIdToCellList.at(x.second) = x.first.list;
	}
	
	m_refinedCellIdToUnrefined.clear();
	for(uint32_t i(0), s((uint32_t) m_cellIdToCellList.size()); i < s; ++i) {
		m_refinedCellIdToUnrefined.push_back(i);
	}
	
	m_cs &= ~CS_HAVE_REFINED_CELLS;
	m_cs |= CS_HAVE_CELLS;
	
	SSERIALIZE_EXPENSIVE_ASSERT(selfTest());
}

void OsmTriangulationRegionStore::printStats(std::ostream& out) {
	if (cellCount() <= 1)
		return;
	std::vector<uint32_t> triangCountOfCells(cellCount(), 0);
	for(Finite_faces_iterator it(finite_faces_begin()), end(finite_faces_end()); it != end; ++it) {
		uint32_t fid = cellId(it);
		triangCountOfCells.at(fid) += 1;
	}
	//skip cell 0 since that one is not created by any region and thus should not contain many or any items
	
	std::vector<uint32_t>::const_iterator maxElem = std::max_element(triangCountOfCells.begin()+1, triangCountOfCells.end());
	std::vector<uint32_t>::const_iterator minElem = std::min_element(triangCountOfCells.begin()+1, triangCountOfCells.end());
#ifdef SSERIALIZE_HAS_LIB_DTS2
	out << "ExtendedInt64q allocation stats: " << std::endl;
	out << "# extended allocations: " << CGAL::Epeceik_ft::number_of_extended_allocations << '\n';
	out << "# allocations: " << CGAL::Epeceik_ft::number_of_allocations << '\n';
#endif
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
		Face_handle fh = m_grid.locate(lat, lon);
		SSERIALIZE_CHEAP_ASSERT(fh->info().hasCellId());
		cellId = fh->info().cellId();
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

sserialize::UByteArrayAdapter& OsmTriangulationRegionStore::append(sserialize::UByteArrayAdapter& dest, sserialize::ItemIndexFactory & idxFactory, sserialize::Static::spatial::Triangulation::GeometryCleanType gct) {
	CGAL::Unique_hash_map<Face_handle, sserialize::Static::spatial::Triangulation::FaceId> face2FaceId;
	dest.putUint8(1); //version
	m_grid.append(dest, face2FaceId, gct);
	{ //serialize the region lists
		std::vector<uint32_t> tmp;
		for(uint32_t i(0), s((uint32_t) m_cellIdToCellList.size()); i < s; ++i) {
			tmp.push_back( idxFactory.addIndex(m_cellIdToCellList.at(i)) );
		}
		sserialize::BoundedCompactUintArray::create(tmp, dest);
	}
	sserialize::BoundedCompactUintArray::create(m_refinedCellIdToUnrefined, dest);
	{//faceId->cellId
		std::vector<uint32_t> faceId2CellId(m_grid.tds().number_of_faces());
		uint32_t numFiniteFaces = 0;
		for(Finite_faces_iterator fit(m_grid.tds().finite_faces_begin()), fend(m_grid.tds().finite_faces_end()); fit != fend; ++fit) {
			++numFiniteFaces;
			SSERIALIZE_CHEAP_ASSERT(face2FaceId.is_defined(fit));
			faceId2CellId.at(face2FaceId[fit].ut()) = fit->info().cellId();
		}
		faceId2CellId.resize(numFiniteFaces);
		sserialize::BoundedCompactUintArray::create(faceId2CellId, dest);
	}
	return dest;
}

sserialize::UByteArrayAdapter& OsmTriangulationRegionStore::append(sserialize::UByteArrayAdapter& dest, const std::unordered_map< cellid_type, cellid_type >& myIdsToGhCellIds, sserialize::Static::spatial::Triangulation::GeometryCleanType gct) {
#ifdef SSERIALIZE_EXPENSIVE_ASSERT_ENABLED
	sserialize::UByteArrayAdapter::OffsetType initialOffset = dest.tellPutPtr();
#endif

	CGAL::Unique_hash_map<Face_handle, sserialize::Static::spatial::Triangulation::FaceId> face2FaceId;
	cellid_type myNullCellId = sserialize::narrow_check<cellid_type>( myIdsToGhCellIds.size() );
	dest.putUint8(2); //version
	m_grid.append(dest, face2FaceId, gct);
	
	std::vector<cellid_type> faceId2CellId(m_grid.tds().number_of_faces());
	std::size_t numFiniteFaces = 0;
	for(Finite_faces_iterator fit(m_grid.tds().finite_faces_begin()), fend(m_grid.tds().finite_faces_end()); fit != fend; ++fit) {
		++numFiniteFaces;
		SSERIALIZE_CHEAP_ASSERT(face2FaceId.is_defined(fit));
		faceId2CellId.at(face2FaceId[fit].ut()) = fit->info().cellId();
	}
	faceId2CellId.resize(numFiniteFaces);
	for(cellid_type & x : faceId2CellId) {
		if (myIdsToGhCellIds.count(x)) {
			x = myIdsToGhCellIds.at(x);
		}
		else {
			x = myNullCellId;
		}
	}
	sserialize::BoundedCompactUintArray::create(faceId2CellId, dest);
	//create the cell->face mappings
	std::vector<FaceId::underlying_type> cellId2FaceId(myNullCellId);
	for(std::size_t faceId(0), s(faceId2CellId.size()); faceId < s; ++faceId) {
		cellid_type cellId = faceId2CellId.at(faceId);
		if (cellId != myNullCellId) {
			cellId2FaceId.at(cellId) = FaceId(faceId).ut();
		}
	}
	sserialize::BoundedCompactUintArray::create(cellId2FaceId, dest);
	
#ifdef SSERIALIZE_EXPENSIVE_ASSERT_ENABLED
	{
		sserialize::UByteArrayAdapter tmp(dest);
		tmp.setPutPtr(initialOffset);
		tmp.shrinkToPutPtr();
		sserialize::Static::spatial::TriangulationGeoHierarchyArrangement ra(tmp);
		SSERIALIZE_EXPENSIVE_ASSERT(ra.cellCount() == myIdsToGhCellIds.size());
		for(Finite_faces_iterator fit(m_grid.tds().finite_faces_begin()), fend(m_grid.tds().finite_faces_end()); fit != fend; ++fit) {
			SSERIALIZE_EXPENSIVE_ASSERT(face2FaceId.is_defined(fit));
			auto sfaceId = face2FaceId[fit];
			auto myCellId = fit->info().cellId();
			cellid_type remappedCellId = myNullCellId;
			if (myIdsToGhCellIds.count(myCellId)) {
				remappedCellId = myIdsToGhCellIds.at(myCellId);
			}
			SSERIALIZE_EXPENSIVE_ASSERT(remappedCellId == ra.cellIdFromFaceId(sfaceId));
		}
		SSERIALIZE_EXPENSIVE_ASSERT(cellId2FaceId.size() == ra.cellCount());
		for(cellid_type cellId(0), s(cellId2FaceId.size()); cellId != s; ++cellId) {
			SSERIALIZE_EXPENSIVE_ASSERT(FaceId(cellId2FaceId.at(cellId)) == ra.faceIdFromCellId(cellId));
		}
	}
#endif
	return dest;
}

bool OsmTriangulationRegionStore::selfTest() {
	if (m_cs & CS_HAVE_CELLS) {
		std::unordered_set<cellid_type> cellIds;
		for(All_faces_iterator it(m_grid.tds().all_faces_begin()), end(m_grid.tds().all_faces_end()); it != end; ++it) {
			if (!it->info().hasCellId()) {
				return false;
			}
			else if (!m_grid.tds().is_infinite(it)) {
				cellid_type cellId = it->info().cellId();
				cellIds.insert(cellId);
			}
			else {
				if (it->info().cellId() != InfiniteFacesCellId) {
					return false;
				}
			}
		}
		
		//now check for missing cellIds, skip cellId=0 since there are not neccessarily faces that are not in any region
		bool allOk = true;
		for(cellid_type i(1), s(sserialize::narrow_check<cellid_type>(cellIds.size())); i < s; ++i) {
			if (!cellIds.count(i)) {
				std::cout << "OsmTriangulationRegionStore::selfTest: missing cellId=" << i << " out of " << cellIds.size() <<'\n';
				allOk = false;
			}
		}
		std::cout << std::flush;
		if (!allOk) {
			return false;
		}
		
		for(const Face_handle & fh : m_grid.grid().storage()) {
			if (!fh->info().hasCellId()) {
				return false;
			}
		}
	}
	
	if (m_isConnected) {
		std::vector<Face_handle> cellRepresentatives;
		std::vector<cellsize_type> cellSizes;
		cellInfo(cellRepresentatives, cellSizes);
		SSERIALIZE_CHEAP_ASSERT(cellRepresentatives.size() == cellSizes.size() && cellSizes.size() == cellCount());
		CTGraph cg;
		for(cellid_type cellId(1), s(sserialize::narrow_check<cellid_type>(cellRepresentatives.size())); cellId < s; ++cellId) {
			ctGraph(cellRepresentatives.at(cellId), cg);
			if (cg.size() != cellSizes.at(cellId)) {
				return false;
			}
		}
	}
	return true;
}


}//end namespace
