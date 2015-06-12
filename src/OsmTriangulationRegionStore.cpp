#include <osmtools/OsmTriangulationRegionStore.h>
#include <sserialize/utility/ThreadPool.h>

namespace osmtools {
namespace detail {
namespace OsmTriangulationRegionStore {

void CTGraphBase::calcMaxHopDistance(std::vector< std::pair<uint32_t, uint32_t> > & bfsTree) {
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
	
	w.calc(wctx.maxHopDistRoot);
	bfsTree = std::move(w.bfsTree);
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
	
	

	cgFaces.push_back(rfh);
	faceToNodeId[rfh] = 0;
	for(uint32_t i(0); i < cgFaces.size(); ++i) {
		Face_handle fh = cgFaces[i];
		for(int j(0); j < 3; ++j) {
			Face_handle nfh = fh->neighbor(j);
			if (m_faceToCellId.is_defined(nfh) && m_faceToCellId[nfh] == myCellId && !faceToNodeId.is_defined(nfh)) {
				faceToNodeId[nfh] = cgFaces.size();
				cgFaces.emplace_back(nfh);
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
}

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
	
	std::cout << "Cell Triangle stats: \n";
	std::cout << "\tmin: " << *minElem << " at " << minElem - triangCountOfCells.begin() << "\n";
	std::cout << "\tmax: " << *maxElem << " at " << maxElem - triangCountOfCells.begin() << "\n";
	std::cout << "\tmedian: " << sserialize::statistics::median(triangCountOfCells.begin()+1, triangCountOfCells.end(), 0) << "\n";
	std::cout << "\tmean: " << sserialize::statistics::mean(triangCountOfCells.begin()+1, triangCountOfCells.end(), 0) << "\n";
	std::cout << std::flush;
}


}//end namespace