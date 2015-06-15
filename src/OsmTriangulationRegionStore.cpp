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

void OsmTriangulationRegionStore::clearRefinement() {
	for(Finite_faces_iterator it(finite_faces_begin()), end(finite_faces_end()); it != end; ++it) {
		if (m_faceToCellId.is_defined(it)) {
			m_faceToCellId[it] = m_refinedCellIdToUnrefined.at(m_faceToCellId[it]);
		}
	}
	m_refinedCellIdToUnrefined.clear();
	for(uint32_t i(0), s(m_cellIdToCellList.size()); i < s; ++i) {
		m_refinedCellIdToUnrefined.push_back(i);
	}
}

void OsmTriangulationRegionStore::refineCells(uint32_t cellSizeTh, uint32_t threadCount) {
	//first clear the old refinement
	for(Finite_faces_iterator it(finite_faces_begin()), end(finite_faces_end()); it != end; ++it) {
		if (m_faceToCellId.is_defined(it)) {
			m_faceToCellId[it] = m_refinedCellIdToUnrefined.at(m_faceToCellId[it]);
		}
	}
	m_refinedCellIdToUnrefined.clear();
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
		std::vector<uint32_t> cellSizes;
		std::vector<Face_handle> cellRep;
		cellInfo(cellRep, cellSizes);

		//Stuff needed to handle the explicit dual-graph
		CTGraph cg;
		std::vector< std::pair<uint32_t, uint32_t> > bfsTree;
		std::vector<uint32_t> leftNodeBfs, rightNodeBfs; //(id, hopdist)
		SimpleBitVector processedTreeNodes;
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
				
				ctGraph(cellRep.at(cellId), cg);
				assert(cg.size() == cellSizes.at(cellId));
				
				processedTreeNodes.resize(cg.size());
				processedTreeNodes.reset();
				
				//get the face with maximum hop-distance, this is currently O(n^2)
				cg.calcMaxHopDistance(bfsTree);
				
				uint32_t leftNode = bfsTree.front().first;
				uint32_t rightNode = bfsTree.back().first;
				
				uint32_t newCellId = m_refinedCellIdToUnrefined.size();
				m_refinedCellIdToUnrefined.push_back(m_refinedCellIdToUnrefined.at(cellId));
				
				//set new cellReps
				cellRep.at(cellId) = cg.face(leftNode);
				cellRep.push_back(cg.face(rightNode));
				
				//set cellid of right node
				m_faceToCellId[cg.face(rightNode)] = newCellId;
				
				//now assign the nodes their ids
				processedTreeNodes.set(leftNode);
				processedTreeNodes.set(rightNode);
				
				leftNodeBfs.clear();
				rightNodeBfs.clear();
				
				leftNodeBfs.emplace_back(leftNode);
				rightNodeBfs.emplace_back(rightNode);
				
				uint32_t oldCellSize = 1;
				uint32_t newCellSize = 1;
				uint32_t processedNodesCount = 2;
				
				for(uint32_t s(cellSizes.at(cellId)), leftNodeI(0), rightNodeI(0); processedNodesCount < s; ) {
					//iterate over all neighbors and set their hop-distances and their repective cellid
					//we only need to set new cellIds for the right-neighbor-cells
					if (leftNodeI < leftNodeBfs.size()) {
						uint32_t lfbnId = leftNodeBfs.at(leftNodeI);
						const CTGraph::FaceNode & lfn = cg.node(lfbnId);
						for(int i(0); i < 3; ++i) {
							uint32_t nid = lfn.neighbours[i];
							if (nid != CTGraph::FaceNode::NullNeighbor && !processedTreeNodes.isSet(nid)) {
								++processedNodesCount;
								++oldCellSize;
								processedTreeNodes.set(nid);
								leftNodeBfs.emplace_back(nid);
							}
						}
						++leftNodeI;
					}
					if (rightNodeI < rightNodeBfs.size()) {
						uint32_t rfbnId = rightNodeBfs.at(rightNodeI);
						const CTGraph::FaceNode & rfn = cg.node(rfbnId);
						for(int i(0); i < 3; ++i) {
							uint32_t nid = rfn.neighbours[i];
							if (nid != CTGraph::FaceNode::NullNeighbor && !processedTreeNodes.isSet(nid)) {
								++processedNodesCount;
								++newCellSize;
								processedTreeNodes.set(nid);
								rightNodeBfs.emplace_back(nid);
								m_faceToCellId[cg.face(nid)] = newCellId;
							}
						}
						++rightNodeI;
					}
				}
				assert(newCellSize + oldCellSize == cellSizes.at(cellId));
				cellSizes.at(cellId) = oldCellSize;
				cellSizes.push_back(newCellSize);
				
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
	
	out << "Cell Triangle stats: \n";
	out << "\tmin: " << *minElem << " at " << minElem - triangCountOfCells.begin() << "\n";
	out << "\tmax: " << *maxElem << " at " << maxElem - triangCountOfCells.begin() << "\n";
	out << "\tmedian: " << sserialize::statistics::median(triangCountOfCells.begin()+1, triangCountOfCells.end(), 0) << "\n";
	out << "\tmean: " << sserialize::statistics::mean(triangCountOfCells.begin()+1, triangCountOfCells.end(), 0) << "\n";
	out << std::flush;
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