#include <osmtools/OsmTriangulationRegionStore.h>

namespace osmtools {

OsmTriangulationRegionStore::Point OsmTriangulationRegionStore::centroid(const OsmTriangulationRegionStore::Face_handle& fh) {
	return CGAL::centroid(fh->vertex(0)->point(), fh->vertex(1)->point(), fh->vertex(2)->point());
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