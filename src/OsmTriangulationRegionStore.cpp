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
	std::vector<uint32_t> triangCountOfCells(cellCount(), 0);
	for(Finite_faces_iterator it(finite_faces_begin()), end(finite_faces_end()); it != end; ++it) {
		uint32_t fid = cellId(it);
		triangCountOfCells.at(fid) += 1;
	}
	std::sort(triangCountOfCells.begin(), triangCountOfCells.end());
	std::cout << "Cell Triangle stats: \n";
	std::cout << "\tmin: " << triangCountOfCells.front() << "\n";
	std::cout << "\tmax: " << triangCountOfCells.back() << "\n";
	std::cout << "\tmedian: " << triangCountOfCells.at(triangCountOfCells.size()/2) << "\n";
	std::cout << "\tmean: " << sserialize::statistics::mean(triangCountOfCells.begin(), triangCountOfCells.end(), 0) << "\n";
	std::cout << std::flush;
}


}//end namespace