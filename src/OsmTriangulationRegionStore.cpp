#include <osmtools/OsmTriangulationRegionStore.h>

namespace osmtools {

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


}//end namespace