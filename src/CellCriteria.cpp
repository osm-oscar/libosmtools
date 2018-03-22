#include <osmtools/CellCriteria.h>

namespace osmtools {
namespace CellCriteria {

CellTriangleCountCriteria::CellTriangleCountCriteria(uint32_t cellSizeTh) :
m_cellSizeTh(cellSizeTh)
{}

bool CellTriangleCountCriteria::init(const ::osmtools::OsmTriangulationRegionStore & store) {
	return m_cellSizeTh < store.grid().tds().number_of_faces();
}

void CellTriangleCountCriteria::begin() {
	std::cout << "Splitting cells larger than " << m_cellSizeTh << " triangles" << std::endl;
}

void CellTriangleCountCriteria::end() {}

bool CellTriangleCountCriteria::refine(uint32_t cellId, const State & state) {
	return state.cellSizes.at(cellId) > m_cellSizeTh;
}

CellTriangleCountCriteria::CellCriteriaInterface * CellTriangleCountCriteria::copy() {
	return new CellTriangleCountCriteria(m_cellSizeTh);
}

CellDiagonalCriteria::CellDiagonalCriteria(double maxCellDiameter) :
m_maxCellDiameter(maxCellDiameter),
m_dc(sserialize::spatial::DistanceCalculator::DCT_GEODESIC_ACCURATE)
{}

bool CellDiagonalCriteria::init(const ::osmtools::OsmTriangulationRegionStore & store) {
	return store.grid().grid().rect().diagInM() > m_maxCellDiameter;
}

void CellDiagonalCriteria::begin() {
	std::cout << "Splitting cells larger than " << m_maxCellDiameter << "m in diameter" << std::endl;
}
void CellDiagonalCriteria::end() {}

bool CellDiagonalCriteria::refine(const State & state) {
	std::map<uint32_t, sserialize::spatial::GeoRect> bounds;
	for(uint32_t i(0), s(state.cg.size()); i < s; ++i) {
		uint32_t cellId = state.newFaceCellIds[i];
		sserialize::spatial::GeoRect & bound = bounds[cellId];
		auto fh = state.cg.face(i);
		for(int i(0); i < 3; ++i) {
			const auto & p = fh->vertex(i)->point();
			double lat = CGAL::to_double(p.x());
			double lon = CGAL::to_double(p.y());
			bound.enlarge(lat, lon);
		}
	}
	
	for(const auto & bound : bounds) {
		if (bound.second.diagInM() > m_maxCellDiameter) {
			return true;
		}
	}
	return false;
}

bool CellDiagonalCriteria::refine(uint32_t, const State &) {
	throw sserialize::UnimplementedFunctionException("RefineBySize");
	return false;
}

CellDiagonalCriteria::CellCriteriaInterface * CellDiagonalCriteria::copy() {
	return new CellDiagonalCriteria(m_maxCellDiameter);
}

}} //end namespace osmtools::CellCriteria
