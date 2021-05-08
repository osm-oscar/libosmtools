#include <osmtools/CellCriteria.h>

namespace osmtools {
namespace CellCriteria {

CellTriangleCountCriteria::CellTriangleCountCriteria(uint32_t cellSizeTh) :
m_cellSizeTh(cellSizeTh)
{}

int CellTriangleCountCriteria::dataDependence() const {
	return DD_CELLSIZES;
}

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

int CellDiagonalCriteria::dataDependence() const {
	return  DD_CELL_GRAPH | DD_NEW_FACE_CELL_IDS;
}

bool CellDiagonalCriteria::init(const ::osmtools::OsmTriangulationRegionStore & store) {
	return store.grid().grid().rect().diagInM() > m_maxCellDiameter;
}

void CellDiagonalCriteria::begin() {
	std::cout << "Splitting cells larger than " << m_maxCellDiameter << "m in diameter" << std::endl;
}
void CellDiagonalCriteria::end() {}

bool CellDiagonalCriteria::refine(const State & state) {
	struct Data {
		std::size_t triangCount = 0;
		sserialize::spatial::GeoRect bound;
	};
	std::map<cellid_type, Data> data;
	for(std::size_t i(0), s(state.cg.size()); i < s; ++i) {
		cellid_type cellId = state.newFaceCellIds[i];
		Data & d = data[cellId];
		auto fh = state.cg.face(FaceId(i));
		for(int i(0); i < 3; ++i) {
			const auto & p = fh->vertex(i)->point();
			double lat = CGAL::to_double(p.x());
			double lon = CGAL::to_double(p.y());
			d.bound.enlarge(lat, lon);
			d.triangCount += 1;
		}
	}
	
	for(const auto & d : data) {
		if (d.second.triangCount > 1 && d.second.bound.diagInM() > m_maxCellDiameter) {
			return true;
		}
	}
	return false;
}

bool CellDiagonalCriteria::refine(cellid_type cellId, const State & state) {
	sserialize::spatial::GeoRect bound;
	std::size_t triangCount = 0;
	for(std::size_t i(0), s(state.cg.size()); i < s; ++i) {
		if (cellId == state.newFaceCellIds.at(i)) {
			auto fh = state.cg.face(FaceId(i));
			for(int i(0); i < 3; ++i) {
				const auto & p = fh->vertex(i)->point();
				double lat = CGAL::to_double(p.x());
				double lon = CGAL::to_double(p.y());
				bound.enlarge(lat, lon);
				++triangCount;
			}
		}
	}
	return triangCount > 1 && bound.diagInM() > m_maxCellDiameter;
}

CellDiagonalCriteria::CellCriteriaInterface * CellDiagonalCriteria::copy() {
	return new CellDiagonalCriteria(m_maxCellDiameter);
}

}} //end namespace osmtools::CellCriteria
