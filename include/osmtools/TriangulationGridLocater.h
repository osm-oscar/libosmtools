#ifndef LIBOSMTOOLS_TRIANGLUATION_GRID_LOCATOR_H
#define LIBOSMTOOLS_TRIANGLUATION_GRID_LOCATOR_H
#include <sserialize/spatial/RWGeoGrid.h>
#include <sserialize/utility/utilmath.h>
#include <CGAL/number_utils.h>

namespace osmtools {

template<typename TDs>
class GridLocator {
public:
	typedef TDs TriangulationDataStructure;
	typedef typename TDs::Face_handle Face_handle;
private:
	typedef typename TDs::Triangulation::Geom_traits::Point_2 Point_2;
private:
	TDs m_tds;
	sserialize::spatial::RWGeoGrid<Face_handle> m_grid;
public:
	GridLocator() {}
	void initGrid(uint32_t latCount, uint32_t lonCount);
	TDs & tds() { return m_tds; }
	const TDs & tds() const { return m_tds; }
	Face_handle locate(double lat, double lon) const;
};


template<typename TDs>
void
GridLocator<TDs>::initGrid(uint32_t latCount, uint32_t lonCount) {
	{
		sserialize::AtomicMinMax<double> lat, lon;
		for(auto it(m_tds.finite_vertices_begin()), end(m_tds.finite_vertices_end()); it != end; ++it) {
			auto p = it->point();
			lat.update(CGAL::to_double(p.x()));
			lon.update(CGAL::to_double(p.y()));
		}
		m_grid = decltype(m_grid)(sserialize::spatial::GeoRect(lat.min(), lat.max(), lon.min(), lon.max()), latCount, lonCount);
	}
	Face_handle fh;
	for(uint32_t lat(0); lat < latCount; ++lat) {
		for(uint32_t lon(0); lon < lonCount; ++lon) {
			sserialize::spatial::GeoRect cellRect = m_grid.cellBoundary(lat, lon);
			double midLat = 0.5*(cellRect.maxLat()+cellRect.minLat());
			double midLon = 0.5*(cellRect.maxLon()+cellRect.minLon());
			fh = m_tds.locate(Point_2(midLat, midLon), fh);
			m_grid.at(midLat, midLon) = fh;
		}
	}
}

template<typename TDs>
typename GridLocator<TDs>::Face_handle
GridLocator<TDs>::locate(double x, double y) const {
	if (!m_grid.contains(x,y)) {
		return Face_handle();
	}
	return m_tds.locate(Point_2(x,y), m_grid.at(x,y));
}


}//end namespace

#endif