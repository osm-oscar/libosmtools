#ifndef LIBOSMTOOLS_TRIANGLUATION_GRID_LOCATOR_H
#define LIBOSMTOOLS_TRIANGLUATION_GRID_LOCATOR_H
#include <sserialize/spatial/RWGeoGrid.h>
#include <sserialize/utility/utilmath.h>
#include <sserialize/Static/Triangulation.h>
#include <CGAL/number_utils.h>
#include <CGAL/Unique_hash_map.h>

namespace osmtools {

template<typename TDs>
class GridLocator {
public:
	typedef TDs TriangulationDataStructure;
	typedef typename TDs::Face_handle Face_handle;
	typedef typename TDs::Vertex_handle Vertex_handle;
	typedef sserialize::spatial::RWGeoGrid<Face_handle> Grid;
private:
	typedef typename TriangulationDataStructure::Triangulation::Geom_traits::Point_2 Point_2;
private:
	TriangulationDataStructure m_tds;
	sserialize::spatial::RWGeoGrid<Face_handle> m_grid;
public:
	GridLocator() {}
	void initGrid(uint32_t latCount, uint32_t lonCount);
	inline TriangulationDataStructure & tds() { return m_tds; }
	inline const TriangulationDataStructure & tds() const { return m_tds; }
	inline Grid & grid() { return m_grid; }
	inline const Grid & grid() const { return m_grid; }
	///@thread-safety NO
	Face_handle locate(double x, double y) const;
	inline bool contains(double lat, double lon) { return m_grid.contains(lat, lon); }
	///serialize this to sserialize::Static::spatial::TriangulationGridLocator
	sserialize::UByteArrayAdapter & append(sserialize::UByteArrayAdapter& dest, CGAL::Unique_hash_map< GridLocator<TDs>::Face_handle, uint32_t > & face2FaceId) const;
};


template<typename TDs>
void
GridLocator<TDs>::initGrid(uint32_t latCount, uint32_t lonCount) {
	assert(tds().number_of_vertices());
	if (!tds().number_of_vertices()) {
		return;
	}
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
template<typename TDs>
sserialize::UByteArrayAdapter&
GridLocator<TDs>::append(sserialize::UByteArrayAdapter& dest, CGAL::Unique_hash_map< GridLocator<TDs>::Face_handle, uint32_t > & face2FaceId) const {
	dest.putUint8(1);//version
	CGAL::Unique_hash_map<Vertex_handle, uint32_t> vertex2VertexId;
	sserialize::Static::spatial::Triangulation::append(m_tds, face2FaceId, vertex2VertexId, dest);
	vertex2VertexId.clear();
	
	dest << static_cast<const sserialize::spatial::GeoGrid&>(m_grid);
	sserialize::Static::ArrayCreator<uint32_t> ac(dest);
	for(uint32_t i(0), s(m_grid.tileCount()); i < s; ++i) {
		const Face_handle & fh = m_grid.binAt(i);
		assert(face2FaceId.is_defined(fh));
		ac.put(face2FaceId[fh]);
	}
	ac.flush();
	return dest;
}



}//end namespace

#endif