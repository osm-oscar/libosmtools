#ifndef LIBOSMTOOLS_TRIANGLUATION_GRID_LOCATOR_H
#define LIBOSMTOOLS_TRIANGLUATION_GRID_LOCATOR_H
#include <sserialize/spatial/RWGeoGrid.h>
#include <sserialize/algorithm/utilmath.h>
#include <sserialize/Static/Triangulation.h>
#include <sserialize/Static/TriangulationGridLocator.h>
#include <CGAL/number_utils.h>
#include <CGAL/Unique_hash_map.h>

namespace osmtools {

template<typename TDs, bool TNumberTypeIsThreadSafe=false>
class GridLocator {
public:
	typedef TDs TriangulationDataStructure;
	typedef typename TDs::Face_handle Face_handle;
	typedef typename TDs::Vertex_handle Vertex_handle;
	typedef sserialize::spatial::RWGeoGrid<Face_handle> Grid;
public:
	GridLocator() {}
	GridLocator(GridLocator && other) = default;
	GridLocator & operator=(GridLocator && other) {
		m_tds = std::move(other.m_tds);
		m_grid = std::move(other.m_grid);
		return *this;
	}
	void initGrid(uint32_t latCount, uint32_t lonCount);
	inline TriangulationDataStructure & tds() { return m_tds; }
	inline const TriangulationDataStructure & tds() const { return m_tds; }
	inline Grid & grid() { return m_grid; }
	inline const Grid & grid() const { return m_grid; }
	///@thread-safety YES
	Face_handle locate(double x, double y) const;
	inline bool contains(double lat, double lon) { return m_grid.contains(lat, lon); }
	///serialize this to sserialize::Static::spatial::TriangulationGridLocator
	///@thread-safety NO
	sserialize::UByteArrayAdapter & append( sserialize::UByteArrayAdapter& dest,
											CGAL::Unique_hash_map<Face_handle, uint32_t >& face2FaceId,
											sserialize::Static::spatial::Triangulation::GeometryCleanType gct
	);
private:
	typedef typename TriangulationDataStructure::Triangulation::Geom_traits::Point_2 Point_2;
private:
	TriangulationDataStructure m_tds;
	sserialize::spatial::RWGeoGrid<Face_handle> m_grid;
	mutable std::mutex m_lock;
};


template<typename TDs, bool TNumberTypeIsThreadSafe>
void
GridLocator<TDs, TNumberTypeIsThreadSafe>::initGrid(uint32_t latCount, uint32_t lonCount) {
	SSERIALIZE_CHEAP_ASSERT(tds().number_of_vertices());
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
	std::vector<std::pair<double, double>> cellPts;
	for(uint32_t lat(0); lat < latCount; ++lat) {
		for(uint32_t lon(0); lon < lonCount; ++lon) {
			sserialize::spatial::GeoRect cellRect = m_grid.cellBoundary(lat, lon);
			double midLat = cellRect.midLat();
			double midLon = cellRect.midLon();
			fh = m_tds.locate(Point_2(midLat, midLon), fh);
			if (m_tds.is_infinite(fh)) {
				//face is infinite, this either means that there are no triangle within this cell or 
				//the midpoint of the cell is outside the triangulation
				//we have to explore the neighborhood of this cell
				//or the fast way: try the end points
				cellPts.emplace_back(cellRect.maxLat(), cellRect.maxLon());
				cellPts.emplace_back(cellRect.maxLat(), cellRect.minLon());
				cellPts.emplace_back(cellRect.minLat(), cellRect.maxLon());
				cellPts.emplace_back(cellRect.minLat(), cellRect.minLon());
				for(const auto & x : cellPts) {
					fh = m_tds.locate(Point_2(x.first, x.second), fh);
					if (!m_tds.is_infinite(fh)) {
						break;
					}
				}
				cellPts.clear();
			}
			m_grid.at(midLat, midLon) = fh;
		}
	}
	//now check all cells for valid entries.
	//if a cell has an infinite face as entry, try to use the neighbors face
	for(uint32_t lat(0); lat < latCount; ++lat) {
		for(uint32_t lon(0); lon < lonCount; ++lon) {
			sserialize::spatial::GeoGrid::GridBin bin = m_grid.select( m_grid.selectBin(lat, lon) );
			if (!m_tds.is_infinite(m_grid.at(bin.tile))) {
				continue;
			}
			bool notFound = false;
			for(int i(-1); i <= 1 && notFound; ++i) {
				for(int j(-1); j <= 1 && notFound; ++j) {
					int32_t x = (int32_t)bin.x+i;
					int32_t y = (int32_t)bin.y+j;
					if (x >= 0 && y >= 0) {
						fh = m_grid.at((uint32_t)x, (uint32_t)y);
						if (!m_tds.is_infinite(fh)) {
							m_grid.at(bin.tile) = fh;
							notFound = false;
						}
					}
				}
			}
		}
	}
}

template<typename TDs, bool TNumberTypeIsThreadSafe>
typename GridLocator<TDs, TNumberTypeIsThreadSafe>::Face_handle
GridLocator<TDs, TNumberTypeIsThreadSafe>::locate(double x, double y) const {
	if (!m_grid.contains(x,y)) {
		return Face_handle();
	}
	const Face_handle & hint = m_grid.at(x,y);
	Point_2 p(x,y);
	
	if (TNumberTypeIsThreadSafe) {
		return m_tds.locate(p, hint);
	}
	else {
		std::lock_guard<std::mutex> lck(m_lock);
		return m_tds.locate(p, hint);
	}
}
template<typename TDs, bool TNumberTypeIsThreadSafe>
sserialize::UByteArrayAdapter&
GridLocator<TDs, TNumberTypeIsThreadSafe>::append(sserialize::UByteArrayAdapter& dest, CGAL::Unique_hash_map<Face_handle, uint32_t > & face2FaceId, sserialize::Static::spatial::Triangulation::GeometryCleanType gct) {
#ifdef SSERIALIZE_EXPENSIVE_ASSERT_ENABLED
	sserialize::UByteArrayAdapter::OffsetType initialPutPtr = dest.tellPutPtr();
#endif
	dest.putUint8(1);//version
	{
		CGAL::Unique_hash_map<Vertex_handle, uint32_t> vertex2VertexId;
		sserialize::Static::spatial::Triangulation::append(m_tds, face2FaceId, vertex2VertexId, dest, gct);
	}
	
	dest << static_cast<const sserialize::spatial::GeoGrid&>(m_grid);
	sserialize::Static::ArrayCreator<uint32_t> ac(dest);
	for(uint32_t i(0), s(m_grid.tileCount()); i < s; ++i) {
		const Face_handle & fh = m_grid.binAt(i);
		SSERIALIZE_CHEAP_ASSERT(face2FaceId.is_defined(fh) || m_tds.is_infinite(fh));
		uint32_t faceId = sserialize::Static::spatial::Triangulation::NullFace;
		if (face2FaceId.is_defined(fh)) {
			faceId = face2FaceId[fh];
		}
		ac.put(faceId);
	}
	ac.flush();
	
#ifdef SSERIALIZE_EXPENSIVE_ASSERT_ENABLED
	{
		sserialize::UByteArrayAdapter tmp(dest);
		tmp.setPutPtr(initialPutPtr);
		tmp.shrinkToPutPtr();
		sserialize::Static::spatial::TriangulationGridLocator str(tmp);
		SSERIALIZE_EXPENSIVE_ASSERT_EQUAL(str.grid().tileCount(), m_grid.tileCount());
		for(uint32_t i(0), s(m_grid.tileCount()); i < s; ++i) {
			Face_handle fh = m_grid.at(i);
			uint32_t myFaceId = sserialize::Static::spatial::Triangulation::NullFace;
			if (face2FaceId.is_defined(fh)) {
				myFaceId = face2FaceId[fh];
			}
			SSERIALIZE_EXPENSIVE_ASSERT_EQUAL(str.grid().at(i), myFaceId);
			SSERIALIZE_EXPENSIVE_ASSERT(face2FaceId.is_defined(fh) || m_tds.is_infinite(fh));
		}
	}
#endif
	return dest;
}



}//end namespace

#endif