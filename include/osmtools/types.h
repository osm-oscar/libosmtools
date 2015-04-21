#ifndef OSM_TOOLS_TYPES_H
#define OSM_TOOLS_TYPES_H
#include <sserialize/containers/IFLArray.h>
#include <sserialize/spatial/GeoPolygon.h>
#include <sserialize/spatial/GeoMultiPolygon.h>

namespace osmtools {

typedef sserialize::IFLArray<sserialize::spatial::GeoPoint*, const sserialize::spatial::GeoPoint*> PolygonPointsContainer;

typedef sserialize::spatial::detail::GeoPolygon<PolygonPointsContainer> OsmGeoPolygon;

typedef sserialize::IFLArray<OsmGeoPolygon*, const OsmGeoPolygon*> OsmGeoPolygonsContainer;

typedef sserialize::spatial::detail::GeoMultiPolygon<OsmGeoPolygonsContainer> OsmGeoMultiPolygon;

inline OsmGeoPolygon toOsmGeoPolygon(const sserialize::spatial::GeoPolygon & p) {
	return OsmGeoPolygon(PolygonPointsContainer(p.cbegin(), p.cend()));
}

inline sserialize::spatial::GeoPolygon toGeoPolygon(const OsmGeoPolygon & p) {
	return sserialize::spatial::GeoPolygon(sserialize::spatial::GeoPolygon::PointsContainer(p.cbegin(), p.cend()));
}

}//end namespace osmtools



namespace sserialize {
namespace spatial {
namespace detail {

template<>
sserialize::spatial::detail::GeoPolygon<osmtools::PolygonPointsContainer>
sserialize::spatial::detail::GeoPolygon<osmtools::PolygonPointsContainer>::fromRect(const sserialize::spatial::GeoRect & rect);

}}}//end namespace sserialize::spatial::detail

#endif