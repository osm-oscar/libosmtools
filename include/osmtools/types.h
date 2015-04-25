#ifndef OSM_TOOLS_TYPES_H
#define OSM_TOOLS_TYPES_H
#include <sserialize/containers/CFLArray.h>
#include <sserialize/spatial/GeoPolygon.h>
#include <sserialize/spatial/GeoMultiPolygon.h>

namespace osmtools {

typedef sserialize::CFLArray<sserialize::spatial::GeoPoint, sserialize::detail::CFLArray::DefaultPointerGetter<sserialize::spatial::GeoPoint, false> > PolygonPointsContainer;

typedef sserialize::spatial::detail::GeoPolygon<PolygonPointsContainer> OsmGeoPolygon;

typedef sserialize::CFLArray<OsmGeoPolygon, sserialize::detail::CFLArray::DefaultPointerGetter<OsmGeoPolygon, false> > OsmGeoPolygonsContainer;

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