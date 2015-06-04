#include <osmtools/types.h>

namespace sserialize {
namespace spatial {
namespace detail {

template<>
sserialize::spatial::detail::GeoPolygon<osmtools::PolygonPointsContainer>
sserialize::spatial::detail::GeoPolygon<osmtools::PolygonPointsContainer>::fromRect(const sserialize::spatial::GeoRect & rect) {
	osmtools::PolygonPointsContainer c(5);
	c[0] = GeoPoint(rect.minLat(), rect.minLon());
	c[1] = GeoPoint(rect.minLat(), rect.maxLon());
	c[2] = GeoPoint(rect.maxLat(), rect.maxLon());
	c[3] = GeoPoint(rect.maxLat(), rect.minLon());
	c[4] = c[0];
	return GeoPolygon(rect, std::move(c));
}

}}}//end namespace sserialize::spatial::detail