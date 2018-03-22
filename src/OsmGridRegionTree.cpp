#include <osmtools/OsmGridRegionTree.h>

namespace osmtools {

OsmGridRegionTreeBase::FixedSizeDiagRefiner::FixedSizeDiagRefiner(double minDiag, uint32_t latCount, uint32_t lonCount) :
m_dc( std::shared_ptr<sserialize::spatial::detail::DistanceCalculator>(new sserialize::spatial::detail::GeodesicDistanceCalculator()) ),
m_minDiag(minDiag),
m_latCount(latCount),
m_lonCout(lonCount)
{}
OsmGridRegionTreeBase::FixedSizeDiagRefiner::~FixedSizeDiagRefiner() {}

bool OsmGridRegionTreeBase::FixedSizeDiagRefiner::operator()(
	const sserialize::spatial::GeoRect& maxBounds,
	const std::vector<sserialize::spatial::GeoRegion*>& rId2Ptr,
	const std::vector<uint32_t> & sortedRegions,
	sserialize::spatial::GeoGrid & newGrid) const
{
	if (sortedRegions.size() && m_dc.calc(maxBounds.minLat(), maxBounds.minLon(), maxBounds.maxLat(), maxBounds.maxLon()) > m_minDiag) {
		std::vector< uint32_t >::const_iterator it(sortedRegions.cbegin());
		sserialize::spatial::GeoRect myBounds(rId2Ptr[*it]->boundary());
		++it;
		for(std::vector<uint32_t>::const_iterator end(sortedRegions.cend()); it != end; ++it) {
			myBounds.enlarge(rId2Ptr[*it]->boundary());
		}
		newGrid = sserialize::spatial::GeoGrid( myBounds / maxBounds, m_latCount, m_lonCout);
		return true;
	}
	return false;
}

//BEGIN: OsmGridRegionTreeBase

void OsmGridRegionTreeBase::push_back(const sserialize::spatial::GeoRegion * p) {
	sserialize::spatial::GeoRegion * r = 0;
	switch (p->type()) {
	case sserialize::spatial::GS_POLYGON:
		if (dynamic_cast<const sserialize::spatial::GeoPolygon*>(p)) {
			r = ConvertGP<sserialize::spatial::GeoPolygon>::conv(m_polygonPoints, p);
		}
		else if (dynamic_cast<const osmtools::OsmGeoPolygon*>(p)) {
			r = ConvertGP<osmtools::OsmGeoPolygon>::conv(m_polygonPoints, p);
		}
		else {
			throw sserialize::TypeMissMatchException("OsmGridRegionTree");
		}
		break;
	case sserialize::spatial::GS_MULTI_POLYGON:
		if (dynamic_cast<const sserialize::spatial::GeoMultiPolygon*>(p)) {
			r = ConvertGMP<sserialize::spatial::GeoMultiPolygon>::conv(m_polygonPoints, m_polygonsContainer, p);
		}
		else if (dynamic_cast<const osmtools::OsmGeoMultiPolygon*>(p)) {
			r = ConvertGMP<osmtools::OsmGeoMultiPolygon>::conv(m_polygonPoints, m_polygonsContainer, p);
		}
		else {
			throw sserialize::TypeMissMatchException("OsmGridRegionTree");
		}
		break;
	default:
		throw sserialize::TypeMissMatchException("OsmGridRegionTree");
		break;
	}
	SSERIALIZE_CHEAP_ASSERT_EQUAL(r->size(), p->size());
	r->recalculateBoundary();
	SSERIALIZE_CHEAP_ASSERT_EQUAL(r->size(), p->size());
	m_regions.push_back(r);
}


OsmGridRegionTreeBase::OsmGridRegionTreeBase() : m_polygonPoints(sserialize::MM_SHARED_MEMORY), m_polygonsContainer(sserialize::MM_SHARED_MEMORY), m_latRefineCount(2), m_lonRefineCount(2), m_refineMinDiag(250) {}

OsmGridRegionTreeBase::~OsmGridRegionTreeBase() {
	for(std::vector<sserialize::spatial::GeoRegion*>::iterator it(m_regions.begin()), end(m_regions.end()); it != end; ++it) {
		delete *it;
	}
}

const sserialize::spatial::GridRegionTree & OsmGridRegionTreeBase::grt() const {
	return m_grt;
}

void OsmGridRegionTreeBase::clearGRT() {
	m_grt = sserialize::spatial::GridRegionTree();
}

void OsmGridRegionTreeBase::clear() {
	clearGRT();
	m_polygonPoints = PointDataContainer();
	m_polygonsContainer = PolygonsContainer();
	for(std::vector<sserialize::spatial::GeoRegion*>::iterator it(m_regions.begin()), end(m_regions.end()); it != end; ++it) {
		delete *it;
	}
	m_regions = RegionsContainer();
}

void OsmGridRegionTreeBase::snapPoints() {
	for(Point & p : m_polygonPoints) {
		p.snap();
	}
	for(auto & r : regions()) {
		r->recalculateBoundary();
	}
}

void OsmGridRegionTreeBase::printStats(std::ostream & out) {
	out << "OsmGridRegionTree::printStats--BEGIN\n";
	m_grt.printStats(out);
	std::cout << "#points: " << m_polygonPoints.size() << "=" << sserialize::prettyFormatSize(m_polygonPoints.size()*sizeof(PolygonPointsContainer::value_type)) << "\n";
	std::cout << "#GeoMultiPolygons: " << m_polygonsContainer.size() << "=" << sserialize::prettyFormatSize(m_polygonsContainer.size()*sizeof(PolygonsContainer::value_type)) << "\n";
	out << "OsmGridRegionTree::printStats--END\n";
}

std::size_t OsmGridRegionTreeBase::size() const {
	return regions().size();
}

const OsmGridRegionTreeBase::RegionsContainer &
OsmGridRegionTreeBase::regions() const {
	return m_regions;
}

void OsmGridRegionTreeBase::setRefinerOptions(uint32_t latRefineCount, uint32_t lonRefineCount, double minDiag) {
	m_latRefineCount = latRefineCount;
	m_lonRefineCount = lonRefineCount;
	m_refineMinDiag = minDiag;
}
///you can still add more regions after calling this but they will not be part of the tree
///if you want them in the tree aswell then you should call this again (which will rebuild the tree)
void OsmGridRegionTreeBase::addPolygonsToRaster(unsigned int gridLatCount, unsigned int gridLonCount) {
	FixedSizeDiagRefiner refiner(m_refineMinDiag, m_latRefineCount, m_lonRefineCount);
	sserialize::spatial::GeoRect initRect( sserialize::spatial::GeoShape::bounds(m_regions.cbegin(), m_regions.cend()) );
	sserialize::spatial::GeoGrid initGrid(initRect, gridLatCount, gridLonCount);
	typedef sserialize::spatial::GridRegionTree::TypeTraits<FixedSizeDiagRefiner, osmtools::OsmGeoPolygon, osmtools::OsmGeoMultiPolygon> MyTypeTraits;
	m_grt = sserialize::spatial::GridRegionTree(initGrid, m_regions.begin(), m_regions.end(), MyTypeTraits(), refiner);
	m_grt.shrink_to_fit();
}


} //end namespace osmtools
