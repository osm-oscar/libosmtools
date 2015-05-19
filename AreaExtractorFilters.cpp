#include <osmtools/AreaExtractorFilters.h>
#include <vector>
#include <osmpbf/iprimitive.h>

namespace osmtools {
namespace detail {
namespace AreaExtractor {

generics::RCPtr<osmpbf::AbstractTagFilter> Base::createExtractionFilter(ExtractionTypes extractionTypes) {
	if (! (extractionTypes & (~(ET_PRIMITIVE_RELATIONS|ET_PRIMITIVE_WAYS)))) { //extract all
		
	}
	generics::RCPtr<osmpbf::AbstractTagFilter> mainFilter;
	generics::RCPtr<osmpbf::OrTagFilter> areaFilter(new osmpbf::OrTagFilter());

	generics::RCPtr<osmpbf::AndTagFilter> wayFilter(new osmpbf::AndTagFilter());
	generics::RCPtr<osmpbf::AndTagFilter> relationFilter(new osmpbf::AndTagFilter());
	
	//setup primitive type filters
	wayFilter->addChild(new osmpbf::PrimitiveTypeFilter(osmpbf::WayPrimitive));
	relationFilter->addChild(new osmpbf::PrimitiveTypeFilter(osmpbf::RelationPrimitive));
	
	if (extractionTypes & ET_BUILDING) {
		areaFilter->addChild(new osmpbf::KeyOnlyTagFilter("building"));
	}
	if (extractionTypes & ET_BOUNDARIES) {
		areaFilter->addChild(new osmpbf::KeyOnlyTagFilter("boundary"));
	}
	if (extractionTypes & ET_LANDUSE) {
		areaFilter->addChild(new osmpbf::KeyOnlyTagFilter("landuse"));
	}
	if (extractionTypes & ET_NATURAL) {
		areaFilter->addChild(new osmpbf::KeyOnlyTagFilter("natural"));
	}
	if (extractionTypes & ET_AREA) {
		osmpbf::AbstractTagFilter * areaTagFilter = 0;
		std::vector<osmpbf::AbstractTagFilter*> areaExclusions;
		if ((extractionTypes & ET_BUILDING) != ET_BUILDING) {
			areaExclusions.push_back(new osmpbf::BoolTagFilter("building", false));
		}
		if ((extractionTypes & ET_BOUNDARIES) != ET_BOUNDARIES) {
			areaExclusions.push_back(new osmpbf::BoolTagFilter("boundary", false));
		}
		if ((extractionTypes & ET_LANDUSE) != ET_LANDUSE) {
			areaExclusions.push_back(new osmpbf::BoolTagFilter("landuse", false));
		}
		if ((extractionTypes & ET_NATURAL) != ET_NATURAL) {
			areaExclusions.push_back(new osmpbf::BoolTagFilter("natural", false));
		}
		if (areaExclusions.size()) {
			osmpbf::AndTagFilter * andTagFilter = new osmpbf::AndTagFilter();
			andTagFilter->addChild(new osmpbf::KeyOnlyTagFilter("area"));
			andTagFilter->addChildren(areaExclusions.begin(), areaExclusions.end());
			areaTagFilter = andTagFilter;
		}
		else {
			areaTagFilter = new osmpbf::KeyOnlyTagFilter("area");
		}
		areaFilter->addChild(areaTagFilter);
	}
	
	{//set the way filter
		wayFilter->addChild(areaFilter.get());
	}
	
	{//set the relation filter
		osmpbf::AbstractTagFilter * multiPolyFilter;
		if (extractionTypes & ET_ALL_MULTIPOLYGONS) {
			multiPolyFilter = new osmpbf::MultiStringTagFilter("type", {"multipoly", "multipolygon"});
		}
		else {
			multiPolyFilter = osmpbf::newAnd(new osmpbf::MultiStringTagFilter("type", {"multipoly", "multipolygon"}), areaFilter.get());
		}
		if (extractionTypes & ET_BOUNDARIES) {
			relationFilter->addChild( osmpbf::newOr(new osmpbf::StringTagFilter("type", "boundary"), multiPolyFilter) );
		}
		else {
			relationFilter->addChild(multiPolyFilter);
		}
	}
	if ((extractionTypes & ET_PRIMITIVE_WAYS) && extractionTypes & ET_PRIMITIVE_RELATIONS) {
		mainFilter.reset( osmpbf::newOr(wayFilter.get(), relationFilter.get()) );
	}
	else if(extractionTypes & ET_PRIMITIVE_WAYS) {
		mainFilter.reset(wayFilter.get());
	}
	else if (extractionTypes & ET_PRIMITIVE_RELATIONS) {
		mainFilter.reset(relationFilter.get());
	}
	return mainFilter;
}

}}}//end namespace osmtools