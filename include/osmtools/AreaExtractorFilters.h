#ifndef LIBOSMTOOLS_AREA_EXTRACTOR_FILTERS_H
#define LIBOSMTOOLS_AREA_EXTRACTOR_FILTERS_H
#include <osmpbf/filter.h>
/**
  * This is the base of the AreaExtractor. It's only dependency is the osmpbf library
  */

namespace osmtools {
namespace detail {
namespace AreaExtractor {

struct Base {
	typedef enum {
		ET_NONE=0x0,
		ET_PRIMITIVE_WAYS=0x1, ET_PRIMITIVE_RELATIONS=0x2, //the primitive types to extract
		//special regions to extract
		ET_BOUNDARIES=0x4, ET_LANDUSE=0x8, ET_NATURAL=0x10, ET_AREA=0x20, ET_BUILDING=0x20|0x40,
		ET_MULTIPOLYGONS=0x80,
		ET_ALL_SPECIAL_BUT_BUILDINGS=ET_BOUNDARIES|ET_LANDUSE|ET_NATURAL|ET_AREA|ET_PRIMITIVE_WAYS|ET_PRIMITIVE_RELATIONS,
		ET_ALL_SPECIAL=ET_ALL_SPECIAL_BUT_BUILDINGS|ET_BUILDING,
		ET_ALL_MULTIPOLYGONS=ET_MULTIPOLYGONS|ET_PRIMITIVE_RELATIONS,
// 		ET_ALL=ET_ALL_SPECIAL|ET_ALL_MULTIPOLYGONS
	} ExtractionTypes;
	static generics::RCPtr<osmpbf::AbstractTagFilter> createExtractionFilter(ExtractionTypes extractionTypes);
};



}}}//end namespace osmtools

#endif