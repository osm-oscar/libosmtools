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
		ET_NONE=0, ET_BOUNDARIES=0x1, ET_LANDUSE=0x2, ET_AREA=0x4, ET_BUILDING=0x4|0x8,
		ET_ALL_BUT_BUILDINGS=ET_BOUNDARIES|ET_LANDUSE|ET_AREA, ET_ALL=ET_ALL_BUT_BUILDINGS|ET_BUILDING} ExtractionTypes;
	generics::RCPtr<osmpbf::AbstractTagFilter> createExtractionFilter(ExtractionTypes extractionTypes, bool needsName);
};



}}}//end namespace osmtools

#endif