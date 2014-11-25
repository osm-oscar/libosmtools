#ifndef OSM_TOOLS_AREA_EXTRACTOR_H
#define OSM_TOOLS_AREA_EXTRACTOR_H
#include <string>
#include <cstdlib>
#include <omp.h>
#include <sserialize/spatial/GeoMultiPolygon.h>
#include <sserialize/spatial/GeoRegionStore.h>
#include <osmpbf/osmfile.h>
#include <osmpbf/primitiveblockinputadaptor.h>
#include <osmpbf/filter.h>
#include <osmpbf/irelation.h>
#include <osmpbf/iway.h>
#include <osmpbf/inode.h>

namespace osmtools {
namespace detail {
namespace AreaExtractor {

struct MultiPolyResolver {
	typedef std::vector<int64_t> RawWay;
	
	struct RawMultiPoly {
		RawWay innerWays;
		RawWay outerWays;
	};
	
	static bool closedPolysFromWays(const std::vector<RawWay> & ways, std::vector<RawWay> & closedWays);
	static bool multiPolyFromWays(const std::vector<RawWay> & innerIn, const std::vector<RawWay> & outerIn, std::vector<RawWay> & innerOut, std::vector<RawWay> & outerOut);
	
	///@param mapper maps from node-id to sserialize::spatial::GeoPoint @return either a GeoMultiPolygon or a GeoPolygon
	template<typename T_ID_GEO_MAPPER>
	static sserialize::spatial::GeoRegion * multiPolyFromClosedWays(const std::vector<RawWay> & inner, const std::vector<RawWay> & outer, const T_ID_GEO_MAPPER & mapper) {
		sserialize::spatial::GeoMultiPolygon * mupoptr = new sserialize::spatial::GeoMultiPolygon();
		sserialize::spatial::GeoMultiPolygon & mupo = *mupoptr;
		for(std::vector<RawWay>::const_iterator polyit(inner.begin()), end(inner.end()); polyit != end; ++polyit) {
			mupo.innerPolygons().resize(mupo.innerPolygons().size()+1);
			sserialize::remap(*polyit, mupo.innerPolygons().back().points(), mapper);
		}
		for(std::vector<RawWay>::const_iterator polyit(outer.begin()), end(outer.end()); polyit != end; ++polyit) {
			mupo.outerPolygons().resize(mupo.outerPolygons().size()+1);
			sserialize::remap(*polyit, mupo.outerPolygons().back().points(), mapper);
		}
		if (mupoptr->innerPolygons().size() == 0 && mupoptr->outerPolygons().size() == 1) {
			using std::swap;
			sserialize::spatial::GeoPolygon * gpo = new sserialize::spatial::GeoPolygon();
			swap(*gpo, mupoptr->outerPolygons().front());
			gpo->recalculateBoundary();
			delete mupoptr;
			return gpo;
		}
		mupo.recalculateBoundary();
		return mupoptr;
	}
	
	///@param mapper maps from node-id to sserialize::spatial::GeoPoint
	template<typename T_ID_GEO_MAPPER>
	static sserialize::spatial::GeoMultiPolygon multiPolyFromClosedWays(const RawMultiPoly & p, const T_ID_GEO_MAPPER & mapper) {
		return multiPolyFromClosedWays(p.innerWays, p.outerWays, mapper);
	}
	
};

}}//end namespace detail::AreaExtractor

class AreaExtractor {
public:
	struct ValueType {
		std::unordered_map<std::string, std::string> kv;
		std::shared_ptr<sserialize::spatial::GeoRegion> * region;
	};
	typedef enum { ET_NONE=0, ET_BOUNDARIES=0x1, ET_LANDUSE=0x2, ET_AREA=0x4, ET_BUILDING=0x4|0x8, ET_ALL_BUT_BUILDINGS=ET_BOUNDARIES|ET_LANDUSE|ET_AREA, ET_ALL=ET_ALL_BUT_BUILDINGS|ET_BUILDING} ExtractionTypes;
public:
	AreaExtractor() {}
	virtual ~AreaExtractor() {}
	///@param processor (const std::shared_ptr<sserialize::spatial::GeoRegion> & region, osmpbf::IPrimitive & primitive), MUST be thread-safe
	template<typename TProcessor>
	bool extract(const std::string & inputFileName, TProcessor processor, ExtractionTypes extractionTypes = ET_ALL_BUT_BUILDINGS, bool needsName = true);
	static generics::RCPtr<osmpbf::AbstractTagFilter> createExtractionFilter(ExtractionTypes extractionTypes, bool needsName);
};


template<typename TProcessor>
bool AreaExtractor::extract(const std::string & inputFileName, TProcessor processor, ExtractionTypes extractionTypes, bool needsName) {
	if (! (extractionTypes & ET_ALL)) {
		return false;
	}

	osmpbf::OSMFileIn inFile(inputFileName, false);
	osmpbf::PrimitiveBlockInputAdaptor pbi;

	if (!inFile.open()) {
		std::cout << "Failed to open " <<  inputFileName << std::endl;
		return -1;
	}
	std::unordered_map<int64_t, detail::AreaExtractor::MultiPolyResolver::RawWay > rawWays;
	std::unordered_map<int64_t, sserialize::spatial::GeoPoint> nodes;
	
	std::unordered_set<int64_t> relevantWays;
	std::unordered_set<int64_t> relevantRelations;

	generics::RCPtr<osmpbf::AbstractTagFilter> mainFilter(createExtractionFilter(extractionTypes, needsName));
	
	mainFilter->assignInputAdaptor(&pbi);

	sserialize::ProgressInfo progress;
	
	progress.begin(inFile.dataSize(), "AreaExtractor: Nodes/Way-Refs");
	while (inFile.parseNextBlock(pbi)) {
		if (pbi.isNull())
			continue;


		progress(inFile.dataPosition());

		if (!mainFilter->buildIdCache()) {
			continue;
		}
		
		if (pbi.waysSize()) {
			for(osmpbf::IWayStream way(pbi.getWayStream()); !way.isNull(); way.next()) {
				if (way.refsSize() > 4 && *way.refBegin() == *(way.refBegin()+(way.refsSize()-1)) && mainFilter->matches(way)) {
					for(osmpbf::IWayStream::RefIterator it(way.refBegin()), end(way.refEnd()); it != end; ++it) {
						nodes[*it];
					}
					relevantWays.insert(way.id());
				}
			}
		}
		
		if (pbi.relationsSize()) {

			for (osmpbf::IRelationStream rel = pbi.getRelationStream(); !rel.isNull(); rel.next()) {
				if (mainFilter->matches(rel)) {
					for(osmpbf::IMemberStream mem = rel.getMemberStream(); !mem.isNull(); mem.next()) {
						if (mem.type() == osmpbf::WayPrimitive) {
							rawWays[mem.id()];
						}
					}
					relevantRelations.insert(rel.id());
				}
			}
		}
	}
	progress.end();
	
	uint32_t readBlobCount = omp_get_num_procs();
	bool processedFile = false;
	
	inFile.dataSeek(0);
	progress.begin(inFile.dataSize(), "AreaExtractor: Ways of relations");
	while (inFile.parseNextBlock(pbi)) {
		if (pbi.isNull())
			continue;
		progress(inFile.dataPosition());

		if (pbi.waysSize()) {
			for (osmpbf::IWayStream way = pbi.getWayStream(); !way.isNull(); way.next()) {
				int64_t wid = way.id();
				if (rawWays.count(wid)) {
					detail::AreaExtractor::MultiPolyResolver::RawWay & rw = rawWays[wid];
					rw.reserve(way.refsSize());
					for(osmpbf::IWayStream::RefIterator it(way.refBegin()), end(way.refEnd()); it != end; ++it) {
						int64_t wr = *it;
						rw.push_back(wr);
						nodes[wr];
					}
				}
			}
		}
	}
	progress.end();
	

	inFile.dataSeek(0);
	progress.begin(inFile.dataSize(), "AreaExtractor: Gathering nodes");
	while (!processedFile) {
		
		std::vector<osmpbf::BlobDataBuffer> pbiBuffers;
		inFile.getNextBlocks(pbiBuffers, readBlobCount);
		uint32_t pbiCount = pbiBuffers.size();
		processedFile = (pbiCount < readBlobCount);
		#pragma omp parallel for
		for(uint32_t i = 0; i < pbiCount; ++i) {
			osmpbf::PrimitiveBlockInputAdaptor pbi(pbiBuffers[i].data, pbiBuffers[i].availableBytes);
			pbiBuffers[i].clear();
			if (pbi.isNull()) {
				continue;
			}

			if (pbi.nodesSize()) {
				for(osmpbf::INodeStream node(pbi.getNodeStream()); !node.isNull(); node.next()) {
					int64_t nid = node.id();
					if (nodes.count( nid ) ) {
						sserialize::spatial::GeoPoint & gp = nodes.at(nid);
						gp.lat() = node.latd();
						gp.lon() = node.lond();
					}
				}
			}
		}
		progress(inFile.dataPosition());
	}
	progress.end();
	
	std::size_t relevantWaysSize = relevantWays.size();
	std::size_t relevantRelationsSize = relevantRelations.size();
	std::size_t assembledRelevantWays = 0;
	std::size_t assembledRelevantRelations = 0;
	
	inFile.dataSeek(0);
	progress.begin(inFile.dataSize(), "AreaExtractor: Assembling");
	while (inFile.parseNextBlock(pbi)) {
		if (pbi.isNull())
			continue;
		progress(inFile.dataPosition());

		if (pbi.waysSize()) {
			for(osmpbf::IWayStream way(pbi.getWayStream()); !way.isNull(); way.next()) {
				int64_t wid = way.id();
				if (!relevantWays.count(wid)) {
					continue;
				}
				relevantWays.erase(wid);
				std::shared_ptr<sserialize::spatial::GeoPolygon> gpo( new sserialize::spatial::GeoPolygon() );
				sserialize::spatial::GeoPolygon::PointsContainer & gpops = gpo->points();
				gpops.reserve(way.refsSize());
				for(osmpbf::IWayStream::RefIterator it(way.refBegin()), end(way.refEnd()); it != end; ++it) {
					int64_t wr = *it;
					if (nodes.count(wr)) {
						gpops.push_back(nodes[wr]);
					}
					else {
						std::cout << "Way has missing node" << std::endl;
						break;
					}
				}
				if (gpo->points().size()) {
					gpo->recalculateBoundary();
					processor(gpo, way);
					++assembledRelevantWays;
				}
			}
		}
		if (pbi.relationsSize()) {
			for(osmpbf::IRelationStream rel(pbi.getRelationStream()); !rel.isNull(); rel.next()) {
				int64_t relId = rel.id();
				
				if (!relevantRelations.count(relId)) {
					continue;
				}
				relevantRelations.erase(relId);
				std::vector< detail::AreaExtractor::MultiPolyResolver::RawWay > outerWays, innerWays;
				
				bool allWaysAvailable = true;
				
				for(osmpbf::IMemberStream mem = rel.getMemberStream(); !mem.isNull(); mem.next()) {
					if (mem.type() == osmpbf::WayPrimitive) {
						if (rawWays.count(mem.id()) && rawWays[mem.id()].size()) {
							const detail::AreaExtractor::MultiPolyResolver::RawWay & rw = rawWays[ mem.id() ];
							if (mem.role() == "outer" || mem.role() == "" || mem.role() == "exclave" || mem.role() == "Outer" || mem.role() == "outer:FIXME") {
								outerWays.push_back(rw);
							}
							else if (mem.role() == "inner" || mem.role() == "enclave") {
								innerWays.push_back(rw);
							}
							else {
								std::cout << "Illegal role in relation " << rel.id() << ": " << mem.role() << std::endl;
							}
						}
						else {
							allWaysAvailable = false;
						}
					}
				}
				if (outerWays.size() && !detail::AreaExtractor::MultiPolyResolver::multiPolyFromWays(innerWays, outerWays, innerWays, outerWays) && allWaysAvailable) {
					std::cout << "Failed to fully create MultiPolygon from multiple ways for relation " << rel.id() << std::endl;
				}
				if (outerWays.size()) {
					std::shared_ptr<sserialize::spatial::GeoRegion> gmpo( detail::AreaExtractor::MultiPolyResolver::multiPolyFromClosedWays(innerWays, outerWays, nodes) );
					processor(gmpo, rel);
					++assembledRelevantRelations;
				}
			}
		}
	}
	progress.end();
	
	std::cout << "Assembled " << assembledRelevantWays << "/" << relevantWaysSize << " boundary ways and " << assembledRelevantRelations << "/" << relevantRelationsSize << " boundary relations" << std::endl;
	
	inFile.close();
	
	return true;
}

}//end namespace

#endif