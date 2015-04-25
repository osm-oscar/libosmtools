#ifndef OSM_TOOLS_AREA_EXTRACTOR_H
#define OSM_TOOLS_AREA_EXTRACTOR_H
#include <string>
#include <cstdlib>
#include <omp.h>
#include <sserialize/spatial/GeoMultiPolygon.h>
#include <sserialize/spatial/GeoRegionStore.h>
#include <osmpbf/osmfile.h>
#include <osmpbf/primitiveblockinputadaptor.h>
#include <osmpbf/irelation.h>
#include <osmpbf/iway.h>
#include <osmpbf/inode.h>
#include <osmtools/AreaExtractorFilters.h>
#include <osmpbf/parsehelpers.h>

/** This is a simple class to extract areas from OpenStreetMap pbf-files
  * It needs the osmpbf and sserialize libraries in global include paths
  */

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

class AreaExtractor: public detail::AreaExtractor::Base {
public:
	struct ValueType {
		std::unordered_map<std::string, std::string> kv;
		std::shared_ptr<sserialize::spatial::GeoRegion> * region;
	};
private:
	typedef std::pair<double, double> Point;
	struct Context {
		std::unordered_map<int64_t, detail::AreaExtractor::MultiPolyResolver::RawWay > rawWays;
		std::mutex rawWaysLock;
		
		std::unordered_map<int64_t, Point> nodes;
		std::mutex nodesLock;
		
		osmpbf::OSMFileIn inFile;
		ExtractionTypes extractionTypes;
		bool needsName;
		sserialize::ProgressInfo pinfo;
		
		std::atomic<uint32_t> relevantWaysSize;
		std::atomic<uint32_t> relevantRelationsSize;
		
		std::atomic<uint32_t> assembledRelevantWays;
		std::atomic<uint32_t> assembledRelevantRelations;
		
		Context(const std::string & filename) : inFile(filename), relevantWaysSize(0), relevantRelationsSize(0), assembledRelevantWays(0), assembledRelevantRelations(0) {}
	};
	struct ExtractorCallBack {
		virtual void operator()(const std::shared_ptr<sserialize::spatial::GeoRegion> & region, osmpbf::IPrimitive & primitive) = 0;
	};
	struct ExtractionFunctorBase {
		Context * ctx;
		osmpbf::PrimitiveBlockInputAdaptor * pbi;
		generics::RCPtr<osmpbf::AbstractTagFilter> mainFilter;
		void assignInputAdaptor(osmpbf::PrimitiveBlockInputAdaptor & pbi);
		ExtractionFunctorBase(Context * ctx);
		ExtractionFunctorBase(const ExtractionFunctorBase & other);
		~ExtractionFunctorBase() {}
	};
	//Private compressor
	struct RelationWaysExtractor: ExtractionFunctorBase {
		RelationWaysExtractor(Context * ctx) : ExtractionFunctorBase(ctx) {}
		RelationWaysExtractor(const RelationWaysExtractor & o) : ExtractionFunctorBase(o) {}
		void operator()(osmpbf::PrimitiveBlockInputAdaptor & pbi);
	};
	//Private compressor
	struct WayRefsExtractor: ExtractionFunctorBase {
		WayRefsExtractor(Context * ctx) : ExtractionFunctorBase(ctx) {}
		WayRefsExtractor(const WayRefsExtractor & o) : ExtractionFunctorBase(o) {}
		void operator()(osmpbf::PrimitiveBlockInputAdaptor & pbi);
	};
	//No private processor
	struct RelationWayNodeRefsExtractor {
		Context * ctx;
		RelationWayNodeRefsExtractor(Context * ctx) : ctx(ctx) {}
		RelationWayNodeRefsExtractor(const RelationWayNodeRefsExtractor & o) : ctx(o.ctx) {}
		void operator()(osmpbf::PrimitiveBlockInputAdaptor & pbi);
	};
	//Private compressor
	struct WayExtractor: ExtractionFunctorBase {
		ExtractorCallBack * cb;
		WayExtractor(Context * ctx, ExtractorCallBack * cb) : ExtractionFunctorBase(ctx), cb(cb) {}
		WayExtractor(const WayExtractor & o) : ExtractionFunctorBase(o), cb(o.cb) {}
		void operator()(osmpbf::PrimitiveBlockInputAdaptor & pbi);
	};
	//Private compressor
	struct RelationExtractor: ExtractionFunctorBase {
		ExtractorCallBack * cb;
		RelationExtractor(Context * ctx, ExtractorCallBack * cb) : ExtractionFunctorBase(ctx), cb(cb) {}
		RelationExtractor(const RelationExtractor & o) : ExtractionFunctorBase(o), cb(o.cb) {}
		void operator()(osmpbf::PrimitiveBlockInputAdaptor & pbi);
	};
	//No private processor
	struct NodeGatherer {
		Context * ctx;
		NodeGatherer(Context * ctx) : ctx(ctx) {}
		NodeGatherer(const NodeGatherer & o) : ctx(o.ctx) {}
		void operator()(osmpbf::PrimitiveBlockInputAdaptor & pbi);
	};
public:
	AreaExtractor() {}
	virtual ~AreaExtractor() {}
	///@param processor (const std::shared_ptr<sserialize::spatial::GeoRegion> & region, osmpbf::IPrimitive & primitive), MUST be thread-safe
	template<typename TProcessor>
	bool extract(const std::string & inputFileName, TProcessor processor, ExtractionTypes extractionTypes = ET_ALL_BUT_BUILDINGS, bool needsName = true);
};

template<typename TProcessor>
bool AreaExtractor::extract(const std::string & inputFileName, TProcessor processor, ExtractionTypes extractionTypes, bool needsName) {
	if (! (extractionTypes & ET_ALL)) {
		return false;
	}
	
	Context ctx(inputFileName);
	ctx.extractionTypes = extractionTypes;
	ctx.needsName = needsName;

	osmpbf::PrimitiveBlockInputAdaptor pbi;

	if (!ctx.inFile.open()) {
		std::cout << "Failed to open " <<  inputFileName << std::endl;
		return -1;
	}
	
	struct MyCB: ExtractorCallBack {
		TProcessor * processor;
		virtual void operator()(const std::shared_ptr< sserialize::spatial::GeoRegion >& region, osmpbf::IPrimitive& primitive) {
			(*processor)(region, primitive);
		}
	} cb;
	cb.processor = &processor;
	
	NodeGatherer ng(&ctx);
	
	{
		WayRefsExtractor wre(&ctx);
		WayExtractor we(&ctx, &cb);

		ctx.pinfo.begin(ctx.inFile.dataSize(), "AreaExtractor: Ways' node-refs");
		ctx.inFile.reset();
		osmpbf::parseFileCPPThreads(ctx.inFile, wre, 0, 1, true);
		ctx.pinfo.end();

		ctx.pinfo.begin(ctx.inFile.dataSize(), "AreaExtractor: Ways' nodes");
		ctx.inFile.reset();
		osmpbf::parseFileCPPThreads(ctx.inFile, ng, 0, 1, false);
		ctx.pinfo.end();
		
		ctx.pinfo.begin(ctx.inFile.dataSize(), "AreaExtractor: Assembling ways");
		ctx.inFile.reset();
		osmpbf::parseFileCPPThreads(ctx.inFile, WayExtractor(&ctx, &cb), 0, 1, true);
		ctx.pinfo.end();
	}
	
	ctx.nodes.clear();
	{
		RelationWaysExtractor rwe(&ctx);
		RelationWayNodeRefsExtractor rwnr(&ctx);
		RelationExtractor re(&ctx, &cb);

		ctx.pinfo.begin(ctx.inFile.dataSize(), "AreaExtractor: Relation's ways");
		ctx.inFile.reset();
		osmpbf::parseFileCPPThreads(ctx.inFile, rwe, 0, 1, true);
		ctx.pinfo.end();

		ctx.pinfo.begin(ctx.inFile.dataSize(), "AreaExtractor: Relation-ways' node-refs");
		ctx.inFile.reset();
		osmpbf::parseFileCPPThreads(ctx.inFile, rwnr, 0, 1, false);
		ctx.pinfo.end();
		
		ctx.pinfo.begin(ctx.inFile.dataSize(), "AreaExtractor: Relation-ways' nodes");
		ctx.inFile.reset();
		osmpbf::parseFileCPPThreads(ctx.inFile, ng, 0, 1, false);
		ctx.pinfo.end();
		
		ctx.pinfo.begin(ctx.inFile.dataSize(), "AreaExtractor: Assembling relations");
		ctx.inFile.reset();
		osmpbf::parseFileCPPThreads(ctx.inFile, re, 0, 1, true);
		ctx.pinfo.end();
	}
	ctx.nodes.clear();

	std::cout << "Assembled " << ctx.assembledRelevantWays << "/" << ctx.relevantWaysSize << " boundary ways and " << ctx.assembledRelevantRelations << "/" << ctx.relevantRelationsSize << " boundary relations" << std::endl;
	
	ctx.inFile.close();
	
	return true;
}

}//end namespace

#endif