#ifndef OSMTOOLS_CELL_CRITERIA_H
#define OSMTOOLS_CELL_CRITERIA_H
#include <osmtools/OsmTriangulationRegionStore.h>

namespace osmtools {
namespace CellCriteria {
	
class CellTriangleCountCriteria: public ::osmtools::OsmTriangulationRegionStore::CellCriteriaInterface {
public:
	CellTriangleCountCriteria(uint32_t cellSizeTh);
	virtual ~CellTriangleCountCriteria() {}
public:
	virtual bool init(const ::osmtools::OsmTriangulationRegionStore & store) override;
	virtual void begin() override;
	virtual void end() override;
	virtual bool refine(uint32_t cellId, const State & state) override;
	virtual CellCriteriaInterface * copy() override;
private:
	uint32_t m_cellSizeTh;
};

class CellDiagonalCriteria: public ::osmtools::OsmTriangulationRegionStore::CellCriteriaInterface {
public:
	CellDiagonalCriteria(double maxCellDiameter);
	virtual ~CellDiagonalCriteria() {}
public:
	virtual bool init(const ::osmtools::OsmTriangulationRegionStore & store) override;
	virtual void begin() override;
	virtual void end() override;
	virtual bool refine(const State & state) override;
	virtual bool refine(uint32_t cellId, const State & state) override;
	virtual CellCriteriaInterface * copy() override;
private:
	double m_maxCellDiameter;
	sserialize::spatial::DistanceCalculator m_dc;
};
	
}}//end namespace osmtools::CellCriteria

#endif
