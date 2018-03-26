#ifndef OSMTOOLS_CELL_CRITERIA_H
#define OSMTOOLS_CELL_CRITERIA_H
#include <osmtools/OsmTriangulationRegionStore.h>

namespace osmtools {
namespace CellCriteria {
	
class CellTriangleCountCriteria: public ::osmtools::OsmTriangulationRegionStore::CellCriteriaInterface {
public:
	CellTriangleCountCriteria(uint32_t cellSizeTh);
	virtual ~CellTriangleCountCriteria() {}
	virtual int dataDependence() const override;
public:
	virtual bool init(const ::osmtools::OsmTriangulationRegionStore & store) override;
	virtual void begin() override;
	virtual void end() override;
	virtual bool refine(uint32_t cellId, const State & state) override;
	virtual CellCriteriaInterface * copy() override;
private:
	uint32_t m_cellSizeTh;
};

//This should be used in conjuction with a triangle refiner that bounds the maximum edge length
class CellDiagonalCriteria: public ::osmtools::OsmTriangulationRegionStore::CellCriteriaInterface {
public:
	CellDiagonalCriteria(double maxCellDiameter);
	virtual ~CellDiagonalCriteria() {}
	virtual int dataDependence() const override;
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
