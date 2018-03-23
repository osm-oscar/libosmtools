#ifndef OSMTOOLS_MESH_CRITERIA_H
#define OSMTOOLS_MESH_CRITERIA_H
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
//for ExtenedInt64q kernel
#include <libratss/CGAL/ExtendedInt64Cartesian.h>

#include <sserialize/spatial/DistanceCalculator.h>
#include <sserialize/utility/exceptions.h>

namespace osmtools {

struct Construct_refine_points {
public:
	typedef enum {T_CENTROID=1, T_LONGEST_EDGE=2, T_ON_EDGES=3} Type;
public:
	Construct_refine_points(Type t, const sserialize::spatial::DistanceCalculator & dc) : m_t(t), m_dc(dc) {}
	template<typename T_TDS, typename T_OUTPUT_ITERATOR>
	void calc(const typename T_TDS::Face_handle fh, T_OUTPUT_ITERATOR out) {
		using TDS = T_TDS;
		using Point = typename TDS::Point;
		switch (m_t) {
		case T_CENTROID:
			*out = CGAL::centroid(fh->vertex(0)->point(), fh->vertex(1)->point(), fh->vertex(2)->point());
			++out;
			break;
		case T_LONGEST_EDGE:
		{
			std::array<Point, 3> pts = {{
				fh->vertex(0)->point(),
				fh->vertex(1)->point(),
				fh->vertex(2)->point()
			}};
			std::array<double, 3> el;
			double longest = 0;
			for(int i(0); i < 3; ++i) {
				double p1x = CGAL::to_double(pts[i].x());
				double p1y = CGAL::to_double(pts[i].y());
				double p2x = CGAL::to_double(pts[TDS::cw(i)].x());
				double p2y = CGAL::to_double(pts[TDS::cw(i)].y());
				el[i] = m_dc.calc(p1x, p1y, p2x, p2y);
				longest = std::max(el[i], longest);
			}
			for(int i(0); i < 3; ++i) {
				if (el[i] == longest) {
					*out = CGAL::midpoint(fh->vertex(i)->point(), fh->vertex(TDS::cw(i))->point());
					++out;
				}
			}
			break;
		}
		case T_ON_EDGES:
			for(int i(0); i < 3; ++i) {
				*out = CGAL::midpoint(fh->vertex(i)->point(), fh->vertex(TDS::cw(i))->point());
				++out;
			}
			break;
		default:
			throw sserialize::InvalidEnumValueException("Construct_refine_points with type="  + std::to_string(m_t));
			break;
		};
	}
private:
	int m_t;
	sserialize::spatial::DistanceCalculator m_dc;
};

template<typename TDS>
class CentroidDistanceBaseMeshCriteria {
public:
	typedef typename TDS::Face_handle Face_handle;
	typedef typename TDS::Point Point;
	struct Is_bad {
		sserialize::spatial::DistanceCalculator m_dc;
		Is_bad(const sserialize::spatial::DistanceCalculator & dc) : m_dc(dc) {}
		double maxCentroidDist(Face_handle fh) const {
			Point p = CGAL::centroid(fh->vertex(0)->point(), fh->vertex(1)->point(), fh->vertex(2)->point());
			double latp = CGAL::to_double(p.x());
			double lonp = CGAL::to_double(p.y());
			double q = 0.0;
			for(int j(0); j < 3; ++j) {
				Point fp( fh->vertex(j)->point() );
				double lat = CGAL::to_double(fp.x());
				double lon = CGAL::to_double(fp.y());
				double tmp = m_dc.calc(latp, lonp, lat, lon);
				q = std::max<double>(tmp, q);
			}
			//clip centroid distance to 1.0 (which is pretty small)
			q = std::max<double>(q, 1.0); 
			return q;
		}
	};
	
	using Construct_refine_points = ::osmtools::Construct_refine_points;
	
public:
	CentroidDistanceBaseMeshCriteria(Construct_refine_points::Type crpt = Construct_refine_points::T_CENTROID) : m_dc(sserialize::spatial::DistanceCalculator::DCT_GEODESIC_ACCURATE),
	m_crpt(crpt)
	{}
public:
	void crpt(Construct_refine_points::Type type) {
		m_crpt = type;
	}
	Construct_refine_points::Type crpt() const { return m_crpt; }
public:
	Construct_refine_points construct_refine_points_object() {
		return Construct_refine_points(m_crpt, dc());
	}
public:
	static bool usesCellIds() { return false; }
protected:
	const sserialize::spatial::DistanceCalculator & dc() const { return m_dc; }
protected:
	sserialize::spatial::DistanceCalculator m_dc;
	Construct_refine_points::Type m_crpt;
};

template<typename TDS>
class CentroidDistanceMeshCriteria: public CentroidDistanceBaseMeshCriteria<TDS> {
public:
	typedef CentroidDistanceBaseMeshCriteria<TDS> MyParentClass;
	typedef typename MyParentClass::Face_handle Face_handle;
	typedef typename MyParentClass::Point Point;
public:
	typedef double Quality;
	struct Is_bad: MyParentClass::Is_bad {
		double m_r;
		sserialize::spatial::DistanceCalculator m_dc;
		Is_bad(double maxDist, const sserialize::spatial::DistanceCalculator & dc) : MyParentClass::Is_bad(dc), m_r(maxDist), m_dc(dc) {}
		
		inline Quality quality(CGAL::Mesh_2::Face_badness fb) const {
			if (fb == CGAL::Mesh_2::NOT_BAD) {
				return m_r;
			}
			else {
				return std::numeric_limits<double>::max();
			}
		}
		
		CGAL::Mesh_2::Face_badness operator()(Quality q) const {
			if (q > m_r) {
				return CGAL::Mesh_2::IMPERATIVELY_BAD;
			}
			else {
				return CGAL::Mesh_2::NOT_BAD;
			}
		}
		
		CGAL::Mesh_2::Face_badness operator()(Face_handle fh, Quality & q) const {
			q = MyParentClass::Is_bad::maxCentroidDist(fh);
			return (*this)(q);
		};
	};
private:
	double m_r;
public:
	CentroidDistanceMeshCriteria(double maxDist) : m_r(maxDist) {}
	Is_bad is_bad_object() const { return Is_bad(m_r, MyParentClass::dc()); }
};

template<typename TDS>
class EdgeLengthMeshCriteria: public CentroidDistanceBaseMeshCriteria<TDS> {
public:
	typedef CentroidDistanceBaseMeshCriteria<TDS> MyParentClass;
	typedef typename MyParentClass::Face_handle Face_handle;
	typedef typename MyParentClass::Point Point;
public:
	typedef double Quality;
	struct Is_bad {
		double m_r;
		sserialize::spatial::DistanceCalculator m_dc;
		Is_bad(double maxLen, const sserialize::spatial::DistanceCalculator & dc) : m_r(maxLen), m_dc(dc) {}
		
		CGAL::Mesh_2::Face_badness operator()(Quality q) const {
			if (q > m_r) {
				return CGAL::Mesh_2::IMPERATIVELY_BAD;
			}
			else {
				return CGAL::Mesh_2::NOT_BAD;
			}
		}
		
		Quality quality(CGAL::Mesh_2::Face_badness fb) const {
			if (fb == CGAL::Mesh_2::NOT_BAD) {
				return m_r;
			}
			else {
				return std::numeric_limits<double>::max();
			}
		}
		
		CGAL::Mesh_2::Face_badness operator()(Face_handle fh, Quality & q) const {
			std::array<Point, 3> pts = {{
				fh->vertex(0)->point(),
				fh->vertex(1)->point(),
				fh->vertex(2)->point()
			}};
			q = 0;
			for(int i(0); i < 3; ++i) {
				double p1x = CGAL::to_double(pts[i].x());
				double p1y = CGAL::to_double(pts[i].y());
				double p2x = CGAL::to_double(pts[TDS::cw(i)].x());
				double p2y = CGAL::to_double(pts[TDS::cw(i)].y());
				double dist = m_dc.calc(p1x, p1y, p2x, p2y);
				q = std::max(dist, q);
			}
			return (*this)(q);
		};
	};
private:
	double m_r;
public:
	///@param maxLength the maxium length of the longest edge
	EdgeLengthMeshCriteria(double maxLength) :
		MyParentClass(Construct_refine_points::T_LONGEST_EDGE),
		m_r(maxLength)
	{}
	Is_bad is_bad_object() const { return Is_bad(m_r, MyParentClass::dc()); }
};

template<typename TDS>
class EdgeLengthRatioMeshCriteria: public CentroidDistanceBaseMeshCriteria<TDS> {
public:
	typedef CentroidDistanceBaseMeshCriteria<TDS> MyParentClass;
	typedef typename MyParentClass::Face_handle Face_handle;
	typedef typename MyParentClass::Point Point;
public:
	typedef double Quality;
	struct Is_bad {
		double m_r;
		sserialize::spatial::DistanceCalculator m_dc;
		Is_bad(double maxRatio, const sserialize::spatial::DistanceCalculator & dc) : m_r(maxRatio), m_dc(dc) {}
		
		CGAL::Mesh_2::Face_badness operator()(Quality q) const {
			if (q > m_r) {
				return CGAL::Mesh_2::IMPERATIVELY_BAD;
			}
			else {
				return CGAL::Mesh_2::NOT_BAD;
			}
		}
		
		Quality quality(CGAL::Mesh_2::Face_badness fb) const {
			if (fb == CGAL::Mesh_2::NOT_BAD) {
				return m_r;
			}
			else {
				return std::numeric_limits<double>::max();
			}
		}
		
		CGAL::Mesh_2::Face_badness operator()(Face_handle fh, Quality & q) const {
			std::array<Point, 3> pts = {{
				fh->vertex(0)->point(),
				fh->vertex(1)->point(),
				fh->vertex(2)->point()
			}};
			double longest = 0;
			double shortest = std::numeric_limits<double>::max();
			for(int i(0); i < 3; ++i) {
				double p1x = CGAL::to_double(pts[i].x());
				double p1y = CGAL::to_double(pts[i].y());
				double p2x = CGAL::to_double(pts[TDS::cw(i)].x());
				double p2y = CGAL::to_double(pts[TDS::cw(i)].y());
				double dist = m_dc.calc(p1x, p1y, p2x, p2y);
				longest = std::max(dist, longest);
				shortest = std::min(dist, shortest);
			}
			shortest = std::max<double>(shortest, std::numeric_limits<double>::epsilon());
			longest = std::max<double>(longest, std::numeric_limits<double>::epsilon());
			if (longest < 1.0 && shortest < 1.0) { //very small triangle, refining these would be useless (even 1 meter may be to big)
				q = 1.0;
			}
			else {
				q = longest/shortest;
			}
			return (*this)(q);
		};
	};
private:
	double m_r;
public:
	///@param maxRatio the maxium ratio between the shortest and longest edge
	EdgeLengthRatioMeshCriteria(double maxRatio) :
		MyParentClass(Construct_refine_points::T_LONGEST_EDGE),
		m_r(maxRatio)
	{}
	Is_bad is_bad_object() const { return Is_bad(m_r, MyParentClass::dc()); }
};


template<typename TDS>
class LipschitzMeshCriteria: public CentroidDistanceBaseMeshCriteria<TDS> {
public:
	typedef CentroidDistanceBaseMeshCriteria<TDS> MyParentClass;
	typedef typename MyParentClass::Face_handle Face_handle;
	typedef typename MyParentClass::Point Point;
	typedef double Quality;
	struct Is_bad: MyParentClass::Is_bad {
		double m_s;
		TDS * m_tds;
		Is_bad(double maxDist, const sserialize::spatial::DistanceCalculator & dc, TDS * tds) : MyParentClass::Is_bad(dc), m_s(maxDist), m_tds(tds) {}
		
		inline Quality quality(CGAL::Mesh_2::Face_badness fb) const {
			if (fb == CGAL::Mesh_2::NOT_BAD) {
				return m_s;
			}
			else {
				return std::numeric_limits<double>::max();
			}
		}
		
		CGAL::Mesh_2::Face_badness operator()(Quality q) const {
			if (q > m_s) {
				return CGAL::Mesh_2::IMPERATIVELY_BAD;
			}
			else {
				return CGAL::Mesh_2::NOT_BAD;
			}
		}
		
		CGAL::Mesh_2::Face_badness operator()(Face_handle fh, Quality & q) {
			double maxSlope = 0.0;
			double myCD = MyParentClass::Is_bad::maxCentroidDist(fh);
			//by definition: slope goes from nfh to fh
			//=> if this is the smallest triangle, then the slope is negative for all neighbors, if not it's positive
			for(int j(0); j < 3; ++j) {
				Face_handle nfh(fh->neighbor(j));
				if (!m_tds->is_infinite(nfh)) {
					double nfhD = MyParentClass::Is_bad::maxCentroidDist(nfh);
					if (myCD > nfhD) {
						maxSlope = std::max<double>(myCD / nfhD, maxSlope);
					}
				}
			}
			q = maxSlope;
			return (*this)(q);
		};
	};
private:
	double m_s;
	TDS * m_tds;
public:
	LipschitzMeshCriteria(double maxSlope, TDS * tds) : m_s(maxSlope), m_tds(tds) {}
	Is_bad is_bad_object() const { return Is_bad(m_s, MyParentClass::dc(), m_tds); }
};

template<typename T_BASE_MESH_CRITERIA>
class RefineTrianglesWithCellIdMeshCriteria: public T_BASE_MESH_CRITERIA {
public:
	typedef T_BASE_MESH_CRITERIA MyParentClass;
	typedef typename MyParentClass::Face_handle Face_handle;
	typedef typename MyParentClass::Point Point;
	typedef typename MyParentClass::Quality Quality;
	struct Is_bad: MyParentClass::Is_bad {
		Is_bad(const typename MyParentClass::Is_bad & base) : MyParentClass::Is_bad(base) {}
		
		CGAL::Mesh_2::Face_badness operator()(Quality q) const {
			return MyParentClass::Is_bad::operator()(q);
		}
		
		CGAL::Mesh_2::Face_badness operator()(Face_handle fh, Quality & q) {
			if (fh->info().hasCellId()) {
				return MyParentClass::Is_bad::operator()(fh, q);
			}
			else {
				q = MyParentClass::Is_bad::quality(CGAL::Mesh_2::NOT_BAD);
				return CGAL::Mesh_2::NOT_BAD;
			}
		};
	};
public:
	RefineTrianglesWithCellIdMeshCriteria(const MyParentClass & base) : MyParentClass(base) {}
	Is_bad is_bad_object() const { return Is_bad(MyParentClass::is_bad_object()); }
	static bool usesCellIds() { return true; }
};

}//end namespace osmtools

#endif
