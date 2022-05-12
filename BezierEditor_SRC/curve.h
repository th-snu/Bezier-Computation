#ifndef _CURVE_H_
#define _CURVE_H_

#define PRECISION   1e-5
#define EPS         1e-6        /* data type is float */
#define INF   FLT_MAX

#include <vector>

typedef float REAL;
typedef REAL  Point[2];

typedef struct CubicBezierCurve 
{
	Point control_pts[4];
} CubicBezierCurve;

typedef struct Arc
{
	Point center;
	REAL radius;
	REAL begin;
	REAL end;
} Arc;

#ifdef DEBUG
void PRINT_CTRLPTS(CubicBezierCurve* crv);
#else
#   define PRINT_CTRLPTS(X)
#endif

#define SET_PT2(V, V1, V2) do { (V)[0] = (V1); (V)[1] = (V2); } while (0)

void evaluate(const CubicBezierCurve *curve, const REAL t, Point value);

void middle_point(const Point p1, const Point p2, const Point p_out);

void sum_point(const Point p1, const Point p2, Point &p_out);

void subtract_point(const Point p1, const Point p2, Point &p_out);

void copy_point(const Point input, Point &output);

void subdivide(const CubicBezierCurve *curve, const CubicBezierCurve *output);

void subdivide(const CubicBezierCurve *curve, std::vector<CubicBezierCurve> &output, int power);

REAL norm(Point p);

void normalize(Point &p);

void get_tangent(const CubicBezierCurve *curve, Point &tan_begin, Point &tan_end);

void get_bisection(const Point &p1, const Point &p2, Point &l_tan, Point &l_center);

void get_line_segs(const CubicBezierCurve *curve, Point &l1_tan, Point &l1_center, Point &l2_tan, Point &l2_center);

void get_line_intersection(const Point &l1_tan, const Point &l1_center, const Point &l2_tan, const Point &l2_center, Point &center);

void set_arc_center(const CubicBezierCurve *curve, const Point &arc1_point, const Point &arc2_point, const Point &inflect, Arc *arc1, Arc *arc2);

void get_biarc_inflect(const CubicBezierCurve *curve, Point &inflect, bool mode);

void to_biarc(const CubicBezierCurve *curve, Point &inflect, Arc *arc1, Arc *arc2);

#endif /* _CURVE_H_ */
