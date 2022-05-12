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

class Arc
{
public:
	Point center;
	REAL radius;
	REAL begin;
	REAL end;

	bool is_line() const;
};

#ifdef DEBUG
void PRINT_CTRLPTS(CubicBezierCurve* crv);
#else
#   define PRINT_CTRLPTS(X)
#endif

#define SET_PT2(V, V1, V2) do { (V)[0] = (V1); (V)[1] = (V2); } while (0)
#define IS_NAN(V) (V != V)

void evaluate(const CubicBezierCurve *curve, const REAL t, Point value);

void middle_point(const Point p1, const Point p2, Point &p_out);

void division_point(const Point p1, const Point p2, const REAL t, Point &p_out);

void sum_point(const Point p1, const Point p2, Point &p_out);

void subtract_point(const Point p1, const Point p2, Point &p_out);

void copy_point(const Point input, Point &output);

REAL distance(const Point &p1, const Point &p2);

REAL norm(const Point &p);

void normalize(Point &p);

void get_bisection(const Point &p1, const Point &p2, Point &l_tan, Point &l_center);

void get_line_intersection(const Point &l1_tan, const Point &l1_center, const Point &l2_tan, const Point &l2_center, Point &center);

REAL atan(Point &p);

#endif /* _CURVE_H_ */
