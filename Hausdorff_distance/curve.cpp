#include "curve.h"
#include "utils.h"

#ifdef DEBUG
void PRINT_CTRLPTS(CubicBezierCurve* crv) {
	int i;
	printf("curve %p\n[\n", crv);
	for (i=0; i<4; ++i)
		printf("[%f, %f]\n", crv->control_pts[i][X], crv->control_pts[i][Y]);
	printf("]\n");
}
#endif

bool Arc::is_line() const{
	return this->radius != this->radius;
}

void evaluate(const CubicBezierCurve *curve, const REAL t, Point value)
{
	const REAL t_inv = 1.0f - t;
	const REAL t_inv_sq = t_inv * t_inv;
	const REAL t_sq = t * t;
	const REAL b0 = t_inv_sq * t_inv;
	const REAL b1 = 3 * t_inv_sq * t;
	const REAL b2 = 3 * t_inv * t_sq;
	const REAL b3 = t_sq * t;
	SET_VECTOR2(value, 0, 0);
	VECTOR2_X_SCALA_ADD(value, curve->control_pts[0], b0);
	VECTOR2_X_SCALA_ADD(value, curve->control_pts[1], b1);
	VECTOR2_X_SCALA_ADD(value, curve->control_pts[2], b2);
	VECTOR2_X_SCALA_ADD(value, curve->control_pts[3], b3);
}

void division_point(const Point p1, const Point p2, const REAL t, Point &p_out)
{
	p_out[0] = (p1[0] * t + p2[0] * (1.0 - t));
	p_out[1] = (p1[1] * t + p2[1] * (1.0 - t));
}

void middle_point(const Point p1, const Point p2, Point &p_out)
{
	p_out[0] = (p1[0] + p2[0]) / 2;
	p_out[1] = (p1[1] + p2[1]) / 2;
}

void sum_point(const Point p1, const Point p2, Point &p_out)
{
	p_out[0] = p1[0] + p2[0];
	p_out[1] = p1[1] + p2[1];
}

void subtract_point(const Point p1, const Point p2, Point &p_out)
{
	p_out[0] = p1[0] - p2[0];
	p_out[1] = p1[1] - p2[1];
}

void copy_point(const Point input, Point &output)
{
	output[0] = input[0];
	output[1] = input[1];
}

REAL distance(const Point &p1, const Point &p2){
	Point d;
	subtract_point(p1, p2, d);
	return norm(d);
}

REAL norm(const Point &p)
{
	return std::sqrt(p[0] * p[0] + p[1] * p[1]);
}

void normalize(Point &p)
{
	REAL norm_p = norm(p);
	p[0] /= norm_p;
	p[1] /= norm_p;
}

void get_bisection(const Point &p1, const Point &p2, Point &l_tan, Point &l_center)
{
	middle_point(p1, p2, l_center);
	subtract_point(p2, p1, l_tan);
	std::swap(l_tan[0], l_tan[1]);
	if (l_tan[0] < 0) l_tan[0] = -l_tan[0];
	else l_tan[1] = -l_tan[1];
	normalize(l_tan);
}

void get_line_intersection(const Point &l1_tan, const Point &l1_center, const Point &l2_tan, const Point &l2_center, Point &center)
{
	// circle cannot be formed if two tangent vector is parralel to one another
	if (std::abs(l1_tan[0] - l2_tan[0]) < EPS && std::abs(l1_tan[1] - l2_tan[1]) < EPS) {
		center[0] = NAN;
		center[1] = NAN;
	}
	else {
		// t * l1_tan + l1_center = p * l2_tan + l2_center
		REAL p;
		if (l1_tan[0] == 0){
			if (l2_center[0] == l1_center[0]){
				copy_point(l2_center, center);
			}
			else {
				p = (l1_center[0] - l2_center[0]) / l2_tan[0];
				center[0] = l1_center[0];
				center[1] = l2_center[1] + p * l2_tan[1];
			}
		}
		else if (l1_tan[1] == 0){
			if (l2_center[1] == l1_center[1]){
				copy_point(l2_center, center);
			}
			else {
				p = (l1_center[1] - l2_center[1]) / l2_tan[1];
				center[0] = l2_center[0] + p * l2_tan[0];
				center[1] = l1_center[1];
			}
		}
		else {
			if (std::abs(l2_tan[0] * l1_tan[1] - l2_tan[1] * l1_tan[0]) > EPS){
				// l1_tan[1] * (p * l2_tan[0] + l2_center[0] - l1_center[0]) = l1_tan[0] * (p * l2_tan[1] + l2_center[1] - l1_center[1])
				p = l1_tan[0] * (l2_center[1] - l1_center[1]) - l1_tan[1] * (l2_center[0] - l1_center[0]);
				p /= l2_tan[0] * l1_tan[1] - l2_tan[1] * l1_tan[0];
				center[0] = l2_center[0] + p * l2_tan[0];
				center[1] = l2_center[1] + p * l2_tan[1];
			}
			else {
				center[0] = NAN;
				center[1] = NAN;
			}
		}
	}
}

REAL atan(Point &p)
{
	return std::atan2(p[1], p[0]);
}
