#include "curve.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

enum XY { X = 0, Y };

#define SET_VECTOR2(V, V1, V2)          do { (V)[X] = (V1); (V)[Y] = (V2); } while (0)
#define COPY_PT(DST, SRC)               do { (DST)[X] = SRC[X]; (DST)[Y] = SRC[Y]; } while (0)
#define VECTOR2_X_SCALA_ADD(O, V, S)    do { O[X] += (S) * (V)[X]; O[Y] += (S) * (V)[Y]; } while (0)

#ifdef DEBUG
void PRINT_CTRLPTS(CubicBezierCurve* crv) {
	int i;
	printf("curve %p\n[\n", crv);
	for (i=0; i<4; ++i)
		printf("[%f, %f]\n", crv->control_pts[i][X], crv->control_pts[i][Y]);
	printf("]\n");
}
#endif

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

void subdivide(const CubicBezierCurve *curve, CubicBezierCurve *output1, CubicBezierCurve *output2)
{
	Point inter1[3];
	Point inter2[2];
	Point inter3[1];

	for(int i = 0; i < 3; i++)
		middle_point(curve->control_pts[i], curve->control_pts[i+1], inter1[i]);
	for(int i = 0; i < 2; i++)
		middle_point(inter1[i], inter1[i+1], inter2[i]);
	for(int i = 0; i < 1; i++)
		middle_point(inter2[i], inter2[i+1], inter3[i]);

	copy_point(curve->control_pts[0], output1->control_pts[0]);
	copy_point(inter1[0], output1->control_pts[1]);
	copy_point(inter2[0], output1->control_pts[2]);
	copy_point(inter3[0], output1->control_pts[3]);

	copy_point(curve->control_pts[3], output2->control_pts[0]);
	copy_point(inter1[2], output2->control_pts[1]);
	copy_point(inter2[1], output2->control_pts[2]);
	copy_point(inter3[0], output2->control_pts[3]);
}

void subdivide(const CubicBezierCurve *curve, std::vector<CubicBezierCurve> &output, int power)
{
	std::vector<CubicBezierCurve> prev, curr;
	curr.push_back(*curve);

	for(int i = 0; i < power; i++){
		prev = curr;
		curr.clear();
		for(auto subcurve: prev){
			CubicBezierCurve curve1, curve2;
			subdivide(&subcurve, &curve1, &curve2);
			curr.push_back(curve1);
			curr.push_back(curve2);
		}
	}

	std::copy(curr.begin(), curr.end(), back_inserter(output));
}

REAL norm(Point p)
{
	return std::sqrt(p[0] * p[0] + p[1] * p[1]);
}

void normalize(Point &p)
{
	REAL norm_p = norm(p);
	p[0] /= norm_p;
	p[1] /= norm_p;
}

void get_tangent(const CubicBezierCurve *curve, Point &tan_begin, Point &tan_end)
{
	for (int i = 1; i <= 4; i++){
		if (i == 4){
			tan_begin[0] = 0.0;
			tan_begin[1] = 1.0;
			break;
		}
		subtract_point(curve->control_pts[i], curve->control_pts[i-1], tan_begin);
		if (norm(tan_begin) > EPS) break;
	}
	normalize(tan_begin);

	for (int i = 2; i >= -1; i++){
		if (i == -1){
			tan_end[0] = 1.0;
			tan_end[1] = 0.0;
			break;
		}
		subtract_point(curve->control_pts[i+1], curve->control_pts[i], tan_end);
		if (norm(tan_end) > EPS) break;
	}
	normalize(tan_end);
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

void get_line_segs(const CubicBezierCurve *curve, Point &l1_tan, Point &l1_center, Point &l2_tan, Point &l2_center)
{
	Point tan_begin, tan_end;
	get_tangent(curve, tan_begin, tan_end);

	// Find vertical bisector of two control points summated with tangent
	Point p1, p2;
	sum_point(curve->control_pts[0], tan_begin, p1);
	sum_point(curve->control_pts[3], tan_end, p2);

	get_bisection(p1, p2, l1_tan, l1_center);
	get_bisection(curve->control_pts[0], curve->control_pts[3], l2_tan, l2_center);
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

void get_biarc_inflect(const CubicBezierCurve *curve, Point &inflect, bool mode)
{
	Point l1_center, l1_tan, l2_center, l2_tan, center;
	get_line_segs(curve, l1_tan, l1_center, l2_tan, l2_center);
	get_line_intersection(l1_tan, l1_center, l2_tan, l2_center, center);

	Point r_vec;
	subtract_point(center, curve->control_pts[0], r_vec);
	REAL r = norm(r_vec);

	if (!mode){
		// two lines does not intersect, a line is the best approximation available instead of a circle
		if (center[0] != center[0]) {
			copy_point(l2_center, inflect);
			return;
		}

		Point r_tan;
		r_tan[0] = r * l2_tan[0];
		r_tan[1] = r * l2_tan[1];

		Point i_plus, i_minus;
		Point d_plus_vec, d_minus_vec;
		REAL d_plus, d_minus;

		sum_point(center, r_tan, i_plus);
		subtract_point(center, r_tan, i_minus);
		subtract_point(curve->control_pts[0], i_plus, d_plus_vec);
		subtract_point(curve->control_pts[0], i_minus, d_minus_vec);

		d_plus = norm(d_plus_vec);
		d_minus = norm(d_minus_vec);

		if (d_plus < d_minus) copy_point(i_plus, inflect);
		else copy_point(i_minus, inflect);
	}

	// use actual intersection with simple algorithm
	else {
		// does not work yet
		REAL min_d = INF;
		for (int i = 0; i < 5000; i++){
			Point p, d_vec;
			evaluate(curve, (double)i / 5000, p);
			subtract_point(center, p, d_vec);
			REAL d = norm(d_vec);
			REAL e = std::abs(r - d);
			if (e < min_d) copy_point(p, inflect);
		}
	}
}

REAL atan(Point &p)
{
	return std::atan2(p[1], p[0]);
}

void set_arc_center(const CubicBezierCurve *curve, const Point &arc1_point, const Point &arc2_point, const Point &inflect, Arc *arc1, Arc *arc2)
{
	Point arc1_tan, arc2_tan;
	get_tangent(curve, arc1_tan, arc2_tan);

	// find middle point of each arc's endpoint
	Point bi1_tan, bi1_center;
	Point bi2_tan, bi2_center;

	get_bisection(arc1_point, inflect, bi1_tan, bi1_center);
	get_bisection(arc2_point, inflect, bi2_tan, bi2_center);

	// find normal vector for endpoint of each arc
	Point arc1_norm, arc2_norm;
	
	arc1_norm[0] = arc1_tan[1];
	arc1_norm[1] = arc1_tan[0];
	if (arc1_norm[0] < 0) arc1_norm[0] = -arc1_norm[0];
	else arc1_norm[1] = -arc1_norm[1];

	arc2_norm[0] = arc2_tan[1];
	arc2_norm[1] = arc2_tan[0];
	if (arc2_norm[0] < 0) arc2_norm[0] = -arc2_norm[0];
	else arc2_norm[1] = -arc2_norm[1];

	// find center of two arcs
	Point center1, center2;
	get_line_intersection(bi1_tan, bi1_center, arc1_norm, arc1_point, center1);
	get_line_intersection(bi2_tan, bi2_center, arc2_norm, arc2_point, center2);

	if(center1[0] != center1[0]) {
		center1[0] = NAN;
		center1[1] = NAN;
	}
	if(center2[0] != center2[0]) {
		center2[0] = NAN;
		center2[1] = NAN;
	}

	copy_point(center1, arc1->center);
	copy_point(center2, arc2->center);
}

void to_biarc(const CubicBezierCurve *curve, Point &inflect, Arc *arc1, Arc *arc2)
{
	const Point &arc1_point = curve->control_pts[0];
	const Point &arc2_point = curve->control_pts[3];
	set_arc_center(curve, arc1_point, arc2_point, inflect, arc1, arc2);

	const Point &center1 = arc1->center;
	const Point &center2 = arc2->center;

	Point r_vec1, r_vec2;
	subtract_point(arc1_point, center1, r_vec1);
	subtract_point(arc2_point, center2, r_vec2);
	arc1->radius = norm(r_vec1);
	arc2->radius = norm(r_vec2);

	Point arc1_end_on_circle, arc1_begin_on_circle;
	Point arc2_end_on_circle, arc2_begin_on_circle;

	subtract_point(arc1_point, center1, arc1_begin_on_circle);
	subtract_point(inflect, center1, arc1_end_on_circle);
	subtract_point(arc2_point, center2, arc2_end_on_circle);
	subtract_point(inflect, center2, arc2_begin_on_circle);

	arc1->begin = atan(arc1_begin_on_circle);
	arc1->end = atan(arc1_end_on_circle);

	arc2->begin = atan(arc2_begin_on_circle);
	arc2->end = atan(arc2_end_on_circle);

	// If arc is going counterclock-wise, swap begin and end
	normalize(r_vec1);
	normalize(r_vec2);
	Point arc1_tan, arc2_tan;
	get_tangent(curve, arc1_tan, arc2_tan);

	if (r_vec1[1] * arc1_tan[0] - r_vec1[0] * arc1_tan[1] > 0.0)
		std::swap(arc1->begin, arc1->end);
	if (arc1->end < arc1->begin) arc1->end += 2 * M_PI;

	if (r_vec2[1] * arc2_tan[0] - r_vec2[0] * arc2_tan[1] > 0.0)
		std::swap(arc2->begin, arc2->end);
	if (arc2->end < arc2->begin) arc2->end += 2 * M_PI;
	
	// process cases where line should be used instead of arc
	if (center1[0] != center1[0]){
		arc1->radius = NAN;
		copy_point(arc1_point, arc1->center);
		Point vec;
		subtract_point(inflect, arc1_point, vec);
		arc1->begin = vec[0];
		arc1->end = vec[1];
	}
	if (center2[0] != center2[0]){
		arc2->radius = NAN;
		copy_point(arc2_point, arc1->center);
		Point vec;
		subtract_point(inflect, arc2_point, vec);
		arc2->begin = vec[0];
		arc2->end = vec[1];
	}
}