#include "biarc_approx.h"
#include "utils.h"
#include <limits>
#include <GL/glut.h>

void subdivide(const CubicBezierCurve *curve, CubicBezierCurve *output1, CubicBezierCurve *output2)
{
	subdivide(curve, 0.5, output1, output2);
}

void subdivide(const CubicBezierCurve *curve, REAL t, CubicBezierCurve *output1, CubicBezierCurve *output2)
{
	Point inter1[3];
	Point inter2[2];
	Point inter3[1];

	t = 1.0 - t;

	for(int i = 0; i < 3; i++){
		inter1[i][0] = curve->control_pts[i][0] * t + curve->control_pts[i+1][0] * (1.0-t);
		inter1[i][1] = curve->control_pts[i][1] * t + curve->control_pts[i+1][1] * (1.0-t);
	}
	for(int i = 0; i < 2; i++){
		inter2[i][0] = inter1[i][0] * t + inter1[i+1][0] * (1.0-t);
		inter2[i][1] = inter1[i][1] * t + inter1[i+1][1] * (1.0-t);
	}
	for(int i = 0; i < 1; i++){
		inter3[i][0] = inter2[i][0] * t + inter2[i+1][0] * (1.0-t);
		inter3[i][1] = inter2[i][1] * t + inter2[i+1][1] * (1.0-t);
	}

	copy_point(curve->control_pts[0], output1->control_pts[0]);
	copy_point(inter1[0], output1->control_pts[1]);
	copy_point(inter2[0], output1->control_pts[2]);
	copy_point(inter3[0], output1->control_pts[3]);

	copy_point(curve->control_pts[3], output2->control_pts[3]);
	copy_point(inter1[2], output2->control_pts[2]);
	copy_point(inter2[1], output2->control_pts[1]);
	copy_point(inter3[0], output2->control_pts[0]);
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

void get_biarc_inflect(const CubicBezierCurve *curve, Point &inflect)
{
	Point l1_center, l1_tan, l2_center, l2_tan, center;
	get_line_segs(curve, l1_tan, l1_center, l2_tan, l2_center);
	get_line_intersection(l1_tan, l1_center, l2_tan, l2_center, center);

	Point r_vec;
	subtract_point(center, curve->control_pts[0], r_vec);
	REAL r = norm(r_vec);

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
		copy_point(inflect, vec);
		arc1->begin = vec[0];
		arc1->end = vec[1];
	}
	if (center2[0] != center2[0]){
		arc2->radius = NAN;
		copy_point(arc2_point, arc2->center);	
		Point vec;
		copy_point(inflect, vec);
		arc2->begin = vec[0];
		arc2->end = vec[1];
	}
}

void draw_arc(const std::shared_ptr<Arc> arc){
	glColor3ub(255, 0, 0);
	glLineWidth(10.0);
	if (arc->is_line()){
		glBegin(GL_LINE_STRIP);
		glVertex2f(arc->center[0], arc->center[1]);
		glVertex2f(arc->center[0] + arc->begin, arc->center[1] + arc->end);
		glEnd();
	}
	else {
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i <= 100; i++)
		{
			Point pt;
			const REAL t = (REAL)i / (REAL)100;
			const REAL angle = t * (arc->end - arc->begin) + arc->begin;
			pt[0] = arc->center[0] + arc->radius * cos(angle);
			pt[1] = arc->center[1] + arc->radius * sin(angle);
			glVertex2f(pt[0], pt[1]);
		}
		glEnd();
	}
	glLineWidth(1.0);
}

REAL distance_line(const Point p, const Point line_begin, const Point line_end){
	Point vec;
	subtract_point(line_end, line_begin, vec);

	Point vec_perp;
	vec_perp[0] = -vec[1];
	vec_perp[1] = vec[0];

	Point vec_to_line;
	subtract_point(p, line_begin, vec_to_line);
	
	// Get distance by projecting A to B
	REAL vec_norm = norm(vec);
	REAL dis = (vec_to_line[0] * vec_perp[0] + vec_to_line[1] * vec_perp[1]) / vec_norm;

	return std::abs(dis);
}

REAL distance(const Point p, const Point line_begin, const Point line_end){
	Point vec;
	subtract_point(line_end, line_begin, vec);

	Point vec_perp;
	vec_perp[0] = -vec[1];
	vec_perp[1] = vec[0];

	Point vec_to_line;
	subtract_point(p, line_begin, vec_to_line);
	
	// Get distance by projecting A to B
	REAL vec_norm = norm(vec);
	REAL dis = (vec_to_line[0] * vec_perp[0] + vec_to_line[1] * vec_perp[1]) / vec_norm;

	// Check validity of distance
	Point perp_intersection = { p[0] - dis / vec_norm * vec_perp[0], p[1] - dis / vec_norm * vec_perp[1] };
	if ((perp_intersection[0] - line_begin[0]) * (perp_intersection[0] - line_end[0]) <= 0 &&
		(perp_intersection[1] - line_begin[1]) * (perp_intersection[1] - line_end[1]) <= 0){
		return abs(dis);
	}
	else {
		Point tmp;
		subtract_point(line_begin, p, tmp);
		REAL begin_dis = norm(tmp);
		subtract_point(line_end, p, tmp);
		REAL end_dis = norm(tmp);
		REAL ep_dis = std::min(begin_dis, end_dis);
		return ep_dis;
	}
}

REAL distance(const std::shared_ptr<Arc> arc1, const std::shared_ptr<Arc> arc2)
{
	REAL d = std::numeric_limits<REAL>::max();
	if (!arc1->is_line()){
		if (arc2->is_line()){
			return distance(arc2, arc1);
		}
		// Two arcs
		else {
			// Intersection check
			REAL center_d = distance(arc1->center, arc2->center);
			if (center_d < arc1->radius + arc2->radius
				&& center_d > std::abs(arc1->radius - arc2->radius)){

				// Find two candidate points, check if arcs contains them in a range
				REAL angle1c1, angle1c2;
				REAL angle2c1, angle2c2;
				Point tangent;
				subtract_point(arc1->center, arc2->center, tangent);
				REAL base_angle = atan2(tangent[1], tangent[0]);

				// Law of Cosines
				REAL offset_angle_1 = acos((arc1->radius * arc1->radius + center_d * center_d - arc2->radius * arc2->radius) / 2.0 / arc1->radius / center_d);
				REAL offset_angle_2 = acos((arc2->radius * arc2->radius + center_d * center_d - arc1->radius * arc1->radius) / 2.0 / arc2->radius / center_d);

				angle1c1 = base_angle - M_PI - offset_angle_1;
				angle1c2 = base_angle - M_PI + offset_angle_1;
				angle2c1 = base_angle + offset_angle_2 - 2 * M_PI;
				angle2c2 = base_angle - offset_angle_2 - 2 * M_PI;

				while (angle1c1 < arc1->begin) angle1c1 += 2 * M_PI;
				while (angle1c2 < arc1->begin) angle1c2 += 2 * M_PI;
				while (angle2c1 < arc2->begin) angle2c1 += 2 * M_PI;
				while (angle2c2 < arc2->begin) angle2c2 += 2 * M_PI;

				if ((angle1c1 < arc1->end && angle2c1 < arc2->end) ||
					(angle1c2 < arc1->end && angle2c2 < arc2->end)){
					draw_arc(arc1);
					draw_arc(arc2);
					return 0.0;
				}
			}

			// Check endpoint of two arcs
			Point arc1_e1 = { arc1->center[0] + arc1->radius * cos(arc1->begin), arc1->center[1] + arc1->radius * sin(arc1->begin) };
			Point arc1_e2 = { arc1->center[0] + arc1->radius * cos(arc1->end), arc1->center[1] + arc1->radius * sin(arc1->end) };
			Point arc2_e1 = { arc2->center[0] + arc2->radius * cos(arc2->begin), arc2->center[1] + arc2->radius * sin(arc2->begin) };
			Point arc2_e2 = { arc2->center[0] + arc2->radius * cos(arc2->end), arc2->center[1] + arc2->radius * sin(arc2->end) };
	
			d = std::min(d, distance(arc1_e1, arc2_e1));
			d = std::min(d, distance(arc1_e1, arc2_e2));
			d = std::min(d, distance(arc1_e2, arc2_e1));
			d = std::min(d, distance(arc1_e2, arc2_e2));

			// Check an endpoint of an arc and an interior point of the other arc, on a line connected to center
			Point da1e1, da1e2, da2e1, da2e2;
			subtract_point(arc2->center, arc1_e1, da1e1);
			subtract_point(arc2->center, arc1_e2, da1e2);
			subtract_point(arc1->center, arc2_e1, da2e1);
			subtract_point(arc1->center, arc2_e2, da2e2);

			normalize(da1e1);
			normalize(da1e2);
			normalize(da2e1);
			normalize(da2e2);

			REAL angle11, angle12, angle21, angle22;
			angle11 = std::atan2(da1e1[1], da1e1[0]);
			angle12 = std::atan2(da1e2[1], da1e2[0]);
			angle21 = std::atan2(da2e1[1], da2e1[0]);
			angle22 = std::atan2(da2e2[1], da2e2[0]);

			if (angle11 < arc2->begin) angle11 += 2 * M_PI;
			if (angle12 < arc2->begin) angle12 += 2 * M_PI;
			if (angle21 < arc1->begin) angle21 += 2 * M_PI;
			if (angle22 < arc1->begin) angle22 += 2 * M_PI;
			if (angle11 < arc2->end){
				Point interior = { arc2->center[0] + arc2->radius * cos(angle11), arc2->center[1] + arc2->radius * sin(angle11) };
				d = std::min(d, distance(interior, da1e1));
			}
			if (angle12 < arc2->end){
				Point interior = { arc2->center[0] + arc2->radius * cos(angle12), arc2->center[1] + arc2->radius * sin(angle12) };
				d = std::min(d, distance(interior, da1e2));
			}
			if (angle21 < arc1->end){
				Point interior = { arc1->center[0] + arc1->radius * cos(angle21), arc1->center[1] + arc1->radius * sin(angle21) };
				d = std::min(d, distance(interior, da2e1));
			}
			if (angle22 < arc1->end){
				Point interior = { arc1->center[0] + arc1->radius * cos(angle22), arc1->center[1] + arc1->radius * sin(angle22) };
				d = std::min(d, distance(interior, da2e2));
			}

			// Check two interior points on a line connecting two center points
			Point dc;
			subtract_point(arc1->center, arc2->center, dc);
			normalize(dc);
			REAL angle1 = std::atan2(-dc[1], -dc[0]);
			REAL angle2 = std::atan2(dc[1], dc[0]);

			if (angle1 < arc1->begin) angle1 += 2 * M_PI;
			if (angle2 < arc2->begin) angle2 += 2 * M_PI;
			if (angle1 < arc1->end && angle2 < arc2->end){
				Point interior1 = { arc1->center[0] + arc1->radius * cos(angle1), arc1->center[1] + arc1->radius * sin(angle1) };
				Point interior2 = { arc2->center[0] + arc2->radius * cos(angle2), arc2->center[1] + arc2->radius * sin(angle2) };
				d = std::min(d, distance(interior1, interior2));
			}

			return d;
		}
	}
	else {
		Point p1;
		p1[0] = arc1->begin;
		p1[1] = arc1->end;

		// An arc and a line
		if (!arc2->is_line()){
			Point &line_end = p1;
			Point &line_begin = arc1->center;

			// Intersection check
			REAL cld = distance(arc2->center, line_begin, line_end);
			if (cld < arc2->radius){
				Point vec;
				subtract_point(line_end, line_begin, vec);

				Point vec_perp;
				vec_perp[0] = -vec[1];
				vec_perp[1] = vec[0];

				Point vec_to_line;
				subtract_point(arc2->center, line_begin, vec_to_line);
				normalize(vec_to_line);
				normalize(vec_perp);
				
				REAL vec_norm = norm(vec);
				REAL dis = (vec_to_line[0] * vec_perp[0] + vec_to_line[1] * vec_perp[1]);
				vec_perp[0] *= dis;
				vec_perp[1] *= dis;

				REAL angle = std::atan2(vec_perp[1], vec_perp[0]);
				REAL offset_angle = std::acos(cld / arc2->radius);

				REAL angle1 = angle - offset_angle, angle2 = angle + offset_angle - 2 * M_PI;
				while (angle1 < arc2->begin) angle1 += 2 * M_PI;
				while (angle2 < arc2->begin) angle2 += 2 * M_PI;

				Point inter1 = { arc2->center[0] + arc2->radius * cos(angle1), arc2->center[1] + arc2->radius * sin(angle1) };
				Point inter2 = { arc2->center[0] + arc2->radius * cos(angle2), arc2->center[1] + arc2->radius * sin(angle2) };
				if ((angle1 < arc2->end && ((inter1[0] - line_begin[0]) * (inter1[0] - line_end[0]) < 0 || (inter1[1] - line_begin[1]) * (inter1[1] - line_end[1]) < 0)) ||
					(angle2 < arc2->end && ((inter2[0] - line_begin[0]) * (inter2[0] - line_end[0]) < 0 || (inter2[1] - line_begin[1]) * (inter2[1] - line_end[1]) < 0))){
					draw_arc(arc1);
					draw_arc(arc2);
					return 0;
				}

				return 0;
			}
			
			// Endpoint of arc to line
			Point arc2_e1 = { arc2->center[0] + arc2->radius * cos(arc2->begin), arc2->center[1] + arc2->radius * sin(arc2->begin) };
			Point arc2_e2 = { arc2->center[0] + arc2->radius * cos(arc2->end), arc2->center[1] + arc2->radius * sin(arc2->end) };

			REAL center_d = distance(arc2->center, line_begin, line_end);
			
			d = std::min(d, distance(arc2_e1, line_begin, line_end));
			d = std::min(d, distance(arc2_e2, line_begin, line_end));

			// Interior point of arc to line
			Point vec;
			subtract_point(line_end, line_begin, vec);

			Point vec_perp;
			vec_perp[0] = -vec[1];
			vec_perp[1] = vec[0];

			Point vec_to_line;
			subtract_point(arc2->center, line_begin, vec_to_line);
			
			REAL vec_norm = norm(vec);
			REAL dis = (vec_to_line[0] * vec_perp[0] + vec_to_line[1] * vec_perp[1]) > 0 ? 1.0 : -1.0;
			REAL angle = atan2(dis * vec_perp[1], dis * vec_perp[0]);

			if (angle < arc2->begin) angle += 2 * M_PI;
			if (angle < arc2->end){
				Point interior = { arc2->center[0] + dis * arc2->radius * vec_perp[0], arc2->center[1] + dis * arc2->radius * vec_perp[1] };
				d = std::min(d, distance(interior, line_begin, line_end));
			}

			return d;
		}
		// Two lines
		else {
			Point p2;
			p2[0] = arc2->begin;
			p2[1] = arc2->end;

			d = std::min(d, distance(arc1->center, arc2->center, p2));
			d = std::min(d, distance(p1, arc2->center, p2));
			d = std::min(d, distance(p2, arc1->center, p1));
			d = std::min(d, distance(arc2->center, arc1->center, p1));

			return d;
		}
	}
}