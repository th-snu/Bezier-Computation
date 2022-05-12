#include "hausdorff.h"

REAL sample_points_distance(const Point &p, const CubicBezierCurve &c, int num_samples){
	REAL min_dist = std::numeric_limits<REAL>::max();
	for (int i = 0; i < num_samples + 1; i++){
		REAL t = (REAL)i / (REAL)num_samples;
		Point pt;
		evaluate(&c, t, pt);
		REAL dist = distance(p, pt);
		if (dist < min_dist)
			min_dist = dist;
	}
	return min_dist;
}

bool is_point_in_tri(const Point &p, const Point &t1, const Point &t2, const Point &t3){
	float d1, d2, d3;
	d1 = (p[0] - t2[0]) * (t1[1] - t2[1]) - (t1[0] - t2[0]) * (p[1] - t2[1]);
	d2 = (p[0] - t3[0]) * (t2[1] - t3[1]) - (t2[0] - t3[0]) * (p[1] - t3[1]);
	d3 = (p[0] - t1[0]) * (t3[1] - t1[1]) - (t3[0] - t1[0]) * (p[1] - t1[1]);
	return !(d1 < 0 || d2 < 0 || d3 < 0) || !(d1 > 0 || d2 > 0 || d3 > 0);
}

REAL distance_lower_bound(const Point &p, const CubicBezierCurve& c){
	// Find convex hull of bezier
	REAL x_min = std::numeric_limits<REAL>::max();
	REAL y_min = std::numeric_limits<REAL>::max();
	REAL x_max = -std::numeric_limits<REAL>::max();
	REAL y_max = -std::numeric_limits<REAL>::max();

	int indices[4];
	for (int i = 0; i < 4; i++){
		if (c.control_pts[i][0] < x_min){
			x_min = c.control_pts[i][0];
			indices[0] =  i;
		}
		if (c.control_pts[i][1] < y_min){
			y_min = c.control_pts[i][1];
			indices[1] =  i;
		}
		if (c.control_pts[i][0] > x_max){
			x_max = c.control_pts[i][0];
			indices[2] =  i;
		}
		if (c.control_pts[i][1] > y_max){
			y_max = c.control_pts[i][1];
			indices[3] =  i;
		}
	}
	
	int res = 1;

	std::vector<int> set_indices = {indices[0]};
 
    // Pick all elements one by one
    for (int i = 1; i < 4; i++) {
        int j = 0;
        for (j = 0; j < i; j++)
            if (indices[i] == indices[j]){
                break;
			}
        if (i == j){
			set_indices.push_back(indices[i]);
            res++;
		}
    }
	std::sort(set_indices.begin(), set_indices.end());

	std::vector<std::pair<REAL, REAL>> hull;

	if (res == 4){
		hull.push_back(std::make_pair(c.control_pts[indices[0]][0], c.control_pts[indices[0]][1]));
		hull.push_back(std::make_pair(c.control_pts[indices[1]][0], c.control_pts[indices[1]][1]));
		hull.push_back(std::make_pair(c.control_pts[indices[2]][0], c.control_pts[indices[2]][1]));
		hull.push_back(std::make_pair(c.control_pts[indices[3]][0], c.control_pts[indices[3]][1]));
	}
	else if (res == 3){
		// Check if other one point is within hull
		int missing_idx = 3;
		for (int i = 0; i < 3; i++){
			if (set_indices[i] != i){
				missing_idx = i;
				break;
			}
		}
		if (is_point_in_tri(c.control_pts[missing_idx], c.control_pts[set_indices[0]], c.control_pts[set_indices[1]], c.control_pts[set_indices[2]])){
			hull.push_back(std::make_pair(c.control_pts[set_indices[0]][0], c.control_pts[set_indices[0]][1]));
			hull.push_back(std::make_pair(c.control_pts[set_indices[1]][0], c.control_pts[set_indices[1]][1]));
			hull.push_back(std::make_pair(c.control_pts[set_indices[2]][0], c.control_pts[set_indices[2]][1]));
		}
		else {
			hull.push_back(std::make_pair(c.control_pts[set_indices[0]][0], c.control_pts[set_indices[0]][1]));
			hull.push_back(std::make_pair(c.control_pts[set_indices[1]][0], c.control_pts[set_indices[1]][1]));
			hull.push_back(std::make_pair(c.control_pts[set_indices[2]][0], c.control_pts[set_indices[2]][1]));
			hull.push_back(std::make_pair(c.control_pts[missing_idx][0], c.control_pts[missing_idx][1]));
		}
	}
	else if (res == 2){
		// Check if other two points are within hull
		std::vector<int> missing_idx;
        int ind1 = set_indices[0], ind2 = set_indices[1];
		for (int i = 0; i < 3; i++){
			if (set_indices[i] != i){
				missing_idx.push_back(i);
				if (missing_idx.size() == 1){
					set_indices.insert(set_indices.begin() + i, i);
				}
				else break;
			}
		}
        if (missing_idx.size() < 2) missing_idx.push_back(3);

		REAL dis1 = distance_line(c.control_pts[missing_idx[0]], c.control_pts[ind1], c.control_pts[ind2]);
		REAL dis2 = distance_line(c.control_pts[missing_idx[1]], c.control_pts[ind1], c.control_pts[ind2]);
    
		if (dis2 > dis1) std::swap(missing_idx[0], missing_idx[1]);
		if (is_point_in_tri(c.control_pts[missing_idx[1]], c.control_pts[missing_idx[0]], c.control_pts[ind1], c.control_pts[ind2])){
			hull.push_back(std::make_pair(c.control_pts[missing_idx[0]][0], c.control_pts[missing_idx[0]][1]));
			hull.push_back(std::make_pair(c.control_pts[ind1][0], c.control_pts[ind1][1]));
			hull.push_back(std::make_pair(c.control_pts[ind2][0], c.control_pts[ind2][1]));
		}
		else {
			hull.push_back(std::make_pair(c.control_pts[missing_idx[0]][0], c.control_pts[missing_idx[0]][1]));
			hull.push_back(std::make_pair(c.control_pts[missing_idx[1]][0], c.control_pts[missing_idx[1]][1]));
			hull.push_back(std::make_pair(c.control_pts[ind1][0], c.control_pts[ind1][1]));
			hull.push_back(std::make_pair(c.control_pts[ind2][0], c.control_pts[ind2][1]));
		}
	}
	else if (res == 1){
		Point pt = { c.control_pts[0][0], c.control_pts[0][1] };
		return distance(p, pt);
	}

	if (hull.size() == 4){
		// Order of points might be changed so that lines wouldn't intersect
        Point inter;
        Point line_tan1, line_tan2;
        Point line_end1 = {hull[1].first, hull[1].second};
        Point line_begin1 = {hull[0].first, hull[0].second};
        Point line_end2 = {hull[3].first, hull[3].second};
        Point line_begin2 = {hull[2].first, hull[2].second};

        subtract_point(line_end1, line_begin1, line_tan1);
        normalize(line_tan1);
        subtract_point(line_end2, line_begin2, line_tan2);
        normalize(line_tan2);
        get_line_intersection(line_tan1, line_begin1, line_tan2, line_begin2, inter);
        if (!((inter[0] - line_begin1[0]) * (inter[0] - line_end1[0]) >= 0 ||
            (inter[1] - line_begin1[1]) * (inter[1] - line_end1[1]) >= 0)){
            std::swap(hull[1].first, hull[2].first);
            std::swap(hull[1].second, hull[2].second);
        }

        line_begin1[0] = hull[2].first, line_begin1[1] = hull[2].second;
        line_begin2[0] = hull[0].first, line_begin2[1] = hull[0].second;
        subtract_point(line_end1, line_begin1, line_tan1);
        normalize(line_tan1);
        subtract_point(line_end2, line_begin2, line_tan2);
        normalize(line_tan2);
        get_line_intersection(line_tan1, line_begin1, line_tan2, line_begin2, inter);
        if (!((inter[0] - line_begin1[0]) * (inter[0] - line_end1[0]) >= 0 ||
            (inter[1] - line_begin1[1]) * (inter[1] - line_end1[1]) >= 0)){
            std::swap(hull[3].first, hull[2].first);
            std::swap(hull[3].second, hull[2].second);
        }
	}

	// Check if point is inside hull, use triangle check
	for (int i = 0; i < hull.size() - 1; i += 2){
		int j = (i + 1) % hull.size();
		int k = (i + 2) % hull.size();
		Point t1 = { hull[i].first, hull[i].second };
		Point t2 = { hull[j].first, hull[j].second };
		Point t3 = { hull[k].first, hull[k].second };
		if (is_point_in_tri(p, t1, t2, t3)){
			return 0.0;
        }
	}

	// If point is not in convex hull, compute distance to convex hull
	REAL min_dist = std::numeric_limits<REAL>::max();
	for (int i = 0; i < hull.size() - 1; i++){
		Point p1 = {hull[i].first, hull[i].second};
		Point p2 = {hull[i + 1].first, hull[i + 1].second};
		REAL dist = distance(p, p1, p2);
		if (dist < min_dist)
			min_dist = dist;
	}
	Point p1 = {hull[hull.size() - 1].first, hull[hull.size() - 1].second};
	Point p2 = {hull[0].first, hull[0].second};
	REAL dist = distance(p, p1, p2);
	if (dist < min_dist)
		min_dist = dist;

	return min_dist;
}

CubicBezierCurve subcurve_by_endpoint(const CubicBezierCurve &c, REAL t1, REAL t2){
	CubicBezierCurve seg1, seg2, dummy;
	subdivide(&c, t1, &dummy, &seg1);

	t2 = (t2 - t1) / (1.0 - t1);
	subdivide(&seg1, t2, &seg2, &dummy);

	return seg2;
}

REAL projection(const Point &p, const CubicBezierCurve &c){
	std::priority_queue<min_pair, std::vector<min_pair>, std::greater<min_pair>> q;
	
	// Use bounding box for bound computation, use biarc for final computation
	REAL lower_bound = distance_lower_bound(p, c);
	Point middle;
	evaluate(&c, 0.5, middle);
	REAL upper_bound = distance(p, middle);
	REAL eps = 1e-5;

	REAL globalt = 0.5;

	q.push(std::make_pair(lower_bound, std::make_pair(0.0, 1.0)));
	while (!q.empty()){
		REAL curr_bound = q.top().first;
		if (curr_bound > upper_bound){
			break;
		}
		auto t1 = q.top().second.first;
		auto t2 = q.top().second.second;
		q.pop();
		lower_bound = curr_bound;
		
		CubicBezierCurve curve = subcurve_by_endpoint(c, t1, t2);

		CubicBezierCurve curve1, curve2;
		subdivide(&curve, &curve1, &curve2);
		auto c1_interv = std::make_pair(t1, (t1 + t2)/2.0);
		auto c2_interv = std::make_pair((t1 + t2)/2.0 , t2);
		REAL c1_bound = distance_lower_bound(p, curve1);
        if (c1_bound < curr_bound) 
            c1_bound = curr_bound;
		REAL c2_bound = distance_lower_bound(p, curve2);
        if (c2_bound < curr_bound) 
            c2_bound = curr_bound;

		Point middle;
		evaluate(&curve1, 0.5, middle);
		REAL local_bound = distance(p, middle);
		if (upper_bound > local_bound) {
			upper_bound = local_bound;
			globalt = (c1_interv.first + c1_interv.second) / 2.0;
		}
        evaluate(&curve2, 0.5, middle);
        local_bound = distance(p, middle);
        if (upper_bound > local_bound) {
            upper_bound = local_bound;
            globalt = (c2_interv.first + c2_interv.second) / 2.0;
        }

		if (upper_bound - lower_bound < eps) break;

		if (c1_bound < upper_bound){
			q.push(std::make_pair(c1_bound, c1_interv));
        }
		if (c2_bound < upper_bound){
			q.push(std::make_pair(c2_bound, c2_interv));
        }
	}
    printf("%f %f\n", lower_bound, upper_bound);

	return globalt;
}

REAL sample_lower_bound(const CubicBezierCurve &c1, const CubicBezierCurve &c2, const int num_samples, Point &pt1, Point &pt2){
	std::vector<REAL> pts1x, pts1y;
	std::vector<REAL> pts2x, pts2y;
	for (int i = 0; i <= num_samples; i++){
		Point pt;
		const REAL t = (REAL)i / (REAL)num_samples;
		evaluate(&c1, t, pt);
		pts1x.push_back(pt[0]);
		pts1y.push_back(pt[1]);
	}

	REAL max_dist = 0.0;
	for (int i = 0; i < pts1x.size(); i++){
		Point pt = { pts1x[i], pts1y[i] };
		REAL t = projection(pt, c2);
		Point proj;
		evaluate(&c2, t, proj);
		REAL dist = distance(pt, proj);
		if (max_dist < dist){
            max_dist = dist;
            copy_point(pt, pt1);
            copy_point(proj, pt2);
        } 
	}

	return max_dist;
}