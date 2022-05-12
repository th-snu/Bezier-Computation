#include <queue>
#include <limits>
#include <bits/stdc++.h>
#include "aabb.h"
#include "biarc_approx.h"

REAL sample_points_distance(const Point &p, const CubicBezierCurve &c, int num_samples);

bool is_point_in_tri(const Point &p, const Point &t1, const Point &t2, const Point &t3);

REAL distance_lower_bound(const Point &p, const CubicBezierCurve& c);

CubicBezierCurve subcurve_by_endpoint(const CubicBezierCurve &c, REAL t1, REAL t2);

typedef std::pair<REAL, std::pair<REAL, REAL>> min_pair;

REAL projection(const Point &p, const CubicBezierCurve &c);

REAL sample_lower_bound(const CubicBezierCurve &c1, const CubicBezierCurve &c2, const int num_samples, Point &pt1, Point &pt2);
