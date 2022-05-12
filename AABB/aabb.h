#include "curve.h"
#include <memory>

class AABB {
public:
    REAL x[2];
    REAL y[2];
};

class Hierarchy {
public:
    CubicBezierCurve curve;
    AABB box;
    std::shared_ptr<Hierarchy> left;
    std::shared_ptr<Hierarchy> right;
};

AABB get_arc_aabb(const Arc *arc);

CubicBezierCurve to_bezier(const Arc *arc);

REAL arc_approx_error_bound(const Arc *arc, const CubicBezierCurve *curve);

REAL bezier_error_bound(const CubicBezierCurve *curve1, const CubicBezierCurve *curve2);

AABB combine(const AABB &box1, const AABB &box2);