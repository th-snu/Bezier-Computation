#include <GL/glut.h>
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
    std::shared_ptr<Hierarchy> left = nullptr;
    std::shared_ptr<Hierarchy> right = nullptr;
    std::shared_ptr<Arc> arc = nullptr;
};

AABB get_arc_aabb(const Arc *arc);

CubicBezierCurve to_bezier(const Arc *arc);

REAL arc_approx_error_bound(const Arc *arc, const CubicBezierCurve *curve);

REAL bezier_error_bound(const CubicBezierCurve *curve1, const CubicBezierCurve *curve2);

AABB combine(const AABB &box1, const AABB &box2);

REAL get_AABB(const CubicBezierCurve &seg, const Arc arc, AABB& box, bool isDraw);

void draw_AABB(const AABB &box);

REAL distance(const AABB &box1, const AABB &box2);

REAL volume(const AABB &box);
