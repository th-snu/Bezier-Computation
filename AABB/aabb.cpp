#include "aabb.h"
#include "utils.h"
#include <algorithm>

AABB get_arc_aabb(const Arc *arc){
    REAL begin, end;
    begin = arc->begin < 0 ? arc->begin + 2 * M_PI : arc->begin;
    end = arc->begin < 0 ? arc->end + 2 * M_PI : arc->end;
    
    std::vector<REAL> cand_x;
    std::vector<REAL> cand_y;

    if (arc->radius != arc->radius){
        cand_x.push_back(arc->center[0]);
        cand_x.push_back(arc->begin);
        cand_y.push_back(arc->center[1]);
        cand_y.push_back(arc->end);
    }
    else {
        if (begin < 0.5 * M_PI && end > 0.5 * M_PI)
            cand_y.push_back(arc->center[1] + arc->radius);
        if (begin < 1.0 * M_PI && end > 1.0 * M_PI)
            cand_x.push_back(arc->center[0] - arc->radius);
        if (begin < 1.5 * M_PI && end > 1.5 * M_PI)
            cand_y.push_back(arc->center[1] - arc->radius);
        if (begin < 2.0 * M_PI && end > 2.0 * M_PI)
            cand_x.push_back(arc->center[0] + arc->radius);

        cand_x.push_back(arc->center[0] + arc->radius * cos(arc->begin));
        cand_x.push_back(arc->center[0] + arc->radius * cos(arc->end));
        cand_y.push_back(arc->center[1] + arc->radius * sin(arc->begin));
        cand_y.push_back(arc->center[1] + arc->radius * sin(arc->end));
    }

    AABB res;
    auto minmax_x = std::minmax_element(cand_x.begin(), cand_x.end());
    auto minmax_y = std::minmax_element(cand_y.begin(), cand_y.end());

    SET_VECTOR2(res.x, *(minmax_x.first), *(minmax_x.second));
    SET_VECTOR2(res.y, *(minmax_y.first), *(minmax_y.second));

    return res;
}

CubicBezierCurve to_bezier(const Arc *arc){
    CubicBezierCurve bezier;
    Point *p;
    p = (bezier.control_pts);
    if (arc->radius != arc->radius){
        SET_VECTOR2(p[0], arc->center[0], arc->center[1]);
        SET_VECTOR2(p[3], arc->begin, arc->end);
        middle_point(p[0], p[3], p[1]);
        middle_point(p[1], p[3], p[2]);
        middle_point(p[0], p[1], p[1]);
    }
    else {
        Point arc_begin, arc_end;
        SET_VECTOR2(arc_begin, arc->center[0] + arc->radius * cos(arc->begin), arc->center[1] + arc->radius * sin(arc->begin));
        SET_VECTOR2(arc_end, arc->center[0] + arc->radius * cos(arc->end), arc->center[1] + arc->radius * sin(arc->end));

        Point tan1, tan2;
        SET_VECTOR2(tan1, sin(arc->begin), -cos(arc->begin));
        SET_VECTOR2(tan2, sin(arc->end), -cos(arc->end));

        Point intersection;
        get_line_intersection(tan1, arc_begin, tan2, arc_end, intersection);

        if (IS_NAN(intersection)) {
            // Control should not reach here
        }
        else {
            copy_point(arc_begin, p[0]);
            copy_point(arc_end, p[3]);

            SET_VECTOR2(p[1], (p[0][0] + intersection[0] * 2.0)/3.0, (p[0][1] + intersection[1] * 2.0)/3.0);
            SET_VECTOR2(p[2], (p[3][0] + intersection[0] * 2.0)/3.0, (p[3][1] + intersection[1] * 2.0)/3.0);
        }
    }

    return bezier;
}

REAL arc_approx_error_bound(const Arc *arc, const CubicBezierCurve *curve){
    CubicBezierCurve approx_arc = to_bezier(arc);
    REAL inter_bezier_error = bezier_error_bound(curve, &approx_arc);

    // No approxmiation error for a line
    if (IS_NAN(arc->radius)){
        return inter_bezier_error;
    }

    REAL mid_angle = (arc->begin + arc->end) / 2.0;
    Point e_q, e_o;
    evaluate(&approx_arc, 0.5, e_q);
    SET_VECTOR2(e_o, arc->center[0] + arc->radius * cos(mid_angle), arc->center[1] + arc->radius * sin(mid_angle));

    REAL arc_bezier_error = distance(e_q, e_o);

    return arc_bezier_error + inter_bezier_error;
}

REAL bezier_error_bound(const CubicBezierCurve *curve1, const CubicBezierCurve *curve2){
    REAL bound1 = 0.0, bound2 = 0.0;
    for (int i = 0; i < 4; i++){
        REAL cand = distance(curve1->control_pts[i], curve2->control_pts[i]);
        if (cand > bound1) bound1 = cand;
    }
    for (int i = 0; i < 4; i++){
        REAL cand = distance(curve1->control_pts[3-i], curve2->control_pts[i]);
        if (cand > bound2) bound2 = cand;
    }

    return bound1 > bound2 ? bound2 : bound1;
}

AABB combine(const AABB &box1, const AABB &box2){
    AABB res;
    res.x[0] = std::min(box1.x[0], box2.x[0]);
    res.y[0] = std::min(box1.y[0], box2.y[0]);
    res.x[1] = std::max(box1.x[1], box2.x[1]);
    res.y[1] = std::max(box1.y[1], box2.y[1]);

    return res;
}