#include "curve.h"
#include <memory>

void subdivide(const CubicBezierCurve *curve, CubicBezierCurve *output1, CubicBezierCurve *output2);

void subdivide(const CubicBezierCurve *curve, REAL t, CubicBezierCurve *output1, CubicBezierCurve *output2);

void subdivide(const CubicBezierCurve *curve, std::vector<CubicBezierCurve> &output, int power);

void get_tangent(const CubicBezierCurve *curve, Point &tan_begin, Point &tan_end);

void get_line_segs(const CubicBezierCurve *curve, Point &l1_tan, Point &l1_center, Point &l2_tan, Point &l2_center);

void set_arc_center(const CubicBezierCurve *curve, const Point &arc1_point, const Point &arc2_point, const Point &inflect, Arc *arc1, Arc *arc2);

void get_biarc_inflect(const CubicBezierCurve *curve, Point &inflect);

void to_biarc(const CubicBezierCurve *curve, Point &inflect, Arc *arc1, Arc *arc2);

REAL distance(const std::shared_ptr<Arc> arc1, const std::shared_ptr<Arc> arc2);

void draw_arc(const std::shared_ptr<Arc> arc);
