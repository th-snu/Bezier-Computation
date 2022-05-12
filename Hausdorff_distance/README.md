# Minimum Distance Computation of Two Cubic Bezier Curves
Tested on Ubuntu 20.04
OpenGL version: 4.6.0 NVIDIA 460.91.03

"make all" for compile

## Implementation

### Biarc Approximation
- Bezier curve is subdivided to power of given level by \+ , \- keyboard input
- Biarc approximation is used to approximate bezier curve
- If error bound is smaller with line approximation, line is used instead of a biarc

### Bounding Box
- Bounding box is computed based on biarc approximation
- Biarc is build only at the leaf nodes
- It is assumed that each biarc is approimxating corresponding subdivided bezier curve segment
- Minimum and maximum bound of two arcs combined is the resulting AABB
- Parent node is build by taking minimum and maximum value of each AABB

## Point to Bezier Curve Distance Computation
- Use heap with lower bound, and also update upper bounds
- Lines surrounding convex hull of the bezier curve is used to compute lower bound
- t = 0.5 point of bezier curve is used to update upper bound
- Lower bound is added to queue, and corresponding curve segment is divided until error is less than epsilon 

### Hausdorff Distance Computation
- Use heap with upper bound, t corresponding to two endpoints of corresponding bezier curve is stored.
- Use t to reconstruct corresponding bezier segment from original bezier curve
- Lower bound is computed by sampling points and computing point to bezier curve distance
- Upper bound is computed by projecting two end points of the bezier curve segment to the other bezier segment
- Distance between control points is used to compute upper bound

## Key Binding

- I : Reset points
- L : Use dotted line for bezier curves (Default: False)
- C : Draw control mesh (Default: True)

- 1 : Save current control points
- 2 : Load last saved control points

- Mouse Wheel : Change magnification of displayed curves
- Mouse Drag : Changes view area

- \+ , \- : Increase and Decrease subdivision level of the bezier curve (Default: 6)