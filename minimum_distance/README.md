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

### Minimum Distance Computation
- Pair of segments of two bezier curves are added to priority queue
- Pair with smallest lower bound is updated and segment corresponding to AABB with larger volume is divided
- Computation of distance between two subdivided curves at leaf is approximated using biarcs
- Distance between biarcs updates upper bound
- First encountered intersecting biarcs are visualized with red lines
- Pair of curves that updated the upper bound for the last time is visualized with red lines as well

## Key Binding

- I : Reset points
- L : Use dotted line for bezier curves (Default: False)
- C : Draw control mesh (Default: True)

- 1 : Draw line or arc used for AABB approximation (Default: False)
- 2 : Draw AABB for the bezier curve (Default: True)
- 3 : Draw entire hierarchy of AABB (Default: False)

- Mouse Wheel : Change magnification of displayed curves
- Mouse Drag : Changes view area

- \+ , \- : Increase and Decrease subdivision level of the bezier curve (Default: 6)