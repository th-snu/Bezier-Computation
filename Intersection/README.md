# Intersection Test of Two Cubic Bezier Curves
Tested on Ubuntu 20.04
OpenGL version: 4.6.0 NVIDIA 460.91.03

"make all" for compile

## Implementation

### Biarc Approximation
- Bezier curve is subdivided to power of given level by \+ , \- keyboard input
- Biarc approximation is used to approximate bezier curve
- If arc has wider than 90 degree angle, line is used to approximate the curve instead

### Bounding Box
- Bounding box is computed based on biarc approximation
- It is assumed that each biarc is approimxating corresponding subdivided bezier curve segment
- Minimum and maximum bound of two arcs combined is the resulting AABB
- Parent node is build by taking minimum and maximum value of each AABB

### Intersection Test
- Intersection tested is done based on BVH tree built above
- Thus, result may not be precise if bezier curve is not subdivided enough
- Intersection of two curves are drawn with thick red lines

## Key Binding

- I : Reset points
- L : Use dotted line for bezier curves (Default: False)
- C : Draw control mesh (Default: True)

- 1 : Draw line or arc used for AABB approximation (Default: False)
- 2 : Draw AABB for the bezier curve (Default: True)
- 3 : Draw entire hierarchy of AABB (Default: False)

- \+ , \- : Increase and Decrease subdivision level of the bezier curve (Default: 6)