# Bezier to Biarcs approximator
Program is compiled in Ubuntu 20.04 using make all command

## Implementation
- Bezier curve is subdivided to power of given level by + - keyboard input
- Biarc approximation is used to approximate bezier curve
- Bounding box is computed based on biarc approximation
- It is assumed that each biarc is approimxating corresponding subdivided bezier curve segment
- Minimum and maximum bound of two arcs combined is the resulting AABB
- If arc has wider than 90 degree angle, line is used to approximate the curve instead

## Key Binding

- I : Reset points
- L : Use dotted line for bezier curves
- C : Draw control mesh

- 1 : Draw line or arc used for AABB approximation
- 2 : Draw AABB for the bezier curve
- 3 : Draw entire hierarchy of AABB

- \+ , \- : Increase and Decrease subdivision level of the bezier curve