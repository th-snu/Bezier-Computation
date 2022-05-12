#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

enum XY { X = 0, Y };

#define SET_VECTOR2(V, V1, V2)          do { (V)[X] = (V1); (V)[Y] = (V2); } while (0)
#define COPY_PT(DST, SRC)               do { (DST)[X] = SRC[X]; (DST)[Y] = SRC[Y]; } while (0)
#define VECTOR2_X_SCALA_ADD(O, V, S)    do { O[X] += (S) * (V)[X]; O[Y] += (S) * (V)[Y]; } while (0)