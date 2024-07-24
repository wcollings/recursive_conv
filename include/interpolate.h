#ifndef __INTERP_H__
#define __INTERP_H__

#include "central.h"

/*
 * Linear interpolation.
 * params:
 * `x1`: the earlier datapt (in time, or sweep variable of some kind)
 * `y1`: the output at x1
 * `x2`: the later datapt
 * `y2`: the output at x2
 * `xn`: the point to interpolate at
 * Returns: the output at xn, to a first approximation
*/
prec_t interpolate(prec_t x1, prec_t y1, prec_t x2, prec_t y2, prec_t xn);

#endif 
