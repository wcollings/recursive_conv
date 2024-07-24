#include "../include/interpolate.h"

prec_t interpolate(prec_t x1, prec_t y1, prec_t x2, prec_t y2, prec_t xn) {
	prec_t slope = (y2-y1)/(x2-x1);
	return xn*slope + y1;
}
