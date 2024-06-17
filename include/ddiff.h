#ifndef __DDIFF_H__
#define __DDIFF_H__

#include "central.h"

/*
 * Creates a divided difference table based on the x and y given.
 * This essentially approximates `n-1` derivatives of y(x).
 * out[i] -> the ith derivative of this function
 *
 * `x` The inputs to the function (must be `num_ele` elements long)
 * `y` the outputs of the function (must be `num_ele` elements long)
 * `num_ele` the number of elements in x and y. Also corresponds to greater than the number of derivatives you want to take
 * Returns:
 * a table of derivatives (`prec_t` s) of length `num_ele`
*/
prec_t * divided_difference(prec_t * x,prec_t *y, int num_ele);

/*
 * Take a list of derivatives, and turn that into a Taylor sequence.
 *
 * T(f(x)) = sum(i=0->N) (f^(i)(a)(x-a))/i!)
 */
void scale_to_taylor(prec_t * terms, int num_ele);
#endif
