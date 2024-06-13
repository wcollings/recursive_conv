#ifndef __DDIFF_H__
#define __DDIFF_H__

/*
 * @brief Creates a divided difference table based on the x and y given.
 * This essentially approximates `n-1` derivatives of y(x).
 * out[i] -> the ith derivative of this function
 *
 * `x` The inputs to the function (must be `num_ele` elements long)
 * `y` the outputs of the function (must be `num_ele` elements long)
 * `num_ele` the number of elements in x and y. Also corresponds to greater than the number of derivatives you want to take
 * Returns:
 * a table of derivatives (doubles) of length `num_ele`
*/
double * divided_difference(double * x,double *y, int num_ele);

/*
 * Take a list of derivatives, and turn that into a Taylor sequence.
 *
 * T(f(x)) = sum(i=0->N) (f^(i)(a)(x-a))/i!)
 */
void scale_to_taylor(double * terms, int num_ele);
#endif
