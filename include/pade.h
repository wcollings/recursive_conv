#ifndef __PADE_H__
#define __PADE_H__

#include "poly.h"
/*
 * Holds the terms for a divided fraction representation. Either in the form:
 *
 * \sum_{i=0}^{M} (a_i s^i)
 * -------------------------
 * \sum_{j=0}^{N} (b_j s^j)
 *
 * or
 *
 * \sum_{i=0}^{M} (a_i)/(s-\sigma_i)
 *
 */
struct Pade_t {
	struct Polynomial_t * num;
	struct Polynomial_t * denom;
};
struct Pade_t * Pade_init(double *A, double *B,int M, int N);
double Pade_eval(struct Pade_t * self, double s);
struct Pade_t * Pade_create_fit(double *C, int num_ele,int m);
#endif
