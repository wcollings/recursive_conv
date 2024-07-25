/*
 * @author Will
 * @brief Holds all the structs, functions, etc. for creating, evaluating, separating 
 * the Pade Approximant
 *
 * @dependents sara
*/
#ifndef __PADE_H__
#define __PADE_H__

#include "central.h"
#include "poly.h"
#include "linear.h"
#include "csv.h"
#include <stdlib.h>
#include <string.h> // just need this for memcpy
#include <math.h>

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
	prec_c_t offset;
	enum values_type_sel vals;
	struct Polynomial_t * num;
	struct Polynomial_t * denom;
};

/*
 * Initialize a Pade Approximant using two arrays. 
 * This object assumed ownership of the arrays
*/
struct Pade_t * pade_init(prec_c_t *A, prec_c_t *B,int M, int N);

/*
 *	If there's a K_0 term that we need to remember, this will do it
 * This object assumed ownership of the arrays
*/
struct Pade_t * pade_init_with_offset(prec_c_t *A, prec_c_t *B,int M, int N,prec_c_t offset);

void pade_free(struct Pade_t * self);
void pade_print(struct Pade_t * self);
void pade_print_roots(struct Pade_t * self);

/*
 * Initialize a Pade Approximant using two polynomials. 
 * This object assumes ownership of the polynomials
*/
struct Pade_t * pade_init_poly(struct Polynomial_t * num, struct Polynomial_t * denom);

/*
 * Evaluate a Pade Approximant at a certain point
*/
prec_c_t pade_eval(struct Pade_t * self, prec_c_t s);

/*
 * Create a Pade approximant given a polynomial that hold the coefficients
 * of a series representation of an object.
 * with N elements in the polynomial, we create a Pade approximant R_{m,N-m}
*/
struct Pade_t * pade_create_fit(struct Polynomial_t * taylor,int m);

/*
 * Perform Partial Fraction Decomposition on a given Pade Approximant
 *
 * `NOTE`: the returned Pade approximant holds two polynomials, but the 
 * polynomials are not representing straight coefficients. If we say a_i is self->num[i],
 * and b_i is self->denom[i], then these now represent (a_0)/(x-b_0) + (a_1)/(x-b_1) + ...
 *
 * This is also set in the field `enum valus_type_sel`, but it's good to note it yourself as well
 */
struct Pade_t * pade_separate(struct Pade_t * self);
#endif
