/*
 * @author Will
 * @brief Holds the representation of a polynomial, and all associated functions you could 
 * want with that
 *
 * @dependents pade
 *
*/
#ifndef __POLY_H__
#define __POLY_H__

#include "central.h"
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>

enum values_type_sel { Roots, Vals};

/*
 * Holds a polynomial, which is represented internally as simply an array of coefficients.
 *
 * The terms are ALWAYS in DECENDING order!
 */
struct Polynomial_t {
	int num_terms;
	enum values_type_sel tp;
	prec_c_t * terms; /* The terms are in highest-order to lowest-order*/
};


/*
 * Initialize a polynomial structure to hold a given number of elements. The elements need to then be added by you! They are _not_ malloc'd
*/
struct Polynomial_t * poly_init_bare(int num_terms);

/* 
 * Clear the memory and clean up the polynomial object
*/
void poly_free(struct Polynomial_t * a);

/*
 * Add two polynomials, return the result. The final polynomial is the same order as the higher of the two
 */
struct Polynomial_t * poly_add(struct Polynomial_t * a, struct Polynomial_t * b);

/*
 * multiply two polynomials, return the result. The final polynomial has order (|a| + |b| - 1)
 */
struct Polynomial_t * poly_mul(struct Polynomial_t * a, struct Polynomial_t * b);

/*
 * Evaluate a polynomial `a` at `val`
*/
prec_c_t poly_eval(struct Polynomial_t * a, prec_c_t val);

/*
 * Pretty print the polynomial
*/
void poly_print(struct Polynomial_t * p);

/*
 * perform sythetic division with a 1st order normalized polynomial (a.k.a. (x-a_0))
 *
 * The second polynomial (which we're dividing by) is just passed by the a_0 term. Assumes that 
 * a_0 is a root of the first polynomial as well, i.e. there will be no remainder. 
 */
struct Polynomial_t * poly_sd_1term(struct Polynomial_t * num, prec_c_t z);

/*
 * Use some root finding algorithm to find all the (real) roots of the given polynomial
 *
 * `NOTE`: the returned values of the polynomial are not coefficients! Call c_i self->values[i],
 * then the returned polynomial is encoded as (x-c_0)(x-c_1)...(x-c_n)
 */
struct Polynomial_t * poly_get_roots(struct Polynomial_t * num);

/*
 * Recenters a polynomial's representation around a given point.
 *
 * `NOTE`: should this just act on the original object, return void?
*/
struct Polynomial_t * poly_recenter(struct Polynomial_t * src, prec_c_t c);

/*
 * Turn a polynomial into a depressed polynomial, where the highest coefficient is 1
 */
void poly_depress(struct Polynomial_t * self);

/*
 * Convenience function to reverse an array
*/
void flip_arr(prec_c_t * arr, int n);
struct Polynomial_t * roots_binomial(struct Polynomial_t * self);
struct Polynomial_t * roots_trinomial(struct Polynomial_t * self);

#endif
