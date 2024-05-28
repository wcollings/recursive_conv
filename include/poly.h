#ifndef __POLY_H__
#define __POLY_H__

/*
 * Holds a polynomial, which is represented internally as simply an array of coefficients.
 */
struct Polynomial_t {
	int num_terms;
	double * terms; /* The terms are in highest-order to lowest-order*/
};

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
double poly_eval(struct Polynomial_t * a, double val);
void poly_print(struct Polynomial_t * p);
#endif
