#ifndef __PADE_H__
#define __PADE_H__
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
	int M; /* The number of terms in the numerator */
	int N; /* The number of terms in the denominator */
	double *A; /* The numerator terms */
	double *B; /* The denominator coefficients */
};
struct Pade_t * Pade_init(double *A, double *B,int M, int N);
double Pade_eval(struct Pade_t * self, double s);
struct Pade_t * Pade_create_fit(double *C, int num_ele,int m);
#endif
