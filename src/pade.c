#include "../include/pade.h"
#include "../include/poly.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct Pade_t * Pade_init(double *A, double *B,int M, int N) {
	struct Pade_t * self = malloc(sizeof(struct Pade_t));
	self->num=malloc(sizeof(struct Polynomial_t));
	self->denom=malloc(sizeof(struct Polynomial_t));
	self->num->num_terms=M;
	self->denom->num_terms=N;
	self->num->terms=malloc(sizeof(double)*M);
	self->denom->terms=malloc(sizeof(double)*N);
	memcpy(self->num->terms,A,sizeof(double)*M);
	memcpy(self->denom->terms,B,sizeof(double)*N);
	return self;
}
double Pade_eval(struct Pade_t * self, double s) {
	double num=poly_eval(self->num, s);
	double denom=poly_eval(self->denom, s);
	return num/denom;
}
struct Pade_t * Pade_create_fit(double *C, int num_ele,int m) {
	// all the logic for matrix creation has already been worked out in pade.py
	// Just reference that
	float ** A=malloc(sizeof(float*)*m);
	for (int i=0; i < m; ++i) A[i]=malloc(sizeof(float)*(m+2));
	for (int i=num_ele-1; i > m-1; --i) {
		for (int j=0; j < m+1; ++j) {
			A[m-i][j]=(i>=j? C[i-j]:0);
		}
	}
	double *a_terms, *b_terms;
	return Pade_init(a_terms, b_terms, m, num_ele-m);
}
