#include "../include/pade.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct Pade_t * Pade_init(double *A, double *B,int M, int N) {
	struct Pade_t * self = malloc(sizeof(struct Pade_t));
	self->M=M;
	self->N=N;
	self->A=malloc(sizeof(double)*M);
	self->B=malloc(sizeof(double)*N);
	memcpy(self->A,A,sizeof(double)*M);
	memcpy(self->B,B,sizeof(double)*N);
	return self;
}
double Pade_eval(struct Pade_t * self, double s) {
	double num=self->A[0];
	for (int i=1; i < self->M; ++i) {
	 	num+=self->A[i]*pow(s,i);
	}
	double denom=self->B[0];
	for (int i=1; i < self->N; ++i) {
		denom+=self->B[i]*pow(s,i);
	}
	return num/denom;
}
struct Pade_t * Pade_create_fit(double *C, int num_ele,int m) {
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
