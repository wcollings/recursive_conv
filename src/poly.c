#include "poly.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

struct Polynomial_t * poly_add(struct Polynomial_t * a, struct Polynomial_t * b) {
	struct Polynomial_t * res=malloc(sizeof(struct Polynomial_t));
	// a should always be the bigger of the two
	if (a->num_terms < b->num_terms) {
		struct Polynomial_t * lhs=a;
		a=b;
		b=lhs;
	}
	res->num_terms=a->num_terms;
	res->terms=malloc(sizeof(double)*res->num_terms);
	int i=0;
	for (; i < b->num_terms; ++i) {
		res->terms[i]=a->terms[i]+b->terms[i];
	}
	for (; i < a->num_terms; ++i) {
		res->terms[i]=a->terms[i];
	}
	return res;
}
struct Polynomial_t * poly_mul(struct Polynomial_t * a, struct Polynomial_t * b) {
	struct Polynomial_t * res=malloc(sizeof(struct Polynomial_t));
	// a should always be the bigger of the two
	if (a->num_terms < b->num_terms) {
		struct Polynomial_t * lhs=a;
		a=b;
		b=lhs;
	}
	res->num_terms=a->num_terms+b->num_terms-1;
	res->terms=calloc(res->num_terms,sizeof(double));
	for (int i=0; i < a->num_terms; ++i) {
		for (int j=0; j < b->num_terms; ++j ) {
			res->terms[i+j]=a->terms[i]*b->terms[i];
		}
	}
	return res;
}
double poly_eval(struct Polynomial_t * a, double val) {
	double res=a->terms[0];
	for (int i=1; i < a->num_terms; ++i) {
		res=a->terms[i]+val*res;
	}
	return res;
}
void poly_print(struct Polynomial_t * p) {
	printf("%-ex^%d",p->terms[0],p->num_terms);
	for (int i=1; i < p->num_terms; ++i) {
		printf("%+ex^%d",p->terms[i],p->num_terms-i);
	}
}
