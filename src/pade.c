#include "../include/pade.h"
#include "../include/poly.h"

struct Pade_t * pade_init(prec_c_t *A, prec_c_t *B,int M, int N) {
	struct Pade_t * self = malloc(sizeof(struct Pade_t));
	self->offset=0;
	self->vals=Vals;
	self->num=malloc(sizeof(struct Polynomial_t));
	self->denom=malloc(sizeof(struct Polynomial_t));
	self->num->coeff=1;
	self->denom->coeff=1;
	self->num->num_terms=M;
	self->num->tp=Vals;
	self->denom->tp=Vals;
	self->num->terms=A;
	self->denom->num_terms=N;
	self->denom->terms=B;
	return self;
}
struct Pade_t * pade_init_with_offset(prec_c_t *A, prec_c_t *B,int M, int N,prec_c_t offset) {
	struct Pade_t * self = pade_init(A,B,M,N);
	self->offset=offset;
	return self;
}

struct Pade_t * pade_init_poly(struct Polynomial_t * num, struct Polynomial_t * denom) {
	struct Pade_t * self = malloc(sizeof(struct Pade_t));
	self->num=num;
	self->denom=denom;
	self->vals=Vals;
	return self;
}

void pade_free(struct Pade_t * self) {
	poly_free(self->num);
	poly_free(self->denom);
	free(self);
}
prec_c_t eval_with_roots(struct Pade_t * self,prec_c_t s) {
	prec_c_t res=self->offset;
	for (int i=0; i < self->denom->num_terms; ++i) {
		res += (self->num->terms[i])/(s-self->denom->terms[i]);
	}
	return res;
}

prec_c_t pade_eval(struct Pade_t * self, prec_c_t s) {
	if (self->vals==Vals) {
		prec_c_t num=poly_eval(self->num, s);
		prec_c_t denom=poly_eval(self->denom, s);
		return self->offset + num/denom;
	}
	else {
		return self->offset + eval_with_roots(self,s);
	}
}

struct Pade_t * pade_separate(struct Pade_t * self) {
	struct Polynomial_t * roots=poly_get_roots(self->denom);
	int end=self->denom->num_terms-1;
	struct Polynomial_t * evaled = malloc(sizeof(struct Polynomial_t));
	evaled->terms = malloc(sizeof(prec_c_t)*self->denom->num_terms);
	evaled->num_terms=roots->num_terms;
	for (int i=0; i < roots->num_terms; ++i) {
		struct Polynomial_t * sd = poly_sd_1term(self->denom,roots->terms[i]);
		prec_c_t temp = poly_eval(sd,roots->terms[i]);
		prec_c_t up=poly_eval(self->num,roots->terms[i]);
		evaled->terms[i]=up/temp;
		poly_free(sd);
	}
	struct Pade_t * ret=pade_init_poly(evaled, roots);
	ret->vals=Roots;
	ret->offset=self->offset;
	return ret;
}
