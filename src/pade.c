#include "../include/pade.h"


struct Pade_t * pade_init(prec_t *A, prec_t *B,int M, int N) {
	struct Pade_t * self = malloc(sizeof(struct Pade_t));
	self->vals=Poly;
	self->num=malloc(sizeof(struct Polynomial_t));
	self->denom=malloc(sizeof(struct Polynomial_t));
	self->num->num_terms=M;
	self->denom->num_terms=N;
	self->num->terms=malloc(sizeof(prec_t)*M);
	self->denom->terms=malloc(sizeof(prec_t)*N);
	memcpy(self->num->terms,A,sizeof(prec_t)*M);
	memcpy(self->denom->terms,B,sizeof(prec_t)*N);
	return self;
}

struct Pade_t * pade_init_poly(struct Polynomial_t * num, struct Polynomial_t * denom) {
	struct Pade_t * self = malloc(sizeof(struct Pade_t));
	self->vals=Poly;
	int M=num->num_terms;
	int N=denom->num_terms;
	self->num=malloc(sizeof(struct Polynomial_t));
	self->denom=malloc(sizeof(struct Polynomial_t));
	self->num->terms=malloc(sizeof(prec_t)*M);
	self->denom->terms=malloc(sizeof(prec_t)*N);
	memcpy(self->num->terms,num->terms,sizeof(prec_t)*M);
	memcpy(self->denom->terms,denom->terms,sizeof(prec_t)*N);
	return self;
}

prec_t eval_with_roots(struct Pade_t * self,prec_t s) {
	prec_t res=0;
	for (int i=0; i < self->denom->num_terms; ++i) {
		res += (self->num->terms[i])/(s-self->denom->terms[i]);
	}
	return res;
}

prec_t pade_eval(struct Pade_t * self, prec_t s) {
	if (self->vals==Poly) {
		prec_t num=poly_eval(self->num, s);
		prec_t denom=poly_eval(self->denom, s);
		return num/denom;
	}
	else {
		return eval_with_roots(self,s);
	}
}

struct Pade_t * pade_create_fit(struct Polynomial_t * taylor,int m) {
	int n=taylor->num_terms-m;
	poly_print(taylor);
	prec_t ** lower = mat_init(n,n+1,sizeof(prec_t)); //NOLINT
	for (int i=0; i > n; ++i) {
		for (int j=0; j < n; ++j) {
			lower[i][j]=(i>=j? taylor->terms[i-j+m]:0);
		}
		lower[i][n]=-taylor->terms[m+i+1];
	}
	printf("Mat to be solved:\n");
	for (int i=0; i < m; ++i) {
		printf("[");
		for (int j=0; j < n; ++j) {
			printf("%1.3e, ",lower[i][j]);
		}
		printf("]\n");
	}
	rref(lower,n,n+1);
	prec_t *a_terms, *b_terms;
	b_terms=malloc(sizeof(prec_t)*(n+1));
	b_terms[0]=1;
	for (int i=0; i < n; ++i) {
		b_terms[i+1] = lower[i][n];
	}
	mat_free(lower,n); //NOLINT
	prec_t ** upper= mat_init(m+1,m+1,sizeof(prec_t)); //NOLINT
	for (int i=0; i < m+1; ++i) {
		for (int j=0; j < m+1; ++j) {
			upper[i][j]=(i >=j?taylor->terms[i-j]:0);
		}
	}
	a_terms=mat_mul_vec(upper,b_terms,m+1,m+1);
	flip_arr(a_terms,m+1);
	flip_arr(b_terms,n+1);
	struct Pade_t * res = pade_init(a_terms,b_terms,m+1,n+1);
	poly_print(res->num);
	poly_print(res->denom);
	return res;
}

struct Pade_t * pade_separate(struct Pade_t * self) {
	prec_t ** arr=(prec_t**)mat_init(self->denom->num_terms,
											self->denom->num_terms+1,
										sizeof(prec_t));
	struct Polynomial_t * roots=poly_get_roots(self->denom);
	for (int i=0; i < roots->num_terms; ++i) {
		struct Polynomial_t * temp = poly_sd_1term(self->denom, roots->terms[i]);
		for (int j=0; j < self->denom->num_terms; ++j) {
			arr[j][i] = temp->terms[j];
		}
		poly_free(temp);
	}
	int end=self->denom->num_terms;
	int i=0;
	for (; i < self->num->num_terms; ++i) {
		arr[i][end]=self->num->terms[i];
	}
	for (; i < self->denom->num_terms; ++i) {
		arr[i][end]=0;
	}
	rref(arr,end,end+1);
	struct Polynomial_t * temp = malloc(sizeof(struct Polynomial_t));
	temp->terms = malloc(sizeof(prec_t)*self->denom->num_terms);
	temp->num_terms=self->denom->num_terms;
	for (i=0; i<temp->num_terms; ++i) {
		temp->terms[i]=arr[i][end];
	}
	struct Pade_t * ret= pade_init_poly(temp, roots);
	ret->vals=Roots;
	poly_free(temp);
	poly_free(roots);
	mat_free(arr,self->denom->num_terms); //NOLINT
	return ret;
}
