#include "../include/pade.h"

void mat_free(void ** A, int m) {
	for (int i=0; i < m; ++i ) {
		free(A[i]);
	}
	free(A);
}

void ** mat_init(int m, int n,size_t size) {
	void ** inner = malloc(sizeof(void*)*m);
	for (int i=0; i < m; ++i) {
		inner[i] = malloc(size*n);
	}
	return inner;
}
struct Pade_t * pade_init(double *A, double *B,int M, int N) {
	struct Pade_t * self = malloc(sizeof(struct Pade_t));
	self->vals=Poly;
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
struct Pade_t * pade_init_poly(struct Polynomial_t * num, struct Polynomial_t * denom) {
	struct Pade_t * self = malloc(sizeof(struct Pade_t));
	self->vals=Poly;
	int M=num->num_terms;
	int N=denom->num_terms;
	self->num=malloc(sizeof(struct Polynomial_t));
	self->denom=malloc(sizeof(struct Polynomial_t));
	self->num->terms=malloc(sizeof(double)*M);
	self->denom->terms=malloc(sizeof(double)*N);
	memcpy(self->num->terms,num->terms,sizeof(double)*M);
	memcpy(self->denom->terms,denom->terms,sizeof(double)*N);
	return self;
}

double eval_with_roots(struct Pade_t * self,double s) {
	double res=0;
	for (int i=0; i < self->denom->num_terms; ++i) {
		res += (self->num->terms[i])/(s-self->denom->terms[i]);
	}
	return res;
}
double pade_eval(struct Pade_t * self, double s) {
	if (self->vals==Poly) {
		double num=poly_eval(self->num, s);
		double denom=poly_eval(self->denom, s);
		return num/denom;
	}
	else {
		return eval_with_roots(self,s);
	}
}
struct Pade_t * pade_create_fit(struct Polynomial_t * taylor,int m) {
	int n=taylor->num_terms-m;
	double ** lower = mat_init(n,n+1,sizeof(double)); //NOLINT
	for (int i=0; i > n; ++i) {
		for (int j=0; j < n; ++j) {
			lower[i][j]=(i>=j? taylor->terms[i-j+m]:0);
		}
		lower[i][n]=-taylor->terms[m+i+1];
	}
	rref(lower,n+1,n);
	double *a_terms, *b_terms;
	b_terms=malloc(sizeof(double)*(n+1));
	b_terms[0]=1;
	for (int i=0; i < n; ++i) {
		b_terms[i+1] = lower[i][n];
	}
	mat_free(lower,n); //NOLINT
	double ** upper= mat_init(m+1,m+1,sizeof(double)); //NOLINT
	for (int i=0; i < m+1; ++i) {
		for (int j=0; j < m+1; ++j) {
			upper[i][j]=(i >=j?taylor->terms[i-j]:0);
		}
	}
	a_terms=mat_mul_vec(upper,b_terms,m+1,m+1);
	flip_arr(a_terms,m+1);
	flip_arr(b_terms,n+1);
	return pade_init(a_terms, b_terms, m+1, n+1);
}

struct Pade_t * pade_separate(struct Pade_t * self) {
	double ** arr=(double**)mat_init(self->denom->num_terms,
											self->denom->num_terms+1,
										sizeof(double));
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
	temp->terms = malloc(sizeof(double)*self->denom->num_terms);
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
