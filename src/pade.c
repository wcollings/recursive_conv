#include "../include/pade.h"
#include <stdio.h>

void pade_iter(struct Pade_t * self) {
	prec_c_t temp;
	for (int i=0; i < self->num->num_terms; ++i)
		temp = self->num->terms[i];
	for (int i=0; i < self->denom->num_terms; ++i)
		temp = self->denom->terms[i];
}

struct Pade_t * pade_init(prec_c_t *A, prec_c_t *B,int M, int N) {
	struct Pade_t * self = malloc(sizeof(struct Pade_t));
	self->offset=0;
	self->vals=Vals;
	self->num=malloc(sizeof(struct Polynomial_t));
	self->denom=malloc(sizeof(struct Polynomial_t));
	self->num->num_terms=M;
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
	prec_c_t res=0;
	for (int i=0; i < self->denom->num_terms; ++i) {
		res += (self->num->terms[i])/(s-self->denom->terms[i]);
	}
	return res;
}

void pade_print(struct Pade_t * self) {
	char num[200];// = malloc(sizeof(char)*200);
	char denom[200];// = malloc(sizeof(char)*200);
	int num_ne=0, denom_ne=0;
	int i;
	for (i=0; i < self->num->num_terms-1; ++i) {
		int p = self->num->num_terms-i-1;
		int sz = snprintf(&num[num_ne],30,"+(%-"PRNT_SPEC"%+"PRNT_SPEC"i)x^%d",creal(self->num->terms[i]),cimag(self->num->terms[i]),p);
		if (sz > 0) num_ne+=sz;
	}
	int sz = snprintf(&num[num_ne],30,"+(%-"PRNT_SPEC"%+"PRNT_SPEC"i)",creal(self->num->terms[i]),cimag(self->num->terms[i]));
	/* int sz = snprintf(&num[num_ne],15,"%+1.4e",self->num->terms[self->num->num_terms-1]); */
	if (sz > 0) num_ne+=sz;
	for (i=0; i < self->denom->num_terms-1; ++i) {
		int p = self->denom->num_terms-i-1;
		int sz = snprintf(&denom[denom_ne],30,"(%-"PRNT_SPEC"%+"PRNT_SPEC"i)x^%d",creal(self->denom->terms[i]),cimag(self->denom->terms[i]),p);
		/* int sz = snprintf(&denom[denom_ne],15,"%+1.4ex^%d",self->denom->terms[i],p); */
		if (sz > 0) denom_ne+=sz;
	}
	sz = snprintf(&denom[denom_ne],30,"(%-"PRNT_SPEC"%+"PRNT_SPEC"i)",creal(self->denom->terms[i]),cimag(self->denom->terms[i]));
	/* sz = snprintf(&denom[denom_ne],15,"%+1.4e",self->denom->terms[self->denom->num_terms-1]); */
	if (sz > 0) denom_ne+=sz;
	int len = (denom_ne>num_ne?denom_ne:num_ne);
	char center[200];// = malloc(sizeof(char)*(len+6));
	snprintf(center,7,"r(x)= ");
	for (int i=0; i < len; ++i) {
		center[i+6]='-';
	}
	center[len+6]='\0';
	num[num_ne]='\0';
	denom[denom_ne]='\0';
	printf("      %s\n",num);
	printf("%s\n",center);
	printf("      %s\n",denom);
}

void centerline(int num) {
	for (int i=0; i<num; ++i) 
		putchar('-');
}
void pade_print_roots(struct Pade_t * self) {
	char num[300];
	char denom[300];
	char tn[30];
	char td[30];
	int nn, nd;
	printf("       \n");
	printf("r(x) = \n");
	printf("       ");
	printf("\e[2A");
	for (int i=0; i < self->num->num_terms; ++i) {
		printf("\e[1B");
		printf(" + ");
		printf("\e[1A");
		nn=snprintf(tn,30,"%-"PRNT_SPEC"%+"PRNT_SPEC"i",creal(self->num->terms[i]),cimag(self->num->terms[i]));
		nd=snprintf(td,30,"s%+"PRNT_SPEC"%+"PRNT_SPEC"i",-creal(self->denom->terms[i]),-cimag(self->denom->terms[i]));
		int len = (nn>nd? nn:nd);
		printf("%s",tn);
		printf("\e[1B\e[%dD",nn);
		centerline(len);
		printf("\e[1B\e[%dD",len);
		printf("%s",td);
		printf("\e[2A");
		fflush(stdout);
	}
	printf("\e[2B");
	printf("\n");
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

struct Pade_t * pade_create_fit(struct Polynomial_t * taylor,int m) {
	int n=taylor->num_terms-m-1;
	flip_arr(taylor->terms,taylor->num_terms);
	prec_c_t ** lower = mat_init(n,n+1,sizeof(prec_c_t));
	for (int i=0; i < n; ++i) {
		for (int j=0; j < n; ++j) {
			lower[i][j]=(i+m>=j? taylor->terms[i-j+m]:0);
		}
		lower[i][n]=-taylor->terms[m+i+1];
	}
	#if DEBUG_PRINTS
	printf("\nMat to be solved:\n");
	mat_print(lower,n,n+1);
	#endif
	rref(lower,n,n+1);
	prec_c_t *a_terms, *b_terms;
	b_terms=malloc(sizeof(prec_c_t)*(n+1));
	b_terms[0]=1;
	for (int i=0; i < n; ++i) {
		b_terms[i+1] = lower[i][n];
	}
	mat_free(lower,n);
	prec_c_t ** upper= mat_init(m+1,m+1,sizeof(prec_c_t));
	for (int i=0; i < m+1; ++i) {
		for (int j=0; j < m+1; ++j) {
			upper[i][j]=(i >= j?taylor->terms[i-j]:0);
		}
	}
	#if DEBUG_PRINTS
	printf("A-generating matrix:\n");
	mat_print(upper, m+1,m+1);
	#endif
	a_terms=mat_mul_vec(upper,b_terms,m+1,m+1);
	mat_free(upper,m+1);
	flip_arr(a_terms,m+1);
	flip_arr(b_terms,n+1);
	flip_arr(taylor->terms,taylor->num_terms);
	struct Pade_t * res = pade_init(a_terms,b_terms,m+1,n+1);
	return res;
}

struct Pade_t * pade_separate(struct Pade_t * self) {
	#if DEBUG_PRINTS==1
	printf("Pade_Separate\n");
	printf("given pade approximation:\n");
	pade_print(self);
	#endif
	prec_c_t ** arr=mat_init(self->denom->num_terms-1,
								    self->denom->num_terms,
								    sizeof(prec_c_t));
	struct Polynomial_t * roots=poly_get_roots(self->denom);
	int end=self->denom->num_terms-1;
	for (int i=0; i < roots->num_terms; ++i) {
		// This works as long as there are only 2 roots, i.e. if we're using P[1,2]
		arr[0][i]=1;
		arr[1][i]=roots->terms[i];
	}
	int i=0;
	for (; i < self->num->num_terms; ++i) {
		arr[i][end]=self->num->terms[i];
	}
	for (; i < self->denom->num_terms-1; ++i) {
		arr[i][end]=0;
	}
	rref(arr,end,end+1);
	struct Polynomial_t * temp = malloc(sizeof(struct Polynomial_t));
	temp->terms = malloc(sizeof(prec_c_t)*self->denom->num_terms);
	temp->num_terms=roots->num_terms;
	for (i=0; i<temp->num_terms; ++i) {
		temp->terms[i]=arr[i][end];
	}
	flip_arr(temp->terms,temp->num_terms);
	mat_free(arr,self->denom->num_terms-1); 
	struct Pade_t * ret=pade_init_poly(temp, roots);
	ret->vals=Roots;
	return ret;
}
