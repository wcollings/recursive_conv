#include "../include/poly.h"

struct Polynomial_t * poly_init_bare(int num_terms) {
	struct Polynomial_t * out = malloc(sizeof(struct Polynomial_t));
	out->num_terms=num_terms;
	out->terms = calloc(num_terms,sizeof(double));
	return out;
}
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

struct Polynomial_t * poly_sd_1term(struct Polynomial_t * num, double z) {
	struct Polynomial_t * res=malloc(sizeof(struct Polynomial_t));
	res->num_terms=num->num_terms-1;
	res->terms=malloc(sizeof(double)*res->num_terms);
	res->terms[0]=1;
	for (int i=1; i < res->num_terms-1; ++i) {
		res->terms[i]=res->terms[i-1]*z+num->terms[i];
	}
	return res;
}

struct Polynomial_t * poly_get_roots(struct Polynomial_t * p) {
	struct Polynomial_t * roots = malloc(sizeof(struct Polynomial_t));
	roots->num_terms=p->num_terms;
	roots->terms=malloc(sizeof(double)*roots->num_terms);
	
	// find the roots

	return roots;
}
void poly_free(struct Polynomial_t * p) {
	free(p->terms);
	free(p);
}

void flip_arr(double * arr, int n) {
	double temp;
	for (int i=0; i < (int)n/2; ++i) {
		temp=arr[i];
		arr[i]=arr[n-i-1];
		arr[n-i-1]=temp;
	}
}
#define SUMUP(x) (int)(x*(x+1))/2

int * pascal(int row, int invert) {
	int num_ele = SUMUP(row+1)-1;
	int * tri = malloc(sizeof(int)*num_ele);
	int col=0;
	for (int i=0; i < num_ele; ++i) {
		int start=SUMUP(col), end = SUMUP(col+1)-1;
		if (i == start || i == end) {
			tri[i]=1;
			if (i==end) col++;
		}
		else {
			int offset=i-SUMUP(col);
			int idx=SUMUP(col-1)+offset-1;
			tri[i]=tri[idx]+tri[idx+1];
		}
		int offset=i-SUMUP(col);
		if (invert==1 && offset%2==0)
			tri[i]*=-1;
	}
	return tri;
}

struct Polynomial_t * poly_recenter(struct Polynomial_t * src, float c) {
	int N = src->num_terms;
	int * tri = pascal(N, 1);
	struct Polynomial_t * dest = poly_init_bare(N);
	for (int i=0; i < N; ++i) {
		double term = 0;
		for (int j = 0; j <= i; ++i) {
			term += tri[SUMUP(N-j-1)+(i-j)]*src->terms[j]*pow(c,i-j);
		}
		dest->terms[i] = term;
		term=0;
	}
	return dest;
}

void poly_depress(struct Polynomial_t * self) {
	double max_term = self->terms[0];
	for (int i=0; i < self->num_terms; ++i) {
		self->terms[i]/=max_term;
	}
}
