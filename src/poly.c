#include "../include/poly.h"
#include <math.h>

#define SUMUP(x) (int)(x*(x+1))/2

struct Polynomial_t * poly_init_bare(int num_terms) {
	struct Polynomial_t * out = malloc(sizeof(struct Polynomial_t));
	out->num_terms=num_terms;
	out->terms = calloc(num_terms,sizeof(prec_t));
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
	res->terms=malloc(sizeof(prec_t)*res->num_terms);
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
	res->terms=calloc(res->num_terms,sizeof(prec_t));
	for (int i=0; i < a->num_terms; ++i) {
		for (int j=0; j < b->num_terms; ++j ) {
			res->terms[i+j]=a->terms[i]*b->terms[i];
		}
	}
	return res;
}
prec_t poly_eval(struct Polynomial_t * a, prec_t val) {
	prec_t res=a->terms[0];
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

struct Polynomial_t * poly_sd_1term(struct Polynomial_t * num, prec_t z) {
	struct Polynomial_t * res=malloc(sizeof(struct Polynomial_t));
	res->num_terms=num->num_terms-1;
	res->terms=malloc(sizeof(prec_t)*res->num_terms);
	res->terms[0]=1;
	for (int i=1; i < res->num_terms-1; ++i) {
		res->terms[i]=res->terms[i-1]*z+num->terms[i];
	}
	return res;
}

struct Polynomial_t * poly_get_roots(struct Polynomial_t * p) {
	struct Polynomial_t * roots = malloc(sizeof(struct Polynomial_t));
	roots->num_terms=p->num_terms;
	roots->terms=malloc(sizeof(prec_t)*roots->num_terms);
	
	// find the roots

	return roots;
}
void poly_free(struct Polynomial_t * p) {
	free(p->terms);
	free(p);
}

void flip_arr(prec_t * arr, int n) {
	prec_t temp;
	for (int i=0; i < (int)n/2; ++i) {
		temp=arr[i];
		arr[i]=arr[n-i-1];
		arr[n-i-1]=temp;
	}
}

/*
 * Creates and returns an array containing the coefficients of pascal's triangle.
 * This is used primarily in recentering a polynomial. Note that this is a 1d array representing
 * a ragged 2d array. Each row is one element longer than the last. The returned array looks like:
 *
 * [1 | 1 1 | 1 2 1 | 1 3 3 1 |...]
 *
 * If the `invert` flag is set to 1, every even element of each row is inverted. So the result will instead be:
 *
 * [ 1 | 1 -1 | 1 -2 1 | 1 -3 3 -1 |...]
 *
*/
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

struct Polynomial_t * poly_recenter(struct Polynomial_t * src, prec_t c) {
	int N = src->num_terms;
	int * tri = pascal(N, 1);
	struct Polynomial_t * dest = poly_init_bare(N);
	for (int i=0; i < N; ++i) {
		prec_t term = 0;
		for (int j = 0; j <= i; ++i) {
			term += tri[SUMUP(N-j-1)+(i-j)]*src->terms[j]*pow(c,i-j); //NOLINT
		}
		dest->terms[i] = term;
		term=0;
	}
	return dest;
}

void poly_depress(struct Polynomial_t * self) {
	prec_t max_term = self->terms[0];
	for (int i=0; i < self->num_terms; ++i) {
		self->terms[i]/=max_term;
	}
}

struct Polynomial_t * poly_diff(struct Polynomial_t * self) {
	struct Polynomial_t * out = poly_init_bare(self->num_terms-1);
	prec_t term = self->num_terms-1;
	for (int i=0; i < term; ++i) {
		out->terms[i] = self->terms[i]*(term-i);
	}
	return out;
}

struct Polynomial_t * roots_quartic(struct Polynomial_t * self) {
	
}

struct Polynomial_t * roots_binomial(struct Polynomial_t * self) {
	prec_t b,c;
	struct Polynomial_t * roots = poly_init_bare(2);
	if (self->terms[0] != 1) {
		poly_depress(self);
	}
	b=self->terms[1];
	c=self->terms[2];
	if (b > 0) {
		if (pow(b,2) > c) {
			roots->terms[0] = -b-sqrt(pow(b,2)-c);
			roots->terms[1] = c/roots->terms[0];
		}
		// if the roots are complex
		// TODO: handle this more properly, once I decide what that means
		else {
			roots->terms[0] = INFINITY;
			roots->terms[1] = INFINITY;
		}
	}
	else {
		if (pow(b,2) > c) {
			roots->terms[0] = -b+sqrt(pow(b,2)-c);
			roots->terms[1] = c/roots->terms[0];

		}
		// if the roots are complex
		else {
			roots->terms[0] = INFINITY;
			roots->terms[1] = INFINITY;
		}
	}
	return roots;
}

prec_t newton(struct Polynomial_t * f, prec_t xn, prec_t err) {
	struct Polynomial_t * fp = poly_diff(f);
	while (poly_eval(f,xn) > err) {
		xn = xn - (poly_eval(f,xn)/poly_eval(fp,xn));
	}
	return xn;
}

struct Polynomial_t * roots_trinomial(struct Polynomial_t * self) {
	struct Polynomial_t * roots = poly_init_bare(3);
	for (int i=0; i < self->num_terms; ++i) {
		roots->terms[i]=self->terms[i];
	}
	if (self->terms[0] != 1 ) {
		poly_depress(roots);
	}
	prec_t shift_factor = 0;
	prec_t bp,cp;
	if (self->terms[1] !=0) {
		prec_t b = self->terms[1];
		shift_factor=b/3;
		roots = poly_recenter(roots,b/3);
	}
	bp = -roots->terms[2];
	cp = roots->terms[3];
	int num_zeros=0;
	prec_t start_locs[4] = {NAN,NAN,NAN,NAN},
	zeros[4] = {NAN,NAN,NAN,NAN};
	/* 4 because I want a maximum of 3, but the last element will be NULl. */

	// NOTE: see the documentation (online) for wtf is going on here. Got this from a textbook.
	// The number comments refer to which decision box that corresponds to
	// 1.
	if (bp > 0) {
		// 2.
		if (fabs(cp) == 2*pow(bp/3,1.5)) {
			zeros[0] = -sqrt(bp/3)*signbit(cp);
			zeros[1] = -sqrt(bp/3)*signbit(cp);
			start_locs[0] = (cp/2/bp + signbit(cp)*sqrt(bp));
			num_zeros = 2;
		}
		else {
			// 3.
			if (fabs(cp) < 2*pow(bp/3,1.5)) {
				if (cp==0) {
					start_locs[0]=(1./2/bp + sqrt(bp));
					start_locs[1]=(1./2/bp - sqrt(bp));
					start_locs[2]=(0);
				} else {
					start_locs[0]=(cp/2/bp + signbit(cp)*sqrt(bp));
					start_locs[1]=(cp/2/bp - signbit(cp)*sqrt(bp));
					start_locs[2]=(-cp/bp);
				}
			} else {
				// 4.
				if (pow(cp,2) > fabs(pow(bp,3))) {
					start_locs[0] = pow(cp,1./3);
				} else {
					start_locs[0]=(cp/2/bp + signbit(cp)*sqrt(bp));
				}
			}
		}
	} else {
		// 5.
		if (pow(cp,2) > fabs(pow(bp,3))) {
			start_locs[0] = pow(cp,1./3);
		} else {
			start_locs[0]=-cp/2/bp;
		}
	}
	int i=0;
	while (start_locs[i] != NAN) {
		zeros[num_zeros++]=newton(self,start_locs[i],1e-8)+shift_factor;
	}
}
