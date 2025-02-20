#include "../include/poly.h"

void poly_iter(struct Polynomial_t * self) {
	prec_c_t temp;
	for (int i=0; i <self->num_terms; ++i) {
		temp=self->terms[i];
	}
}

//#define SUMUP(x) (int)((x*(x+1.))/2)
int SUMUP(int x) {
	int num = x*(x+1);
	return num/2;
}

struct Polynomial_t * poly_init_bare(int num_terms) {
	struct Polynomial_t * out = malloc(sizeof(struct Polynomial_t));
	out->tp=Vals;
	out->num_terms=num_terms;
	out->coeff=1;
	/* out->terms = calloc(num_terms,sizeof(prec_c_t)); */
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
	res->terms=malloc(sizeof(prec_c_t)*res->num_terms);
	res->coeff=1;
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
	res->coeff=1;
	res->terms=calloc(res->num_terms,sizeof(prec_c_t));
	for (int i=0; i < a->num_terms; ++i) {
		for (int j=0; j < b->num_terms; ++j ) {
			res->terms[i+j]=a->terms[i]*b->terms[i];
		}
	}
	return res;
}

prec_c_t poly_eval(struct Polynomial_t * a, prec_c_t val) {
	prec_c_t res=a->terms[0];
	for (int i=1; i < a->num_terms; ++i) {
		res=a->terms[i]+val*res;
	}
	return res*a->coeff;
}

void poly_print(struct Polynomial_t * p) {
	if (p->tp==Vals) {
		printf("(%+"PRNT_SPEC"%+"PRNT_SPEC"i)x^%d",creal(p->terms[0]),cimag(p->terms[0]),p->num_terms-1);
		for (int i=1; i < p->num_terms; ++i) {
			printf("(%+"PRNT_SPEC"%+"PRNT_SPEC"i)x^%d",creal(p->terms[i]),cimag(p->terms[i]),p->num_terms-i-1);
		}
	} else {
		for (int i=0; i < p->num_terms; ++i) {
			printf("(x%-"PRNT_SPEC"%+"PRNT_SPEC"i)",creal(p->terms[i]),cimag(p->terms[i]));
		}
	}
	printf("\n");
}

struct Polynomial_t * poly_sd_1term(struct Polynomial_t * num, prec_c_t z) {
	struct Polynomial_t * res=poly_init_bare(num->num_terms-1);
	res->terms = malloc((num->num_terms-1)*sizeof(prec_c_t));
	prec_c_t last=num->terms[0];
	res->terms[0]=last;
	int i;
	for (i=1; i < res->num_terms; ++i) {
		res->terms[i]=last*z+num->terms[i];
		last = res->terms[i];
	}
	/* poly_depress(res); */
	return res;
}

struct Polynomial_t * poly_get_roots(struct Polynomial_t * p) {
	struct Polynomial_t * roots;
	switch(p->num_terms-1) {
		case 1:
			roots=poly_init_bare(1);
			roots->terms[0] = p->terms[1]/p->terms[0];
			break;
		case 2:
			roots=roots_binomial(p);
			break;
		case 3:
			roots = roots_trinomial(p);
			break;
		default:
			fprintf(stderr,"Warning: cannot find the roots of %d-degree polynomial\n",p->num_terms-1);
	}
	return roots;
}

void poly_free(struct Polynomial_t * p) {
	free(p->terms);
	free(p);
}

void flip_arr(prec_c_t * arr, int n) {
	prec_c_t temp;
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
	int num_ele = SUMUP(row+1);
	int * tri = malloc(sizeof(int)*num_ele);
	int col=0;
	for (int i=0; i < num_ele; ++i) {
		int start=SUMUP(col);
		int end = SUMUP(col+1)-1;
		if (i == start || i == end) {
			tri[i]=1;
			if (i==end) col++;
		}
		else {
			int offset=i-SUMUP(col);
			int idx=SUMUP(col-1)+offset-1;
			tri[i]=tri[idx]+tri[idx+1];
		}
		/* int offset=i-SUMUP(col); */
	}
	for (int i = 1; i < row+1; ++i) {
		for (int j=1; j < floor((i-1.)/2+2); ++j) {
			tri[SUMUP(i)+2*j-1]*=-1;
		}
	}
	return tri;
}

struct Polynomial_t * poly_recenter(struct Polynomial_t * src, prec_c_t c) {
	int N = src->num_terms;
	int * tri = pascal(N, 1);
	struct Polynomial_t * dest = poly_init_bare(N);
	dest->terms = malloc(N*sizeof(prec_c_t));
	for (int i=0; i < N; ++i) {
		prec_c_t term = 0;
		for (int j = 0; j<=i; ++j) {
			int loc = SUMUP(N-j-1)+(i-j);
			prec_c_t bn=src->terms[j];
			term += tri[loc]*bn*pow(c,i-j); //NOLINT
		}
		dest->terms[i] = term;
	}
	free(tri);
	return dest;
}

void poly_depress(struct Polynomial_t * self) {
	prec_c_t max_term = self->terms[0];
	self->coeff=max_term;
	for (int i=0; i < self->num_terms; ++i) {
		self->terms[i]/=max_term;
	}
}

struct Polynomial_t * poly_diff(struct Polynomial_t * self) {
	struct Polynomial_t * out = poly_init_bare(self->num_terms-1);
	out->terms = calloc(self->num_terms-1,sizeof(prec_c_t));
	int term = self->num_terms-1;
	for (int i=0; i < term; ++i) {
		out->terms[i] = self->terms[i]*(term-i);
	}
	return out;
}

struct Polynomial_t * roots_binomial(struct Polynomial_t * self) {
	prec_c_t a,b,c;
	struct Polynomial_t * roots = poly_init_bare(2);
	roots->terms = malloc((self->num_terms-1)*sizeof(prec_c_t));
	roots->tp=Roots;
	/* poly_depress(self); */
	roots->coeff=self->coeff;
	a=self->terms[0];
	b=self->terms[1]/2/a;
	c=self->terms[2]/a;
	if (cabs(b) > 0) {
		if (cabs(cpow(b,2)) > cabs(c)) {
			roots->terms[0] = -(b+sqrt(pow(b,2)-c));
			roots->terms[1] = c/roots->terms[0];
		}
		else {
			roots->terms[0] = (-b+I*sqrt(c-pow(b,2)));
			roots->terms[1] = (-b-I*sqrt(c-pow(b,2)));
		}
	}
	else {
		if (cabs(cpow(b,2)) > cabs(c)) {
			roots->terms[0] = -b+sqrt(pow(b,2)-c);
			roots->terms[1] = c/roots->terms[0];

		}
		else {
			roots->terms[0] = (-b+I*sqrt(c-pow(b,2)));
			roots->terms[1] = (-b-I*sqrt(c-pow(b,2)));
		}
	}
	return roots;
}

prec_c_t newton(struct Polynomial_t * f, prec_c_t xn, prec_c_t err) {
	struct Polynomial_t * fp = poly_diff(f);
	while (cabs(poly_eval(f,xn)) > cabs(err)) {
		xn = xn - (poly_eval(f,xn)/poly_eval(fp,xn));
	}
	poly_free(fp);
	return xn;
}

prec_c_t ccbrt(prec_c_t z) {
	prec_c_t r = cabs(z);
	prec_c_t theta = carg(z);
	return cbrt(r)*cexp(I*theta/3);
}

struct Polynomial_t * roots_trinomial(struct Polynomial_t * self) {
	if (self->terms[0] != 1) {
		poly_depress(self);
	}
	struct Polynomial_t * inner = self;
	int inner_free = 0;
	if (self->terms[1] != 0) {
		inner = poly_recenter(self,self->terms[1]/3);
		inner_free = 1;
	}
	prec_c_t Q = inner->terms[2]/3;
	prec_c_t R = inner->terms[3]/2;
	prec_c_t D = pow(Q,3) + pow(R,2);
	prec_c_t S = ccbrt(R+sqrt(D));
	prec_c_t T = ccbrt(R-sqrt(D));
	prec_c_t coeff = sqrt(3)/2;
	struct Polynomial_t * ret = poly_init_bare(3);
	ret->tp=Roots;
	ret->terms[0] = S+T;
	ret->terms[1] = (S+T)/-2. + I*coeff*csqrt(S-T);
	ret->terms[2] = (S+T)/-2. - I*coeff*csqrt(S-T);
	if (inner_free) poly_free(inner);
	return ret;
}
