#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "include/pade.h"
#include "include/deriv.h"
#include "include/poly.h"
#include "include/sara.h"

void print_results(struct Solver_t * SOLV) {
	int n=SOLV->eqs->num->num_terms;
	prec_c_t output=SOLV->yy[n][0];
	/* printf("%f,%"PRNT_SPEC"%+"PRNT_SPEC"i\n",SOLV->curr_t,creal(output),cimag(output)); */
	printf("%f,%lf\n",SOLV->curr_t,creal(output));
}

prec_c_t L(prec_c_t f) {
	float a=1e-9,
			b=2.8e-9,
			c=800e-9,
			f0=2e4;
	prec_t res=(0.6366*a)*atan(-c*(f-f0))+b;
	prec_t offset=(0.6366*a)*atan(-c*(1e8-f0))+b;
	return res-offset;
}

int main() {
	float center=2e5;
	int M=1, N=2, x_step=10;
	prec_c_t * derivs = take_derivatives(&L,center,30);
	struct Polynomial_t * taylor = poly_init_bare(N+M+1);
	flip_arr(derivs, N+M+1);
	taylor->terms = derivs;
	struct Polynomial_t * out=poly_recenter(taylor,2e5);
	struct Pade_t * approx = pade_create_fit(out,M);
	approx->offset=L(1e8);
	struct Pade_t * sep=pade_separate(approx);
	sep->offset=L(1e8);
	pade_print(approx);
	pade_print_roots(sep);
	poly_free(taylor);
	poly_free(out);
	pade_free(approx);
	pade_free(sep);
	return 0;

	struct Polynomial_t * num = poly_init_bare(2);
	num->terms=malloc(2*sizeof(prec_c_t));
	num->terms[0]=0.9478;
	num->terms[1]=0.0522;
	struct Polynomial_t * denom = poly_init_bare(2);
	denom->terms=malloc(2*sizeof(prec_c_t));
	denom->terms[0] = -2.105;
	denom->terms[1] = -0.095;
	sep=pade_init_poly(num,denom);
	sep->offset=0;
	sep->vals=Roots;
	struct Solver_t * solv = solver_init(2, sep);
	solv->cb=&print_results;
	float time=0;
	printf("t,v\n");
	for (int i=0; i < 40; ++i) {
		step(solv,1,i);
	}
	/* prec_t res=0; */
	/* prec_t inpts[2]={1,2}; */
	/* printf("TIME\tV\n"); */
	/* for (int i=0; i < 40; ++i) { */
	/* 	inpts[0]+=1; */
	/* 	res=result(2,inpts); */
	/* 	printf("%e\t%f\n",inpts[0],res); */
	/* } */
	return 0;
}
