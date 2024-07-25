#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "include/pade.h"
#include "include/deriv.h"
#include "include/poly.h"
#include "include/sara.h"

void print_results(struct Solver_t * SOLV) {
	prec_t output=SOLV->yy[0][0] + SOLV->yy[1][0];
	printf("%f,%f\n",SOLV->curr_t,output);
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
	free(taylor->terms);
	flip_arr(derivs, N+M+1);
	taylor->terms = derivs;
	struct Polynomial_t * out=poly_recenter(taylor,2e5);
	struct Pade_t * approx = pade_create_fit(out,M);
	pade_separate(approx);
	poly_free(taylor);
	poly_free(out);
	pade_free(approx);
	/* struct Solver_t * solv = solver_init(2, approx); */
	/* solv->cb=&print_results; */
	return 0;
	/* float time=0; */
	/* printf("t,v\n"); */
	/* for (int i=0; i < 40; ++i) { */
	/* 	step(solv,1,i); */
	/* } */
	/* prec_t res=0; */
	/* prec_t inpts[2]={1,2}; */
	/* printf("TIME\tV\n"); */
	/* for (int i=0; i < 40; ++i) { */
	/* 	inpts[0]+=1; */
	/* 	res=result(2,inpts); */
	/* 	printf("%e\t%f\n",inpts[0],res); */
	/* } */
	/* return 0; */
}
