#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "include/pade.h"
#include "include/ddiff.h"
#include "include/poly.h"
#include "include/sara.h"

void print_results(struct Solver_t * SOLV) {
	float output=SOLV->yy[0][0] + SOLV->yy[1][0];
	printf("%f,%f\n",SOLV->curr_t,output);
}

double L(float f) {
	float a=1e-9,
			b=2.8e-9,
			c=800e-9,
			f0=2e4;
	double res=(0.6366*a)*atan(-c*(f-f0))+b;
	return res;
}

int main() {
	float center=2e5;
	int M=3, N=4, x_step=10;
	double * base_vals = malloc(sizeof(double)*(M+N+1));
	double * xs = malloc(sizeof(double)*(M+N+1));
	for (int i=0; i < M+N+1; ++i) {
		xs[i] = center*(i+x_step);
		base_vals[i]=L(center+(i*x_step));
	}
	double * derivs = divided_difference(xs,base_vals,M+N+1);
	scale_to_taylor(derivs, N+M+1);
	struct Polynomial_t * taylor = poly_init_bare(N+M+1);
	flip_arr(taylor->terms, N+M+1);
	poly_recenter(taylor, center);
	struct Pade_t * approx = pade_create_fit(taylor,M);
	struct Solver_t * solv = solver_init(2, approx);
	solv->cb=&print_results;
	float time=0;
	printf("t,v\n");
	for (int i=0; i < 40; ++i) {
		step(solv,1,i);
	}
	return 0;
	/* double res=0; */
	/* double inpts[2]={1,2}; */
	/* printf("TIME\tV\n"); */
	/* for (int i=0; i < 40; ++i) { */
	/* 	inpts[0]+=1; */
	/* 	res=result(2,inpts); */
	/* 	printf("%e\t%f\n",inpts[0],res); */
	/* } */
	/* return 0; */
}
