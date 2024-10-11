#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "include/pade.h"
#include "include/deriv.h"
#include "include/poly.h"
#include "include/sara.h"
#include "include/saber.h"

void print_results(struct Solver_t * SOLV, double output) {
	int n=SOLV->eqs->num->num_terms;
	/* printf("%f,%"PRNT_SPEC"%+"PRNT_SPEC"i\n",SOLV->curr_t,creal(output),cimag(output)); */
	printf("%e,%e\n",SOLV->curr_t,output);
}
void pa(prec_c_t *arr, int ne)
{
	for (int i=0; i < ne; ++i)
		printf("%d->%e\n",i,creal(arr[i]));
}
prec_c_t L(prec_t f) {
	float a=1e-9,
			b=2.8e-9,
			c=800e-9,
			f0=2e4;
	if (f==0)
		return 1/((0.6366*a)*atan(-c*(1e8-f0))+b);
	prec_t res=(0.6366*a)*atan(-c*(f-f0))+b;
	prec_t offset=(0.6366*a)*atan(-c*(1e8-f0))+b;
	return 1/(res-offset);
}
void call_saber_api(struct Pade_t * p) {
	/* double inputs[19]; */
	double * inputs=calloc(sizeof(double),19);
	int * nin=malloc(sizeof(int));
	int * unused=malloc(sizeof(int));
	double * outputs=malloc(sizeof(double));
	double * aundef=malloc(sizeof(double));

	nin[0]=19;
	unused[0]=-1;
	outputs[0]=0.;
	aundef[0]=0.;
	inputs[0]=1;
	inputs[1]=creal(p->offset);
	inputs[2]=cimag(p->offset);
	int i=0;
	for (; i<p->num->num_terms; ++i) {
		int elem_idx=2*i+3;
		inputs[elem_idx]=creal(p->num->terms[i]);
		inputs[elem_idx+1]=cimag(p->num->terms[i]);
	}
	for (; i<4; ++i) {
		int elem_idx=2*(i+1)+1;
		inputs[elem_idx]=0.;
		inputs[elem_idx+1]=0.;
	}
	i=0;
	for (; i<p->denom->num_terms; ++i) {
		int elem_idx=2*(i)+11;
		inputs[elem_idx]=creal(p->denom->terms[i]);
		inputs[elem_idx+1]=cimag(p->denom->terms[i]);
	}
	for (; i<4; ++i) {
		int elem_idx=2*(i+1)+11;
		inputs[elem_idx]=0.;
		inputs[elem_idx+1]=0.;
	}
	IND(inputs,nin,unused,unused,outputs,unused,unused,unused,aundef,unused);

	free(inputs);
	free(nin);
	/* free(unused); */
	/* free(outputs); */
	/* free(aundef); */
}
int main() {
	float center=1e6;
	int M=1, N=2, x_step=10;
	prec_c_t * derivs = take_derivatives(&L,center,40);
	pa(derivs,3);
	struct Polynomial_t * taylor = poly_init_bare(N+M+1);
	flip_arr(derivs, N+M+1);
	taylor->terms = derivs;
	struct Polynomial_t * out=poly_recenter(taylor,center);
	struct Pade_t * approx = pade_create_fit(out,M);
	pade_print(approx);
	approx->offset=L(0);
	struct Pade_t * sep=pade_separate(approx);
	sep->offset=L(0);
	pade_print(sep);
	struct Solver_t * solv = solver_init(2, sep);
	solv->cb=&print_results;
	double time=0;
	double delta_n=1e-9;
	/* printf("t,v\n"); */
	solv->eqs->offset=1e-9;
	for (int i=0; i < 20; ++i) {
		time+=delta_n;
		step(solv,1.,time);
	}
	poly_free(taylor);
	poly_free(out);
	pade_free(approx);
	pade_free(sep);
	return 0;
}
