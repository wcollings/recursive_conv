#include <stdint.h>
#include <stdio.h>
#include <stdlib.h> //NOLINT
#include "../include/sara.h"
#include "../include/poly.h"
#include "../include/saber.h"
/* #include <saberApi.h> */
#define JOB inp[0]
#define SETUP 1
#define STEP 2
#define ACCEPT 3
#define iL out[0]

#define Kstart 2
#define sstart 10

struct init_array {
	prec_c_t k0;
	prec_c_t K1;
	prec_c_t K2;
	prec_c_t K3;
	prec_c_t K4;
	prec_c_t s1;
	prec_c_t s2;
	prec_c_t s3;
	prec_c_t s4;
};
struct work_array {
	double t;
	double vl;
};

struct inpt {
	int64_t sel;
	union input_wrapper{
		struct init_array init;
		struct work_array work;
	} wrap;
};

void pa_saber(double * arr,int ne) {
	for (int i=0; i<ne; ++i) {
		printf("%d->%le\n",i,arr[i]);
	}
}
double do_setup(double * in) {
	int order=4;
	struct Solver_t * SOLV=malloc(sizeof(struct Solver_t));
	prec_t * k_terms;
	k_terms=&in[Kstart];
	prec_t * s_terms=&in[sstart];
	for (int i=0; i < 4; ++i)
	{
		int idx=2*i;
		printf("K_%d=(%1.2le+i%1.2le)\n",i,k_terms[idx],k_terms[idx+1]);
		if (k_terms[idx]==0 && k_terms[idx+1]==0) {
			order=i;
			break;
		}
	}
	printf("Found %d terms\n",order);
	struct Polynomial_t * num=poly_init_bare(order);
	struct Polynomial_t * denom=poly_init_bare(order);
	num->terms=malloc(sizeof(prec_c_t)*order);
	denom->terms=malloc(sizeof(prec_c_t)*order);
	for (int i=0; i < order; ++i) {
		int idx=2*i;
		num->terms[i]=k_terms[idx]+k_terms[idx+1]*I;
		denom->terms[i]=s_terms[idx]+s_terms[idx+1]*I;
	}
	struct Pade_t * pade=pade_init_poly(num,denom);
	pade->offset=in[0]+in[1]*I;
	SOLV->order=order;
	SOLV->eqs=pade;
	SOLV->curr_t=0;
	SOLV->curr_x=0;
	SOLV->tt=malloc(sizeof(prec_t)*order);
	SOLV->xx=malloc(sizeof(prec_t)*order);
	SOLV->yy=malloc(sizeof(prec_c_t)*order);
	switch (order) {
		case 1: SOLV->qq=q1;
				  break;
		case 2: SOLV->qq=q2;
				  break;
		case 3: SOLV->qq=q3;
				  break;
		case 4: SOLV->qq=q4;
				  break;
		default: SOLV->qq=q2;
	}
	write_solver(SOLV,"solv_obj_save.bin");
	return 0;
};

void IND(
        double* inp, 
        int*   ninp,
        int*   ifl,
        int*   nifl,
        double* out,
        int*   nout,
        int*   ofl,
        int*   nofl,
        double* aundef,
        int*   ier)
{
	union{
		double * as_arr;
		struct inpt as_struct;
	} in_arr;
	printf("Got %d inputs\n",ninp[0]);
	printf("Command = %d\n",(int)JOB);
	in_arr.as_arr=inp;
	pa_saber(inp,ninp[0]);
	in_arr.as_struct.sel=(int)in_arr.as_arr[0];
	prec_t t=in_arr.as_struct.wrap.work.t;
	prec_t vl=in_arr.as_struct.wrap.work.vl;
	nout[0]=1;
	switch ((int)JOB) {
		case STEP: iL=do_step(t,vl);
		case ACCEPT: iL=do_accept(t,vl);
		default: iL=do_step(t,vl);
		case SETUP: iL=do_setup(in_arr.as_struct);
	}
	return;
}
