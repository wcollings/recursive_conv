#include <stdint.h>
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

struct init_array {
	double k0;
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

double do_setup(struct inpt inp) {
	int order=4;
	struct Solver_t * SOLV=malloc(sizeof(struct Solver_t));
	prec_c_t * k_terms;
	k_terms=&inp.wrap.init.K1;
	prec_c_t * s_terms;
	s_terms=&inp.wrap.init.s1;
	for (int i=0; i < 4; ++i)
	{
		if (k_terms[i]==0) {
			order=i;
			break;
		}
	}
	struct Polynomial_t * num=poly_init_bare(order);
	struct Polynomial_t * denom=poly_init_bare(order);
	num->terms=malloc(sizeof(prec_c_t)*order);
	denom->terms=malloc(sizeof(prec_c_t)*order);
	for (int i=0; i < order; ++i) {
		num->terms[i]=k_terms[i];
		denom->terms[i]=s_terms[i];
	}
	struct Pade_t * pade=pade_init_poly(num,denom);
	pade->offset=inp.wrap.init.k0;
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

ADDAPI void ADDCALL IND(
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
	in_arr.as_arr=inp;
	in_arr.as_struct.sel=(int)in_arr.as_arr[0];
	prec_t t=in_arr.as_struct.wrap.work.t;
	prec_t vl=in_arr.as_struct.wrap.work.vl;
	switch ((int)JOB) {
		case STEP: iL=do_step(t,vl);
					  break;
		case ACCEPT: iL=do_accept(t,vl);
						 break;
		default: iL= do_step(t,vl);
					break;
		case SETUP: iL = do_setup(in_arr.as_struct);
	}
}
