#include <stdint.h>
#include <stdio.h>
#include <stdlib.h> //NOLINT
#include "../include/sara.h"
#include "../include/poly.h"
#include "../include/saber.h"
#include "../include/log.h"
#include "../include/ll.h"

#define JOB inp[0]
#define SETUP 1
#define STEP 2
#define ACCEPT 3
#define iL out[0]
#define t inp[2]
#define vl inp[1]

#define MODE inp[2]
#define L 1
#define X 0
#define Kstart 4
#define sstart 12

struct Solver_t * do_setup(double * in) {
	int order=4;
	SOLV=malloc(sizeof(struct Solver_t));
	prec_t * k_terms;
	k_terms=&in[Kstart];
	prec_t * s_terms=&in[sstart];
	for (int i=0; i < 4; ++i)
	{
		int idx=2*i;
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
	pade->offset=in[1]+in[2]*I;
	return solver_init(2,pade,(int)in[0]);
	
	SOLV->head.order=order;
	SOLV->head.mode=(uint32_t)in[0];
	SOLV->eqs=pade;
	SOLV->curr_t=0;
	SOLV->curr_x=0;
	SOLV->num_calls=0;
	SOLV->num_steps=0;
	SOLV->tt=calloc(order,sizeof(prec_t));
	SOLV->xx=calloc(order,sizeof(prec_t));
	SOLV->yy=calloc(order,sizeof(prec_c_t));
	SOLV->qq=q2;
	/* switch (order) { */
	/* 	case 1: SOLV->qq=q1; */
	/* 			  break; */
	/* 	case 2: SOLV->qq=q2; */
	/* 			  break; */
	/* 	case 3: SOLV->qq=q3; */
	/* 			  break; */
	/* 	case 4: SOLV->qq=q4; */
	/* 			  break; */
	/* 	default: SOLV->qq=q2; */
	/* } */
	/* write_solver(SOLV,"solv_obj_save.bin"); */
	/* return SOLV; */
};

//SABER_FOREIGN_ROUTINE(IND) {
void IND(
        double* inp, 	//used
        int*   ninp,	//used
        int*   ifl,
        int*   nifl,
        double* out,	//used
        int*   nout,	//used
        int*   ofl,
        int*   nofl,
        double* aundef,	//used
        int*   ier
		  )
{
	int garbage;
	/* char * part_name = C_NAME(cgetstr)(inp[1]); */
	nout[0]=1;
	if (nout[1] < nout[0]) {
		log_msg("nout[1] < nout[0]");
		log_close();
		return;
	}
	switch ((int)JOB) {
		case SETUP: do_setup(&inp[1]);
						/* solvers[0]=*solv; */
						nout[0]=1;
						iL = 1;
						break;
		case STEP: iL=step(SOLV,vl,t);
					garbage=1;
		//printf("step");
						//iL=0;
					  break;
		case ACCEPT: iL=accept(SOLV,vl,t);
					garbage=1;
			//printf("accept");
						//iL=0;
					  break;
		default: iL=step(SOLV,vl,t);
			//printf("default - step");
						//iL=0;
					  break;
	}
	log_msg("Closing log");
	log_close();
	return;
}
