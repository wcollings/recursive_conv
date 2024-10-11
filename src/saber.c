#include <stdint.h>
#include <stdio.h>
#include <stdlib.h> //NOLINT
#include "../include/sara.h"
#include "../include/poly.h"
#include "../include/saber.h"
#include "../include/log.h"

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
#define Kstart 3
#define sstart 11

struct Solver_t * do_setup(double * in) {
	int order=4;
	struct Solver_t * SOLV=malloc(sizeof(struct Solver_t));
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
	SOLV->head.order=order;
	SOLV->head.mode=(int)in[0];
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
	/* write_solver(SOLV,"solv_obj_save.bin"); */
	return SOLV;
};

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
{	char printThis[100];
	char isUndef;
	log_init("log.txt");
	nout[0]=1;
	if (nout[1] < nout[0]) {
		log_msg("nout[1] < nout[0]");
		log_close();
		return;
	}
	if (ninp[0] != 20) {
		snprintf(printThis,100, "ninp[0] != 20, instead = %d", ninp[0]);
		log_msg(printThis);
	}else {
		log_msg("ninp[0] == 20");
	}

	
	log_msg("Checking for undef");
	for(int i=0; i<ninp[0]; i++)
	{	
		isUndef = inp[i] == aundef[0]? 'T':'F';
		snprintf(printThis,100, "inp[%d] undefined: %c", i,isUndef);
		log_msg(printThis);
	}
	
	struct Solver_t * solv=&solvers[0];
	snprintf(printThis,100, "Called with arg %d", (int)JOB);
	log_msg(printThis);
	switch ((int)JOB) {
		case SETUP: solv=do_setup(&inp[1]);
		log_msg("do_setup");
						solvers[0]=*solv;
						nout[0]=1;
						iL = 1;
						break;
		case STEP: //iL=step(solv,vl,t);
		log_msg("step");
						iL=0;
					  break;
		case ACCEPT: //iL=accept(solv,vl,t);
		log_msg("accept");
						iL=0;
					  break;
		default: //iL=step(solv,vl,t);
		log_msg("default - step");
						iL=0;
					  break;
	}
	log_msg("Closing log");
	log_close();
	return;
}
