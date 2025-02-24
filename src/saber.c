#include <stdint.h>
#include <stdio.h>
#include <stdlib.h> //NOLINT
#include "../include/sara.h"
#include "../include/poly.h"
#include "../include/saber.h"
#include "../include/log.h"
#include "../include/ll.h"

#define JOB inp[0]
#define NAME_IDX inp[1]
#define vl inp[2]
#define t inp[3]

#define iL out[0]

#define SOLVER_NAME cgetstr(NAME_IDX)
#define SOLVER find_obj(SOLVER_NAME)


// Within the do_setup, the array is re-indexed to (x-2). This is because the first two elements don't need to 
// be considered within that function, so they're just stripped off.
// The following are array indices consistant with the setup portion
#define MODE (int)in[0]
#define Kstart 3
#define sstart 11

struct Solver_t * do_setup(double * in) {
	int order=4;
	prec_t * k_terms;
	k_terms=&in[Kstart];
	prec_t * s_terms=&in[sstart];
	for (int i=0; i < 4; ++i) {
		int idx=2*i;
		if (k_terms[idx]==0 && k_terms[idx+1]==0) {
			order=i;
			break;
		}
	}
	/* printf("Found %d terms\n",order); */
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
	return solver_init(2,pade,MODE);
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
	nout[0]=1;
	if (nout[1] < nout[0]) {
		return;
	}
	enum call_tp input=(int)JOB;
	switch (input) {
		case INIT: 
						add(SOLVER_NAME,do_setup(&inp[2]));
						solver_reset_time(SOLVER);
						nout[0]=1;
						iL = 1;
						break;
		default:
		case STEP:
						iL=step(SOLVER,vl,t);
						break;
		case ACCEPT: 
						iL=accept(SOLVER,vl,t);
						break;
		case END:
						solver_reset_time(SOLVER);
						break;
	}
	return;
}
