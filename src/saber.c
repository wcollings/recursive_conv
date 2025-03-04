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

void log_params(double in[2]) {
	char * str = malloc(80);
	snprintf(str,80,"v=%e,t=%e\n",in[0],in[1]);
	log_msg(str);
}

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

/*
 * The start point for the whole DLL, and thus the algorithm.
 *
 * First argument in `inp` is always the control flag (see call_tp enum). past that, the arguments required are listed below
 *
 * `INIT`: gets called when the project is first opened
 *
 *  arguments: [ name, MODE, k0, k0i, k1, k1i, k2, k2i, k3, k3i, k4, k4i, sigma1, sigma1i, sigma2, sigma2i, sigma3, sigma3i, sigma4, sigma4i ]
 *  `STEP`: gets called during the time step decision iteration loop.
 *
 *  arguments: [ name, v, t ]
 *  `ACCEPT`: gets called once the time step decision iteration loop finishes, before the next time step iteration loop begins.
 *
 *  arguments: [ name, v, t ]
 *  `START`: called before the start of a transient simulation. Opens the log file.
 *  arguments: [ name ]
 *  `END`: called before the start of a transient simulation. Closes the log file and resets the solver objects.
 *  arguments: [ name ]
 *
 *  The full arguments are as follows:
 *  `inp`: as described above
 *  `ninp`: the size of the given `inp` array.
 *  `ifl`: an array of flags given
 *  `nifl`: the size of the `ifl` array
 *  `out`: the array to hold any and all outputs
 *  `nout`: the number of pre-allocated elements in the output array. If this is less than you need, set it to a larger number and return immediately. The func will be re-called with the right number and only then can you procede.
 *
 *  `olf`: an array of flags to be returned.
 *  `nolf`: the pre-allocated size of the `ofl` array.
 *  `aundef`: the used value that means "undefined" in the simulator. If you need to specify "undefined", just copy the value of this into wherever you need it.
 *  `ier`: if an error has occured.
 */
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
						nout[0]=1;
						iL = 1;
						break;
		case START:
						log_init("sara.log");
						log_msg("Starting transient simulation");
						solver_reset_time(SOLVER);
						break;
		default:
		case STEP:
						iL=step(SOLVER,vl,t);
						break;
		case ACCEPT: 
						iL=accept(SOLVER,vl,t);
						break;
		case END:
						log_msg("Ending transient simulation");
						log_close();
						break;
	}
	return;
}
