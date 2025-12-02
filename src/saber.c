#include <stdint.h>
#include <stdio.h>
#include <stdlib.h> //NOLINT
#include "../include/sara.h"
#include "../include/poly.h"
#include "../include/saber.h"

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

struct Solver_t * SOLV;

/*
 * The input array has the following structure: (square brackets indicate array index)
 * [0] -> Selector for "Pure Inductance" vs "Admittance" mode
 * [1,2] -> real and imaginary components of K0
 * [[3,4]],...[[9,10]] -> real and imaginary components of K1..K4. If there are less than 4 terms to the model, the rest are set to 0
 * [[11,12]],...[17,18]] -> real and imaginary components of sigma1..sigma4. If there are less than 4 terms to the model, the rest are set to 0
 *
 * Setup detects the first (0+0i) element in the Ki sub-vector, and creates the objects it needs with that info.
 * If there is no (0+0i) element, its assumed that all 4 terms are needed.
 * 
 * `num` refers to the Ki vector, and `denom` to the `sigma_i` vector,
 * since they come from (K/(s-sigma))
 */
void do_setup(double * in) {
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
	SOLV=solver_init(order,pade,(uint32_t)in[0]);
};

/* Need all these extra fields to interface correctly with their API
 * `inp` The array of inputs. Note that the length and contents of this will vary with the required action.
 *       The first element of this states what the required action is: 1 for init, 2 for step, 3 for accept.
 * `ninp` The number of inputs given, just in case you need that. Also used for passing strings somehow?
 * `ifl`: Unused
 * `nifl`: Unused, the number of elements in ifl
 * `out`: the output array
 * `nout`: the number of elements in the output array. It's already been decided for you, so if you need a different
 *         number of inputs, you have to change this number, do nothing else, and return. Your function will then get
 *         re-called with the correct number of outputs correctly allocated.
 *`ofl`: unused
 *`nofl`: the number elements in ofl (therefore also unused)
 *`aundef`: the value used to indicate "undefined". You will be given this, so just copy it wherever you need it
 *`ier`: unused
*/
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
        int*   ier
		  )
{
	nout[0]=1;
	if (nout[1] < nout[0]) {
		return;
	}
	switch ((int)JOB) {
		case SETUP: do_setup(&inp[1]);
						nout[0]=1;
						iL = 1;
						break;
		case ACCEPT: iL=accept(SOLV,vl,t);
					  break;
		case STEP:
		default: iL=step(SOLV,vl,t);
					  break;
	}
	return;
}
