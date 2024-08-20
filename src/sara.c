/*
 * SARA - Semi-Analytical Recursive Algorithm
 * Implements Recursive Convolution
 *
*/
#include "../include/sara.h"
#include "../include/pade.h"
#include "../include/interpolate.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MIN(a,b) (a<b?a:b)

struct Solver_t * solver_init(int order,struct Pade_t * eq) {
	struct Solver_t * self = malloc(sizeof(struct Solver_t));
	/* printf("Creating solver obj!\n"); */
	self->cb=NULL;
	self->eqs=eq;
	self->order=order;
	self->curr_t=0;
	self->curr_x=0;
	self->xx = calloc(sizeof(prec_t),order);
	self->tt = calloc(sizeof(prec_t),order);
	self->yy = calloc(sizeof(prec_t),order+1);
	switch (order) {
		case 1: self->qq=q1;
				  break;
		case 2: self->qq=q2;
				  break;
		case 3: self->qq=q3;
				  break;
		case 4: self->qq=q4;
				  break;
		default: self->qq=q2;
	}
	return self;
}

/*
 * Calculates $\zeta_{i,n}$
 * `i`: s_i
 * `n`: Delta_n
 */
float zeta(float i, float n) {
	return -(i*n);
}

/*
 * Calculates $\Phi_{i,n}$
 * `i`: s_i
 * `n`: Delta_n
 */
float Phi(float i, float n) {
	static float res,_i, _n;
	if (i!=_i || n!=_n) {
		_i=i;
		_n=n;
		res= exp(i*n);
	}
	return res;
}

/* Assumes n=1, won't even check otherwise.
 * This is only there so that the function signature lines up with the rest.
*/
prec_c_t q1(float sigma_i,float delta_n,int n) {
	prec_c_t q;
	prec_c_t zi=zeta(sigma_i,delta_n);
	q=(delta_n/zi)*(1-Phi(sigma_i,delta_n));
	return q;
}

prec_c_t q2(float sigma_i,float delta_n, int n) {
	prec_c_t q;
	prec_c_t phi=Phi(sigma_i,delta_n);
	prec_c_t zi=zeta(sigma_i,delta_n);
	if (delta_n==0) {
		q=0;
	}
	else switch(n) {
		case 0: q=(delta_n/cpow(zi,2))*(-1+zi+phi);
				  break;
		case 1: q=(delta_n/cpow(zi,2))*(1-(1+zi)*phi);
				  break;
	}
	return q;

}
prec_c_t q3(float sigma_i,float delta_n,int n) {
	prec_c_t q;
	prec_c_t zi=zeta(sigma_i,delta_n);
	prec_c_t phi = Phi(sigma_i,delta_n);
	if (zi==0) {
		q=0;
	} else switch(n) {
		case 0: q=(delta_n/(2*pow(zi,3)))*(2-(3*zi)+(2*pow(zi,2)) - (2-zi)*phi);
				  break;
		case 1: q=(delta_n/pow(zi,3))*(-2*(1-zi)+(2-pow(zi,2))*phi);
				  break;
		case 2: q=(delta_n/(2*pow(zi,3)))*(2-zi-(2+zi)*phi);
	}
	return q;
}

prec_c_t q4(float i,float delta_n, int n) {
	float q;
	float zi=zeta(i,delta_n);
	return q;
}

void shift(prec_t * arr,int num_ele) {
	for (int i=num_ele-1; i > 0; --i) {
		arr[i]=arr[i-1];
	}
}
void shift_c(prec_c_t * arr,int num_ele) {
	for (int i=num_ele-1; i > 0; --i) {
		arr[i]=arr[i-1];
	}
}

void step(struct Solver_t * SOLV, prec_t inpt, float curr_t) {
	prec_t delta_n = (curr_t-SOLV->curr_t);
	prec_t integ = delta_n * (inpt+SOLV->curr_x)/2;
	// enqueue the newest input
	shift(SOLV->xx,SOLV->order);
	SOLV->xx[0]=integ;
	shift(SOLV->tt,SOLV->order);
	SOLV->tt[0]=curr_t-SOLV->curr_t;
	SOLV->curr_t=curr_t;

	prec_c_t temp, final=SOLV->eqs->offset*inpt;
	for (int i=0; i < SOLV->order; ++i) {
		// calc y_i[n]
		prec_t delta_n=SOLV->tt[i];
		temp=0;
		prec_t sigma=SOLV->eqs->denom->terms[i];
		prec_t a=SOLV->eqs->num->terms[i];
		temp=SOLV->yy[i]*Phi(sigma,SOLV->tt[i]);
		for (int j=0; j <SOLV->order; ++j) {
			// calculate the q terms
			prec_c_t q=SOLV->qq(sigma,delta_n,j);
			temp+=a*q*SOLV->xx[j];
		}
		SOLV->yy[i]=temp;
		final += temp;
	}
	int n=SOLV->eqs->num->num_terms;
	SOLV->yy[n] = final;
	if (SOLV->cb != NULL) {
		(*SOLV->cb)(SOLV);
	}
}

/*
 * These values are almost certainly wrong, and the eqation I'm describing is probably wrong as well.
 * Reference the wiki for more accurate info on the actual equation. But this was my first attempt
 * anyway.
*/
void step_resistance(struct Solver_t * SOLV, prec_t inpt, float curr_t) {
	shift(SOLV->xx,SOLV->order);
	SOLV->xx[0]=inpt;
	shift(SOLV->tt,SOLV->order);
	shift_c(SOLV->yy[0],SOLV->order);
	SOLV->tt[0]=curr_t-SOLV->curr_t;
	prec_t b = 2.04e-9;
	prec_t c = -0.2739;
	prec_t d = -0.03652;
	prec_t c0 = 0.2/exp(b*c);
	prec_t interp_t = curr_t - b;
	prec_t out1 = c0*interpolate(SOLV->tt[1],
										  SOLV->xx[1],
										  SOLV->tt[0],
										  SOLV->xx[0],
										  interp_t);
	prec_t out2 = d*SOLV->xx[0];
	SOLV->yy[0][0] = out1+out2;
}
struct Solver_t * gSOLV=NULL;
prec_t result(int argc, prec_t *argv) {
	if (gSOLV==NULL) {
		prec_c_t A[]={0.5,-0.25};
		prec_c_t B[]={0,0.020833,0,-0.0020833};
		struct Pade_t * p=pade_init(A,B,2,4);
		gSOLV=solver_init(2, p);
	}
	float delta_n=argv[0]-gSOLV->curr_t;
	float temp;
	gSOLV->xx[0]=argv[1];
	// TODO: update this correctly
	// Right now I did just enough so that it would compile, because I'm testing what's above and don't need this yet
	for (int i=0; i < gSOLV->order; ++i) {
		temp=0;
		prec_t sigma=gSOLV->eqs->denom->terms[i];
		prec_t a=gSOLV->eqs->num->terms[i];
		float Q = gSOLV->qq(sigma,delta_n,0);
		temp=gSOLV->yy[i][0]*Phi(sigma,delta_n);
		for (int j=0; j <gSOLV->order; ++j) {
			temp+=a*Q*gSOLV->xx[j];
		}
		shift_c(gSOLV->yy[i],gSOLV->order);
		gSOLV->yy[i][0]=temp;
	}
	gSOLV->curr_t=argv[0];
	return gSOLV->yy[0][0] + gSOLV->yy[1][0];
}
void result_ac(int argc, prec_t *argv) {
	if (gSOLV==NULL) {
		prec_c_t A[]={0.5,-0.25};
		prec_c_t B[]={0,0.020833,0,-0.0020833};
		struct Pade_t * p=pade_init(A,B,2,4);
		gSOLV=solver_init(2, p);
	}
}
