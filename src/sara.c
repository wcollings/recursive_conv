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
	printf("Creating solver obj!\n");
	self->cb=NULL;
	self->eqs=eq;
	self->order=order;
	self->curr_t=0;
	self->xx = calloc(sizeof(prec_t),order);
	self->tt = calloc(sizeof(prec_t),order);
	self->yy=malloc(sizeof(prec_t*)*eq->num->num_terms);
	for (int i=0; i < eq->num->num_terms; ++i) {
		self->yy[i]=calloc(sizeof(prec_t),order);
	}
	switch(order) {
		default:
		case 1: self->q=q1;
				 break;
		case 2: self->q=q2;
				  break;
		case 3: self->q=q3;
				  break;
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

float * q1(float i,float delta_n) {
	static float *q, _i, _n;
	if (_i != i || _n !=delta_n) {
		_i=i;
		_n=delta_n;
		q=malloc(sizeof(float));
		float zi=zeta(i,delta_n);
		q[0]=(delta_n/zi)*(1-Phi(i,delta_n));
	}
	return q;
}

float * q2(float i,float delta_n) {
	static float *q, _i, _n;
	if (_i != i || _n !=delta_n) {
		_i=i;
		_n=delta_n;
		q=malloc(sizeof(float)*2);
		float phi=Phi(i,delta_n);
		float zi=zeta(i,delta_n);
		if (zi==0) {
			q[0]=0;
			q[1]=0;
		}
		else {
			q[0]=(delta_n/pow(zi,2))*(-1+zi+phi);
			q[1]=(delta_n/pow(zi,2))*(1-(1+zi)*phi);
		}
	}
	return q;

}
float * q3(float i,float delta_n) {
	static float *q, _i, _n;
	if (_i != i || _n !=delta_n) {
		_i=i;
		_n=delta_n;
		float zi=zeta(i,delta_n);
		q=malloc(sizeof(float)*3);
		q[0]=(delta_n/(2*pow(zi,3)))*(2-(3*zi)+(2*pow(zi,2)) - (2-zi)*Phi(i,delta_n));
		q[1]=(delta_n/pow(zi,3))*(-2*(1-zi)+(2-pow(zi,2))*Phi(i,delta_n));
		q[2]=(delta_n/(2*pow(zi,3)))*(2-zi-(2+zi)*Phi(i,delta_n));
	}
	return q;
}

float * q4(float i,float delta_n) {
	static int _i,_n;
	static float *q;
	if (_i != i || _n !=delta_n) {
		_i=i;
		_n=delta_n;
		float zi=zeta(i,delta_n);
		q=malloc(sizeof(float)*4);
	}
	return q;
}

void shift(prec_t * arr,int num_ele) {
	for (int i=num_ele-1; i > 0; --i) {
		arr[i]=arr[i-1];
	}
}

void step(struct Solver_t * SOLV, prec_t inpt, float curr_t) {
	// enqueue the newest input
	shift(SOLV->xx,SOLV->order);
	SOLV->xx[0]=inpt;
	shift(SOLV->tt,SOLV->order);
	SOLV->tt[0]=curr_t-SOLV->curr_t;
	prec_t delta_n=SOLV->tt[0];
	prec_t temp;
	for (int i=0; i < SOLV->order; ++i) {
		temp=0;
		prec_t sigma=SOLV->eqs->denom->terms[i];
		prec_t a=SOLV->eqs->num->terms[i];
		float * Q = SOLV->q(sigma,delta_n);
		temp=SOLV->yy[i][0]*Phi(sigma,SOLV->tt[i]);
		for (int j=0; j <SOLV->order; ++j) {
			temp+=a*Q[j]*SOLV->xx[j];
		}
		shift(SOLV->yy[i],SOLV->order);
		SOLV->yy[i][0]=temp;
	}
	if (SOLV->cb != NULL)
		(*SOLV->cb)(SOLV);
	SOLV->curr_t=curr_t;
}

void step_resistance(struct Solver_t * SOLV, prec_t inpt, float curr_t) {
	shift(SOLV->xx,SOLV->order);
	SOLV->xx[0]=inpt;
	shift(SOLV->tt,SOLV->order);
	shift(SOLV->yy[0],SOLV->order);
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
	for (int i=0; i < gSOLV->order; ++i) {
		temp=0;
		prec_t sigma=gSOLV->eqs->denom->terms[i];
		prec_t a=gSOLV->eqs->num->terms[i];
		float * Q = gSOLV->q(sigma,delta_n);
		temp=gSOLV->yy[i][0]*Phi(sigma,delta_n);
		for (int j=0; j <gSOLV->order; ++j) {
			temp+=a*Q[j]*gSOLV->xx[j];
		}
		shift(gSOLV->yy[i],gSOLV->order);
		gSOLV->yy[i][0]=temp;
	}
	gSOLV->curr_t=argv[0];
	return gSOLV->yy[0][0] + gSOLV->yy[1][0];
}
