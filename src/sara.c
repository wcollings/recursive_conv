/*
 * SARA - Semi-Analytical Recursive Algorithm
 * Implements Recursive Convolution
 *
*/
#include "../include/sara.h"
#include "../include/pade.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MIN(a,b) (a<b?a:b)

/*
 * `order` The order of the solver (1, 2, or 3). Least precise (but fastest) to most precise (but slowest)
 * `time` The stop time of the simulation
 * `delta_n` The time step to use for the simulation
 * `num_outputs` The number of output variables to save
 * `returns` a fully initialized Solver object
*/
struct solver_t * init_solver(int order,struct Pade_t * eq) {
	struct solver_t * self = malloc(sizeof(struct solver_t));
	printf("Creating solver obj!\n");
	self->cb=NULL;
	self->eqs=eq;
	self->order=order;
	self->curr_t=0;
	self->xx = calloc(sizeof(double),order);
	self->tt = calloc(sizeof(double),order);
	self->yy=malloc(sizeof(double*)*eq->M);
	for (int i=0; i < eq->M; ++i) {
		self->yy[i]=calloc(sizeof(double),order);
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

/*
 * Shift an array over by one, essentially dequeueing the last element
 * `arr`: the array to shift
 * `num_ele`: the size of the array
*/
void shift(double * arr,int num_ele) {
	for (int i=num_ele-1; i > 0; --i) {
		arr[i]=arr[i-1];
	}
}

/*
 * Take a step in the time domain, advancing the solver by one time step
 * `SOLV`: an instance of the solver
 * `inpt`: The state variable (likey voltage or current) at the next time step
*/
void step(struct solver_t * SOLV, double inpt, float curr_t) {
	// enqueue the newest input
	shift(SOLV->xx,SOLV->order);
	SOLV->xx[0]=inpt;
	shift(SOLV->tt,SOLV->order);
	SOLV->tt[0]=curr_t-SOLV->curr_t;
	double delta_n=SOLV->tt[0];
	double temp;
	for (int i=0; i < SOLV->order; ++i) {
		temp=0;
		float * Q = SOLV->q(SOLV->eqs->B[i],delta_n);
		temp=SOLV->yy[i][0]*Phi(SOLV->eqs->B[i],SOLV->tt[i]);
		for (int j=0; j <SOLV->order; ++j) {
			temp+=SOLV->eqs->A[i]*Q[j]*SOLV->xx[j];
		}
		shift(SOLV->yy[i],SOLV->order);
		SOLV->yy[i][0]=temp;
	}
	if (SOLV->cb != NULL)
		(*SOLV->cb)(SOLV);
	SOLV->curr_t=curr_t;
}
struct solver_t * gSOLV=NULL;
double result(int argc, double *argv) {
	if (gSOLV==NULL) {
		double A[]={0.5,-0.25};
		double B[]={0,0.020833,0,-0.0020833};
		struct Pade_t * p=Pade_init(A,B,2,4);
		gSOLV=init_solver(2, p);
	}
	float delta_n=argv[0]-gSOLV->curr_t;
	float temp;
	gSOLV->xx[0]=argv[1];
	for (int i=0; i < gSOLV->order; ++i) {
		temp=0;
		float * Q = gSOLV->q(gSOLV->eqs->B[i],delta_n);
		temp=gSOLV->yy[i][0]*Phi(gSOLV->eqs->B[i],delta_n);
		for (int j=0; j <gSOLV->order; ++j) {
			temp+=gSOLV->eqs->A[i]*Q[j]*gSOLV->xx[j];
		}
		shift(gSOLV->yy[i],gSOLV->order);
		gSOLV->yy[i][0]=temp;
	}
	gSOLV->curr_t=argv[0];
	return gSOLV->yy[0][0] + gSOLV->yy[1][0];
}
