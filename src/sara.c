/*
 * SARA - Semi-Analytical Recursive Algorithm
 * Implements Recursive Convolution
 *
*/
#include <stdlib.h>
#include <math.h>
#include "../include/sara.h"
#include "../include/pade.h"
#define MIN(a,b) (a<b?a:b)

struct Solver_t * solver_init(int order,struct Pade_t * eq,int mode) {
	struct Solver_t * self = malloc(sizeof(struct Solver_t));
	self->num_calls=0;
	self->eqs=eq;
	self->head.order=order;
	self->head.mode=mode;
	self->curr_t=0;
	self->curr_x=0;
	self->xx = calloc(sizeof(prec_t),order);
	self->tt = calloc(sizeof(prec_t),order);
	self->yy=malloc(sizeof(prec_t*)*order);
	switch (order) {
		case 1: self->qq=q1;
				  break;
		case 2: self->qq=q2;
				  break;
		case 3: self->qq=q3;
				  break;
		default: self->qq=q2;
	}
	return self;
}
void solver_free(struct Solver_t * self) {
	pade_free(self->eqs);
	free(self->tt);
	free(self->xx);
	free(self->yy);
	free(self);
}

/*
 * If the simulation resets, or restarts, but the DLL isn't thrown out,
 * we need to essentially detect that and reset ourselves
*/
void solver_reset_time(struct Solver_t * self) {
	for (int i=0; i < self->head.order; ++i) {
		self->tt[i]=0;
		self->xx[i]=0;
		self->yy[i]=0;
	}
}

/*
 * Calculates $\zeta_{i,n}$
 * `i`: s_i
 * `n`: Delta_n
 */
prec_c_t zeta(prec_c_t i, prec_c_t n) {
	return -(i*n);
}

/*
 * Calculates $\Phi_{i,n}$
 * `i`: s_i
 * `n`: Delta_n
 */
prec_c_t Phi(prec_c_t i, prec_c_t n) {
	return cexp(i*n);
}

/* 
 * Assumes n=1, won't even check otherwise.
 * This is only there so that the function signature lines up with the rest.
*/
prec_c_t q1(prec_c_t sigma_i,prec_c_t delta_n,int n) {
	prec_c_t q;
	prec_c_t zi=zeta(sigma_i,delta_n);
	q=delta_n*(1-zi/2.);
	return q;
}

prec_c_t q2(prec_c_t sigma_i,prec_c_t delta_n, int n) {
	prec_c_t zi=zeta(sigma_i,delta_n);
	if (delta_n==0) {
		return 0;
	} else switch(n) {
		default:
		case 0: return delta_n/2;
		case 1: return (delta_n/2)*(1-zi);
	}
}
prec_c_t q3(prec_c_t sigma_i,prec_c_t delta_n,int n) {
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

prec_t step(struct Solver_t * SOLV, const prec_t inpt, const prec_t curr_t) {
	SOLV->num_steps++;
	if (curr_t < SOLV->curr_t) {
		solver_reset_time(SOLV);
	}
	prec_t delta_n=curr_t-SOLV->curr_t;
	prec_t new_x;
	if (SOLV->head.mode==INDUCTANCE) {
		prec_t last_x=SOLV->curr_x;
		prec_t integ = delta_n*((last_x+inpt)/2);
		prec_t new_x = SOLV->xx[0]+integ;
	} else {
		new_x = inpt;
	}

	prec_c_t temp = 0;
	prec_c_t final=SOLV->eqs->offset*new_x;

	prec_c_t outputs[4]={0,0,0,0};
	// can't save the values to SOLV yet, so have to do this jank instead
	// loop for j=0
	for (int i=0; i < SOLV->head.order; ++i) {
		prec_c_t sigma_i=SOLV->eqs->denom->terms[i];
		prec_c_t Ki=SOLV->eqs->num->terms[i];
		temp=SOLV->yy[i]*Phi(sigma_i,delta_n);
		prec_c_t q=SOLV->qq(sigma_i,delta_n,0);
		temp+=Ki*q*new_x;
		outputs[i]=temp;
	}
	// loop for j=1..n
	for (int i=0; i < SOLV->head.order; ++i) {
		prec_c_t sigma_i=SOLV->eqs->denom->terms[i];
		prec_c_t Ki=SOLV->eqs->num->terms[i];
		temp=0;
		for (int j=1; j <SOLV->head.order; ++j) {
			prec_t delta_n=SOLV->tt[i];
			prec_c_t q=SOLV->qq(sigma_i,delta_n,j);
			temp+=Ki*q*SOLV->xx[j-1];
		}
		outputs[i]+=temp;
		final += outputs[i];
	}
	return creal(final);
}

prec_t accept(struct Solver_t * SOLV, const prec_t inpt, const prec_t curr_t) {
	SOLV->num_calls++;
	shift(SOLV->tt,SOLV->head.order);
	SOLV->tt[0]=curr_t-SOLV->curr_t;
	SOLV->curr_t=curr_t;
	prec_t new_x;
	if (SOLV->head.mode==INDUCTANCE) {
		prec_t last_x=SOLV->curr_x;
		prec_t integ = SOLV->tt[0]*((last_x+inpt)/2);
		prec_t new_x = SOLV->xx[0]+integ;
	} else {
		new_x = inpt;
	}
	shift(SOLV->xx,SOLV->head.order);
	SOLV->xx[0]=new_x;
	SOLV->curr_x=inpt;

	prec_c_t temp;
	prec_c_t final=SOLV->eqs->offset*new_x;
	for (int i=0; i < SOLV->head.order; ++i) {
		prec_t delta_n=SOLV->tt[i];
		temp=0;
		prec_c_t sigma_i=SOLV->eqs->denom->terms[i];
		prec_c_t Ki=SOLV->eqs->num->terms[i];
		temp=SOLV->yy[i]*Phi(sigma_i,delta_n);
		for (int j=0; j <SOLV->head.order; ++j) {
			// calculate the q terms
			prec_c_t q=SOLV->qq(sigma_i,delta_n,j);
			temp+=Ki*q*SOLV->xx[j];
		}
		SOLV->yy[i]=temp;
		final += temp;
	}
	return creal(final);
}
