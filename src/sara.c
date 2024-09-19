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
	self->yy=malloc(sizeof(prec_t*)*order);
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
float zeta(prec_c_t i, prec_c_t n) {
	return -(i*n);
}

/*
 * Calculates $\Phi_{i,n}$
 * `i`: s_i
 * `n`: Delta_n
 */
float Phi(prec_c_t i, prec_c_t n) {
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
prec_c_t q1(prec_c_t sigma_i,prec_c_t delta_n,int n) {
	prec_c_t q;
	prec_c_t zi=zeta(sigma_i,delta_n);
	q=delta_n*(1-zi/2.);
	return q;
}

prec_c_t q2(prec_c_t sigma_i,prec_c_t delta_n, int n) {
	prec_c_t q;
	prec_c_t phi=Phi(sigma_i,delta_n);
	prec_c_t zi=zeta(sigma_i,delta_n);
	if (delta_n==0) {
		q=0;
	}
	else switch(n) {
		case 0: q=delta_n/2;
		/* case 0: q=(delta_n/cpow(zi,2))*(-1+zi+phi); */
				  break;
		case 1: q=(delta_n/2)*(1-zi);
		/* case 1: q=(delta_n/cpow(zi,2))*(1-(1+zi)*phi); */
				  break;
	}
	return q;

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

prec_c_t q4(prec_c_t i,prec_c_t delta_n, int n) {
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

prec_t step(struct Solver_t * SOLV, prec_t inpt, float curr_t) {
	shift(SOLV->tt,SOLV->order);
	SOLV->tt[0]=curr_t-SOLV->curr_t;
	SOLV->curr_t=curr_t;
	prec_t last_x=SOLV->curr_x;
	prec_t integ = SOLV->tt[0]*((last_x+inpt)/2);
	prec_t new_x = SOLV->xx[0]+integ;
	shift(SOLV->xx,SOLV->order);
	SOLV->xx[0]=new_x;
	SOLV->curr_x=inpt;

	prec_c_t temp;
	prec_c_t final=SOLV->eqs->offset*new_x;
	for (int i=0; i < SOLV->order; ++i) {
		prec_t delta_n=SOLV->tt[i];
		temp=0;
		prec_c_t sigma=SOLV->eqs->denom->terms[i];
		prec_c_t a=SOLV->eqs->num->terms[i];
		temp=SOLV->yy[i]*Phi(sigma,SOLV->tt[i]);
		for (int j=0; j <SOLV->order; ++j) {
			// calculate the q terms
			prec_c_t q=SOLV->qq(sigma,delta_n,j);
			temp+=a*q*SOLV->xx[j];
		}
		SOLV->yy[i]=temp;
		final += temp;
	}
	/* int n=SOLV->eqs->num->num_terms; */
	/* SOLV->yy[n] = final; */
	if (SOLV->cb != NULL) {
		(*SOLV->cb)(SOLV);
	}
	return creal(final);
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
	SOLV->yy[0] = out1+out2;
}
prec_t do_step(prec_t inpt, float curr_t) {
	struct Solver_t * solv = read_solver("solv_obj_save.bin");
	prec_t res=step(solv,inpt,curr_t);
	write_solver(solv,"solv_obj_save.bin");
	return res;
}
prec_t do_accept(prec_t inpt, float curr_t) {
	struct Solver_t * solv = read_solver("solv_obj_save.bin");
	prec_t res=step(solv,inpt,curr_t);
	write_solver(solv,"solv_obj_save.bin");
	return res;
}

void write_solver(struct Solver_t * s,char* fname) {
	FILE * fp = fopen(fname,"wb");
	fwrite(&s->order,sizeof(int),1,fp); 
	int num_terms=s->eqs->num->num_terms;
	fwrite(&num_terms,sizeof(int),1,fp); 
	fwrite(&s->eqs->offset,sizeof(prec_c_t),1,fp);
	fwrite(s->eqs->num->terms,sizeof(prec_c_t),num_terms,fp);
	fwrite(s->eqs->denom->terms,sizeof(prec_c_t),num_terms,fp);
	fwrite(&s->curr_t,sizeof(prec_t),1,fp);
	fwrite(&s->curr_x,sizeof(prec_t),1,fp);
	fwrite(s->tt,sizeof(prec_t),num_terms,fp);
	fwrite(s->xx,sizeof(prec_t),num_terms,fp);
	fwrite(s->yy,sizeof(prec_c_t),num_terms,fp);
	fclose(fp);
}
struct Solver_t * read_solver(char * fname) {
	FILE *fp = fopen(fname,"rb");
	struct Solver_t * self=malloc(sizeof(struct Solver_t));
	fread(&self->order,sizeof(int),1,fp);
	prec_c_t K0;
	fread(&K0,sizeof(prec_c_t),1,fp);
	// Read in the equation:
	// Ki terms
	struct Polynomial_t * num = poly_init_bare(self->order);
	num->terms = malloc(sizeof(prec_c_t)*self->order);
	fread(num->terms,sizeof(prec_c_t),self->order,fp);
	//sigma_i terms
	struct Polynomial_t * denom = poly_init_bare(self->order);
	denom->terms = malloc(sizeof(prec_c_t)*self->order);
	fread(denom->terms,sizeof(prec_c_t),self->order,fp);
	// Form the object
	struct Pade_t * eq = pade_init_poly(num,denom);
	eq->offset = K0;
	eq->vals=Roots;
	self->eqs = eq;
	
	fread(&self->curr_t,sizeof(prec_t),1,fp);
	fread(&self->curr_x,sizeof(prec_t),1,fp);
	self->tt = (prec_t*)malloc(sizeof(prec_t)*self->order);
	fread(self->tt,sizeof(prec_t),self->order,fp);
	self->xx = (prec_t*)malloc(sizeof(prec_t)*self->order);
	fread(self->xx,sizeof(prec_t),self->order,fp);
	self->yy = (prec_c_t*)malloc(sizeof(prec_t)*self->order);
	fread(self->yy,sizeof(prec_c_t),self->order,fp);
	switch (self->order) {
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
	self->cb=NULL;
	return self;
	fclose(fp);
}
