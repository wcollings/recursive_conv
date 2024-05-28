#ifndef __SOLVER_H__
#define __SOLVER_H__
#include "pade.h"

struct time_s {
	const float T; /* The final time of the simulation */
	float curr_t;
	const float delta_n; /* The time step of the simulation */
};

void step_time(struct time_s * t);
struct solver_t {
	struct Pade_t * eqs;
	double curr_t;
	int order;
	double * tt;
	double * xx;
	double ** yy;
	float * (*q)(float,float);
	void (*cb)(struct solver_t *); /* A callback function (optional) for printing intermediate results etc.*/
};

float * q1(float,float);
float * q2(float,float);
float * q3(float,float);
float * q4(float,float);
struct solver_t * init_solver(int order,struct Pade_t * eq);
void step(struct solver_t * SOLV, double inpt, float curr_t);
/*
 * Params:
 *  - current time
 *  - voltage across part
 *
*/
double result(int argc, double *argv);

#endif
