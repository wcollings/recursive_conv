#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "pade.h"

struct Time_s {
	const float T; /* The final time of the simulation */
	float curr_t;
	const float delta_n; /* The time step of the simulation */
};
void step_time(struct Time_s * t);

/*
 * Holds all the information needed to compute the next value.
 * This includes history, solver order, current and past time steps, etc.
 *
*/
struct Solver_t {
	struct Pade_t * eqs;
	double curr_t;
	int order;
	double * tt;
	double * xx;
	double ** yy;
	float * (*q)(float,float);
	void (*cb)(struct Solver_t *); /* A callback function (optional) for printing intermediate results etc.*/
};

float * q1(float,float);
float * q2(float,float);
float * q3(float,float);
float * q4(float,float);

/*
 * Initialize a solver object
 * `order` The order of the solver (1, 2, or 3). Least precise (but fastest) to most precise (but slowest)
 * `time` The stop time of the simulation
 * `delta_n` The time step to use for the simulation
 * `num_outputs` The number of output variables to save
 * `returns` a fully initialized Solver object
*/
struct Solver_t * solver_init(int order,struct Pade_t * eq);

/*
 * Shift an array over by one, essentially dequeueing the last element
 * `arr`: the array to shift
 * `num_ele`: the size of the array
*/
void shift(double * arr,int num_ele);

/*
 * Take a step in the time domain, advancing the solver by one time step
 * `SOLV`: an instance of the solver
 * `inpt`: The state variable (likey voltage or current) at the next time step
*/
void step(struct Solver_t * SOLV, double inpt, float curr_t);
/*
 * Params:
 *  - current time
 *  - voltage across part
*/
double result(int argc, double *argv);

#endif
