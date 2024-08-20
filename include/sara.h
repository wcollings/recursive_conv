#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "pade.h"
#include "central.h"

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
 * To solve these, we really need the time deltas, not the times themselves. To save computation,
 * we can compute that once and save it, then save the last actual time so that we can compute
 * the next time delta as well going forward.
 *
*/
struct Solver_t {
	struct Pade_t * eqs;
	prec_t curr_t; /* The last _actual_ time */
	prec_t curr_x; /* The last _actual_ value of x */
	int order;
	prec_t * tt; /* Past time deltas */
	prec_t * xx; /* previous input states */
	prec_c_t * yy; /* previous output states */
	prec_c_t (*qq)(float,float,int);
	void (*cb)(struct Solver_t *); /* A callback function (optional) for printing intermediate results etc.*/
};

prec_c_t q1(float,float,int);
prec_c_t q2(float,float,int);
prec_c_t q3(float,float,int);
prec_c_t q4(float,float,int);

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
void shift(prec_t * arr,int num_ele);

/*
 * Take a step in the time domain, advancing the solver by one time step
 * `SOLV`: an instance of the solver
 * `inpt`: The state variable (likey voltage or current) at the next time step
*/
void step(struct Solver_t * SOLV, prec_t inpt, float curr_t);
/*
 * Params:
 *  - current time
 *  - voltage across part
*/
prec_t result(int argc, prec_t *argv);

#endif
