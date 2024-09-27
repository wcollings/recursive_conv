#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "pade.h"
#include "central.h"
#include <stdint.h>

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
 * `order`: int
 * `eqs`: Pade_t*
 * `curr_t`: prec_t
 * `curr_x`: prec_t
 * `tt`: prec_t[order]
 * `xx`: prec_t[order]
 *	`yy`: prec_c_t[order]
 *	`qq`: prec_c_t(*)(prec_c_t,prec_c_t,int)
 *	`cb`: void(*)(struct Solver_t*)
*/
struct Solver_t {
	int32_t order;
	int16_t num_calls;
	int16_t num_steps;
	struct Pade_t * eqs;
	prec_t curr_t; /* The last _actual_ time */
	prec_t curr_x; /* The last _actual_ value */
	prec_t * tt; /* Past time deltas */
	prec_t * xx; /* previous input states */
	prec_c_t * yy; /* previous output states */
	prec_c_t (*qq)(prec_c_t,prec_c_t,int);
	void (*cb)(struct Solver_t *, double res); /* A callback function (optional) for printing intermediate results etc.*/
};

prec_c_t q1(prec_c_t,prec_c_t,int);
prec_c_t q2(prec_c_t,prec_c_t,int);
prec_c_t q3(prec_c_t,prec_c_t,int);
prec_c_t q4(prec_c_t,prec_c_t,int);

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
prec_t step(struct Solver_t * SOLV, prec_t inpt, prec_t curr_t);
/*
 * Params:
 *  - current time
 *  - voltage across part
*/
prec_t do_step(prec_t inpt, prec_t curr_t);
prec_t do_accept(prec_t inpt, prec_t curr_t);
void write_solver(struct Solver_t * s,char* fname);
struct Solver_t * read_solver(char * fname);

#endif
