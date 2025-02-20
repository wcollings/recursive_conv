#include "include/sara.h"
#include <complex.h>
#include <stdio.h>

int main(int argc, char** argv) {
	if (argc==1) {
		printf("Error - no file name given!\n");
		return -1;
	}
	printf("reading file: %s\n",argv[1]);
	struct Solver_t * solv=read_solver(argv[1]);
	printf("Solver order:%d\n",solv->order);
	printf("Last time step: %lf\n",solv->curr_t);
	printf("Last y components:\n");
/* 	for (int i=0; i < solv->order; ++i) { */
/* 		printf("\t%lf+j%lf\n",creal(solv->yy[i]),cimag(solv->yy[i])); */
/* 	} */
}
