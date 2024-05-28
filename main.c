#include <stdlib.h>
#include <stdio.h>
#include "include/pade.h"
#include "include/sara.h"

void print_results(struct solver_t * SOLV) {
	float output=SOLV->yy[0][0] + SOLV->yy[1][0];
	printf("%f,%f\n",SOLV->curr_t,output);
}

int main() {
	struct Pade_t * p;
	p=malloc(sizeof(struct Pade_t));
	p->M=p->N=2;
	p->A=malloc(sizeof(float)*2);
	p->A[0]=0.9478;
	p->A[1]=0.0522;
	p->B=malloc(sizeof(float)*2);
	p->B[0]=-2.105;
	p->B[1]=-0.095;
	struct solver_t * solv = init_solver(2, p);
	solv->cb=&print_results;
	float time=0;
	printf("t,v\n");
	for (int i=0; i < 40; ++i) {
		step(solv,1,i);
	}
	return 0;
	/* double res=0; */
	/* double inpts[2]={1,2}; */
	/* printf("TIME\tV\n"); */
	/* for (int i=0; i < 40; ++i) { */
	/* 	inpts[0]+=1; */
	/* 	res=result(2,inpts); */
	/* 	printf("%e\t%f\n",inpts[0],res); */
	/* } */
	/* return 0; */
}
