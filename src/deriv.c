#include "../include/deriv.h"
#include <stdio.h>
#define DIFF(h) (*fn)(x+h) - (*fn)(x-h)
#define SUM(h) (*fn)(x+h) + (*fn)(x-h)


prec_t sten_1(prec_t (*fn)(prec_t),prec_t x,prec_t h) {
	int terms[] = {45,9,1};
	prec_t tot = 0;
	for (int i=0; i < 3; ++i) {
		tot += terms[i]*DIFF(h*(i+1));
	}
	return tot/(60*h);
}

prec_t sten_2(prec_t (*fn)(prec_t),prec_t x,prec_t h) {
	int terms[] = {25,44,-9};
	prec_t tot = -120*(*fn)(x);
	for (int i=0; i < 3; ++i) {
		tot += terms[i]*SUM(h*(i+1));
	}
	return tot/(120*pow(h,2));
}

prec_t sten_3(prec_t (*fn)(prec_t),prec_t x,prec_t h) {
	int terms[] = {-13,8,-1};
	prec_t tot = 0;
	for (int i=0; i < 3; ++i) {
		tot += terms[i]*DIFF(h*(i+1));
	}
	return tot/(8*pow(h,3));
}
prec_t sten_4(prec_t (*fn)(prec_t),prec_t x,prec_t h) {
	int terms[] = {899,-1612,837,-124};
	prec_t tot = 0;
	for (int i=0; i < 4; ++i) {
		tot += terms[i]*SUM(h*(i+1));
	}
	return tot/(930*pow(h,4));
}
prec_t sten_5(prec_t (*fn)(prec_t),prec_t x,prec_t h) {
	int terms[] = {5,4,1};
	prec_t tot = 0;
	for (int i=0; i < 3; ++i) {
		tot += terms[i]*DIFF(h*(i+1));
	}
	return tot/(2*pow(h,5));
}
prec_t sten_6(prec_t (*fn)(prec_t),prec_t x,prec_t h) {
	int terms[] = {-217,434,-279,62};
	prec_t tot = 0;
	for (int i=0; i < 4; ++i) {
		tot += terms[i]*SUM(h*(i+1));
	}
	return tot/(217*pow(h,6));
}
prec_t sten_7(prec_t (*fn)(prec_t),prec_t x,prec_t h) {
	int terms[] = {-14,14,-6,1};
	prec_t tot = 0;
	for (int i=0; i < 4; ++i) {
		tot += terms[i]*DIFF(h*(i+1));
	}
	return tot/(2*pow(h,7));
}
prec_t * take_derivatives(prec_t (*fn)(prec_t), prec_t x, prec_t h) {
	prec_t * res = malloc(sizeof(prec_t)*8);
	res[0] = (*fn)(x);
	res[1] = sten_1(fn,x,h);
	res[2] = sten_2(fn,x,h);
	res[3] = sten_3(fn,x,h);
	res[4] = sten_4(fn,x,h);
	res[5] = sten_5(fn,x,h);
	res[6] = sten_6(fn,x,h);
	res[7] = sten_7(fn,x,h);
	for (int i=0; i < 8; ++i) {
		printf("%i->%1.3e\n",i,res[i]);
	}
	return res;
}

int64_t fac(int x) {
	return (x==0?1:x*fac(x-1));
}

void scale_to_taylor(prec_t * terms, int num_ele) {
	for (int i=0; i < num_ele; ++i) {
		terms[i] = terms[i]/fac(i);
		printf("%i->%1.3e\n",i,terms[i]);
	}
}
