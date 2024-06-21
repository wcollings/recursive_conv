#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../include/ddiff.h"
prec_t * divided_difference(prec_t * x,prec_t *y, int num_ele) {
	prec_t * out=malloc(sizeof(prec_t)*num_ele);
	for (int i=0; i < num_ele; ++i)
		out[i]=y[i];
	for (int i=1; i < num_ele+1; ++i) {
		for (int j=num_ele-1; j > i; j--) {
			prec_t num=out[j]-out[j-1];
			prec_t den=x[j]-x[j-i];
			out[j]=num/den;
#if DEBUG_PRINTS
			printf("i=%d j=%d d=%e\n",i,j,out[j]);
#endif
		}
	}
	out[num_ele-1]=(out[num_ele-1]-out[num_ele-2])/(x[num_ele-1]-x[0]);
	return out;
}

int64_t fac(int x) {
	return (x==0?1:x*fac(x-1));
}

void scale_to_taylor(prec_t * terms, int num_ele) {
	for (int i=0; i < num_ele; ++i) {
		terms[i] = terms[i]/fac(i);
	}
}

void print_mac(float *c, int numele) {
	printf("[");
	for (int i=0; i < numele; ++i) {
		float coeff=c[i]/fac(i);
		printf("%e,",coeff);
	}
	printf("]\n");
}

int test_ddiff() {
	prec_t * x=malloc(sizeof(prec_t)*10);
	prec_t * y=malloc(sizeof(prec_t)*10);
	float bot=2e6;
	for (int i=0; i < 10; ++i) {
		x[i]=bot+(i*20);
		y[i]=x[i];
	}
	prec_t * dd=divided_difference(x, y, 10);
	printf("x\t\ty\t\tdd\n");
	for (int i=0; i < 10; ++i) {
		printf("%e\t%e\t%e\n",x[i],y[i],dd[i]);
	}
	return 0;
}
