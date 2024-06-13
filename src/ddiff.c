#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../include/ddiff.h"
double * divided_difference(double * x,double *y, int num_ele) {
	double * out=malloc(sizeof(double)*num_ele);
	for (int i=0; i < num_ele; ++i)
		out[i]=y[i];
	for (int i=1; i < num_ele+1; ++i) {
		for (int j=num_ele-1; j > i; j--) {
			double num=out[j]-out[j-1];
			double den=x[j]-x[j-i];
			out[j]=num/den;
			printf("i=%d j=%d d=%e\n",i,j,out[j]);
		}
	}
	out[num_ele-1]=(out[num_ele-1]-out[num_ele-2])/(x[num_ele-1]-x[0]);
	return out;
}

int64_t fac(int x) {
	return (x==0?1:x*fac(x-1));
}

void scale_to_taylor(double * terms, int num_ele) {
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
	double * x=malloc(sizeof(double)*10);
	double * y=malloc(sizeof(double)*10);
	float bot=2e6;
	for (int i=0; i < 10; ++i) {
		x[i]=bot+(i*20);
		y[i]=L(x[i]);
	}
	double * dd=divided_difference(x, y, 10);
	printf("x\t\ty\t\tdd\n");
	for (int i=0; i < 10; ++i) {
		printf("%e\t%e\t%e\n",x[i],y[i],dd[i]);
	}
}
