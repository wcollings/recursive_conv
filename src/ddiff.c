#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
/*
 * @brief Creates a divided difference table based on the x and y given.
 * This essentially approximates `n-1` derivatives of y(x).
 * out[i] -> the ith derivative of this function
 *
 * `x` The inputs to the function (must be `num_ele` elements long)
 * `y` the outputs of the function (must be `num_ele` elements long)
 * `num_ele` the number of elements in x and y. Also corresponds to greater than the number of derivatives you want to take
 * Returns:
 * a table of derivatives (doubles) of length `num_ele`
*/
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

void print_mac(float *c, int numele) {
	printf("[");
	for (int i=0; i < numele; ++i) {
		float coeff=c[i]/fac(i);
		printf("%e,",coeff);
	}
	printf("]\n");
}


double L(float f) {
	float a=1e-9,
			b=2.8e-9,
			c=800e-9,
			f0=2e4;
	double res=(0.6366*a)*atan(-c*(f-f0))+b;
	return res;
}

float * pade(float *C, int num_ele, int m) {
	float ** A=malloc(sizeof(float*)*m);
	for (int i=0; i < m; ++i) A[i]=malloc(sizeof(float)*(m+2));
	for (int i=num_ele-1; i > m-1; --i) {
		for (int j=0; j < m+1; ++j) {
			A[m-i][j]=(i>=j? C[i-j]:0);
		}
	}
}
int main() {
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
