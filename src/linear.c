#include <stdlib.h>
void rref(double ** mat, int m, int n) {
}

void double_sort(double **mat,int m, int n) {
	int *pivots=malloc(sizeof(double)*n);
	int *temp=malloc(sizeof(double)*n);
	for (int i=0;i < n;++i) {
		temp[i]=i;
	}
	for (int i=0; i < n; ++i) {
		for (int j=0; j < m; ++j) {
			
		}
	}

	free(temp);
	free(pivots);
}
double * mat_mul_vec(double **mat, double * v, int m, int n) {
	double * res = calloc(m,sizeof(double));
	for (int i=0; i < m; ++i) {
		for (int j=0; j < n; ++j) {
			res[i]+=mat[i][j]*v[j];
		}
	}
	return res;
}
