#include "../include/linear.h"

void rref(prec_t ** mat, int m, int n) {
}

void prec_t_sort(prec_t **mat,int m, int n) {
	int *pivots=malloc(sizeof(prec_t)*n);
	int *temp=malloc(sizeof(prec_t)*n);
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

prec_t * mat_mul_vec(prec_t **mat, prec_t * v, int m, int n) {
	prec_t * res = calloc(m,sizeof(prec_t));
	for (int i=0; i < m; ++i) {
		for (int j=0; j < n; ++j) {
			res[i]+=mat[i][j]*v[j];
		}
	}
	return res;
}
