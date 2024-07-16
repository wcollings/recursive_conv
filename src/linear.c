#include "../include/linear.h"

void mat_free(prec_t ** A, int m) {
	if(A == 0) {
			printf("uhhh smth broke");
			return;
		}
	for (int i=0; i < m; ++i ) {
		free(A[i]);
	}
	free(A);
}

prec_t ** mat_init(int m, int n,size_t size) {
	prec_t ** inner = malloc(sizeof(void*)*m);
	for (int i=0; i < m; ++i) {
		inner[i] = malloc(size*n);
	}
	return inner;
}

void double_sort(prec_t ** mat, int m, int n);
prec_t GreatestCommonDivisor (prec_t a, prec_t b);
prec_t LowestCommonMultiple (prec_t a, prec_t b);

void rref(prec_t ** mat, int m, int n) {
	double_sort(mat, m, n);
	prec_t pivotCoeff, otherCoeff, lcm;
	for (int col = 0; col < m; ++col) {
		for (int row = 0; row < m; ++row) {
			if (row != col && mat[row][col] != 0) {
				lcm = LowestCommonMultiple(mat[col][col],mat[row][col]);
				pivotCoeff = lcm/mat[row][col];
				otherCoeff = lcm/mat[col][col];
				for (int k = 0; k < n; ++k) {
					mat[row][k] = (pivotCoeff * mat[row][k]) - (otherCoeff * mat[col][k]);
				}
			}
		}
	}
	for (int row = 0; row < m; ++row) {
		mat[row][n - 1] = mat[row][n - 1] / mat[row][row];
		mat[row][row] = 1;
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 5; j++) {
			printf("%f ",mat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

/*
 * This is supposed to do the sorting of the column vectors, note that the matrix object itself is
 * passed by reference, so any changes made to it in the process of sorting are perminent.
 * Hence why the function returns void.
 *
 * No need for it to be exposed and public, so it's not in the header, but if you think it should
 * be anyways, feel free to add it there.
 *
 * Should sort the matrix such that it will then reduce directly to the identity matrix, i.e.
 * all the variables to isolate are on the diagonal.
 *
 * Parameters:
 * `mat`: the matrix to sort, which should be Mx(M+1) in size - brings it down to the identity matrix
 * and one more column vector with the answers
 * `M`: the number of rows in the matrix
 * `N`: the number of columns in the matri
 *
*/
void double_sort(prec_t ** mat, int m, int n) {

	// I figured it would be easiest to do this one column at a time, so `temp` was supposed to be
	// the column vector we're currently processing
	int * maxColElemIndex = malloc(sizeof(int)*m);						// index = col num; value at index = index of row with highest value of that col
	prec_t maxColElem;							// keep track of largest elem in current column
	prec_t temp;
	int * flagArray = calloc(m, sizeof(int));	// used to rearrange matrix rows by changing pointers
	int nextFreeRow = 0;
	for (int col = 0; col < n - 1; ++col) {
		while(flagArray[nextFreeRow] == 1)
		{
			nextFreeRow++;
		}
		maxColElem = fabs(mat[nextFreeRow][col]);
		maxColElemIndex[col] = nextFreeRow;
		for (int row = nextFreeRow; row < m; ++row) {	
			if(flagArray[row]!=0) {
				continue;
			}
			temp = mat[row][col];
			if (fabs(temp) >= fabs(maxColElem)) {
				flagArray[maxColElemIndex[col]] = 0;
				maxColElem = fabs(mat[row][col]);			// store highest value in column, i.e pivot variable
				maxColElemIndex[col] = row;						// store the row index of the max element, so we can rearrange the matrix accordingly
				flagArray[row] = 1;
			}
		}
		nextFreeRow = 0;
	}
	
	prec_t ** sortedMatrix = mat_init(m,n,sizeof(prec_t)); 
	for (int i = 0; i < m; ++i) {
		sortedMatrix[i] = mat[maxColElemIndex[i]];
	}

	for (int i = 0; i < m; ++i) {
		mat[i] = sortedMatrix[i];
	}

	//free(sortedMatrix);
	free(maxColElemIndex);
	free(flagArray);
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

prec_t GreatestCommonDivisor (prec_t a, prec_t b) {
	/* How do you feel about this code snippet?
	 * It's logically exactly the same, but I can't decide if it's easier or harder to read...
	 *
	return (a ?
					( b ?
						GreatestCommonDivisor(b, fmod(a,b))
					: a)
				: b );
	*/
	if(a == 0) {
		return b;
	} else if(b == 0) {
		return a;
	} else {
		return GreatestCommonDivisor(b,fmod(a,b));
	}
}

prec_t LowestCommonMultiple (prec_t a, prec_t b) {
	return ((a*b)/GreatestCommonDivisor(a,b));
}
