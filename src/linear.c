#include <stdlib.h>

void rref(double ** mat, int m, int n) {
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
 * `N`: the number of columns in the matrix
 *
*/
void double_sort(double **mat,int m, int n) {
	int *pivots=malloc(sizeof(double)*n);

	// I figured it would be easiest to do this one column at a time, so `temp` was supposed to be
	// the column vector we're currently processing
	int *temp=malloc(sizeof(double)*n);
	for (int i=0;i < n;++i) {
		temp[i]=i;
	}
	// sorting goes here? It will probably need to be O(n^3) - an outer loop for each column vector,
	// and whatever sorting algo you want that will go inside there.
	for (int i=0; i < n; ++i) {
		for (int j=0; j < m; ++j) {
			
		}
	}

	free(temp);
	free(pivots);
}
