#include <stdlib.h>
#include <math.h>

void rref(double** mat, int m, int n) {
	int pivotCoeff, otherCoeff;
	for (int col = 0; col < m; ++col)
	{
		for (int row = 0; row < m; ++row)
		{
			if (row != col && mat[row][col] != 0)
			{
				pivotCoeff = mat[col][col];
				otherCoeff = mat[row][col];
				for (int k = 0; k < n; ++k)
				{
					mat[row][k] = (pivotCoeff * mat[row][k]) - (otherCoeff * mat[col][k]);
				}
			}
		}
	}
	for (int row = 0; row < m; ++row)
	{
		mat[row][n - 1] = mat[row][n - 1] / mat[row][row];
		mat[row][row] = 1;
	}
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
 * meow :3
 *
*/
void double_sort(double** mat, int m, int n) {

	// I figured it would be easiest to do this one column at a time, so `temp` was supposed to be
	// the column vector we're currently processing
	int* maxColElemIndex = malloc(sizeof(int)*m);						// index = col num; value at index = index of row with highest value of that col
	double maxColElem;							// keep track of largest elem in current column
	double* temp = malloc(sizeof(double)*m);	// used to rearrange matrix rows by changing pointers
	for (int col = 0; col < n - 1; ++col) {
		maxColElem = fabs(mat[0][col]);
		maxColElemIndex[col] = 0;
		for (int row = 0; row < m; ++row)
		{
			temp[col] = mat[row][col];
			if (fabs(temp[col]) > maxColElem)
			{
				maxColElem = fabs(mat[row][col]);			// store highest value in column, i.e pivot variable
				maxColElemIndex[col] = row;						// store the row index of the max element, so we can rearrange the matrix accordingly
			}
		}
	}

	double ** sortedMatrix = malloc(sizeof(double*)*m);
	for (int i = 0; i < m; ++i)
	{
		sortedMatrix[i] = malloc(sizeof(double) * n);
	}

	for (int i = 0; i < m; ++i) {
		sortedMatrix[i] = mat[maxColElemIndex[i]];
	}

	for (int i = 0; i < m; ++i) {
		mat[i] = sortedMatrix[i];
	}
	//mat = sortedMatrix;


	free(temp);
}
