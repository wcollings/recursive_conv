#ifndef __LINEAR_H__
#define __LINEAR_H__
void rref(double ** mat, int m, int n);

/*
 * multiply a matrix by a vector, with the vector on the right
 * If the matrix is mxn, then the vector should be nx1.
 * If there are extra entries then this will ignore them,
 * though according to math (but who likes that anyways?)
 * it should say "that's not valid!".
 *
 * returns an mx1 column vector 
 *
*/
double * mat_mul_vec(double **mat, double * v, int m, int n);
#endif
