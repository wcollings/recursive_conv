#ifndef __LINEAR_H__
#define __LINEAR_H__

#include <stdlib.h>
#include "central.h"

void mat_free(void ** A, int m);
void ** mat_init(int m, int n,size_t size);
void rref(prec_t ** mat, int m, int n);

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
prec_t * mat_mul_vec(prec_t **mat, prec_t * v, int m, int n);
#endif
