#ifndef __CSV_H__
#define __CSV_H__

#include "poly.h"
#include "central.h"

void write_poly(char* fname,struct Polynomial_t * p);
struct Polynomial_t * read_poly(char* fname);

#endif
