/*
 * Define a bunch of terms that can need to be coordinated between all parts of the lib
 *
 */
#ifndef __CENTRAL_H__
#define __CENTRAL_H__
#include <complex.h>

#define USE_DOUBLE 0
#if USE_DOUBLE
#define prec_t double
#define PRNT_SPEC 1.3le
#else
#define prec_t float
#define PRNT_SPEC "1.3e"
#endif
#define prec_c_t prec_t _Complex
#define DEBUG_PRINTS 1

#endif
