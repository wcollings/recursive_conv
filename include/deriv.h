#include "../include/central.h"
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

/* Compute the first derivative of the function `fn` */
prec_c_t sten_1(prec_c_t (*fn)(prec_t),prec_t x,prec_t h);
prec_c_t sten_2(prec_c_t (*fn)(prec_t),prec_t x,prec_t h);
prec_c_t sten_3(prec_c_t (*fn)(prec_t),prec_t x,prec_t h);
prec_c_t sten_4(prec_c_t (*fn)(prec_t),prec_t x,prec_t h);
prec_c_t sten_5(prec_c_t (*fn)(prec_t),prec_t x,prec_t h);
prec_c_t sten_6(prec_c_t (*fn)(prec_t),prec_t x,prec_t h);
prec_c_t sten_7(prec_c_t (*fn)(prec_t),prec_t x,prec_t h);
/* Take all 7 derivatives. Returns [f,f',f'',...] (8 elements in all) */
prec_c_t * take_derivatives(prec_c_t (*fn)(prec_t), prec_t x, prec_t h);

/*
 * Take a list of derivatives, and turn that into a Taylor sequence.
 *
 * T(f(x)) = sum(i=0->N) (f^(i)(a)(x-a))/i!)
 */
void scale_to_taylor(prec_c_t * terms, int num_ele);
