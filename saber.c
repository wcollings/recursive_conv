#include <stdint.h>
#include "include/sara.h"
#include "include/poly.h"
#include <saberApi.h>
#define JOB inp[0]
#define SETUP 1
#define CALL 2
/* #define K0  inp[1] */
/* #define K1  inp[2] */
/* #define K1i inp[3] */
/* #define K2  inp[4] */
/* #define K2i inp[5] */
/* #define K3  inp[6] */
/* #define K3i inp[7] */
/* #define K4  inp[8] */
/* #define K4i inp[9] */

/* #define s1  inp[10] */
/* #define s1i inp[11] */
/* #define s2  inp[12] */
/* #define s2i inp[13] */
/* #define s3  inp[14] */
/* #define s3i inp[15] */
/* #define s4  inp[16] */
/* #define s4i inp[17] */
/* #define t   inp[1] */
/* #define vl  inp[2] */

struct init_array {
	double k0;
	prec_c_t K1;
	prec_c_t K2;
	prec_c_t K3;
	prec_c_t K4;
	prec_c_t s1;
	prec_c_t s2;
	prec_c_t s3;
	prec_c_t s4;
};
struct work_array {
	double t;
	double vl;
};

struct inpt {
	int64_t sel;
	union input_wrapper{
		struct init_array init;
		struct work_array work;
	} wrap;
};


#if defined(_MSC_VER)
#define ind IND
#endif
SABER_FOREIGN_ROUTINE(ind){
	printf("The size of the first elem is %od",sizeof(inp[0]));
	printf("The size of the second elem is %od",sizeof(inp[1]));
	if (JOB==SETUP) {
		int order = 0;
		for (int i=2; i < 9; i+=2) {
			if (inp[i]==0 && inp[i+1]==0) {
				order = (i-2)>>1;
				break;
			}
		}
		struct Polynomial_t * num = poly_init_bare(order);
		num->terms=malloc(sizeof(prec_c_t)*order);
		struct Polynomial_t * denom = poly_init_bare(order);
		denom->terms=malloc(sizeof(prec_c_t)*order);
		for (int i=2; i<(order+2); i+=2) {
			int index = (i-2)>>1;
		}

	}

}
