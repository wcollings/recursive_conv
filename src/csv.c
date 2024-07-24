#include <stdlib.h>
#include <stdio.h>
#include "../include/csv.h"
int _fac(int i) {
	if (i==0)
		return 1;
	return i*_fac(i-1);
}

void write_poly(char* fname,struct Polynomial_t * p) {
	FILE * fp = fopen(fname,"w");
	for (int i=0; i < p->num_terms; ++i) {
		fprintf(fp,"%e,",p->terms[i]);
	}
	fclose(fp);
}

/*
 * Read a file that contains a list of derivatives,
 * convert them into a taylor series
 * and return the taylor representation
*/
struct Polynomial_t * read_poly(char* fname) {
	double res=0;
	prec_t * arr = malloc(sizeof(prec_t)*20); // just a large-ish buffer. We know it'll be at max 8 elements
	char proc_buff[80];
	char * endp;
	int i;
	FILE * fp = fopen(fname,"r");
	if (fp==NULL) {
		printf("Could not open file!\n");
		return NULL;
	}
	for (i=0; i < 20; ++i) {
		int read=0;
		#if USE_DOUBLE
			read = fscanf(fp,"%le,",&arr[i]);
		#else
			read = fscanf(fp,"%e,",&arr[i]);
		#endif
		if (read == EOF) {
			break;
		}
		#if DEBUG_PRINTS
			printf("read %le\n",arr[i]/_fac(i));
		#endif
		arr[i]/=(prec_t)_fac(i);
		/* arr[i] = strtof(proc_buff,&endp); ///fac(i); */
	}
	/* i++; */
	flip_arr(arr,i);
	struct Polynomial_t * ret = poly_init_bare(i);
	ret->terms = arr;
	fclose(fp);
	return ret;
}
