#include "si_conv.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

char units[]="qryzafpnum KMGTPEZYRQ";

int index_of_char(char *src,char to_find) {
	int i=0;
	while (src[i] != to_find && src[i] != NULL) {
		i++;
	}
	return (src[i]==to_find?i:-1);
}
/*
 * `str` The converted string to pass in. CANNOT BE NULL, MUST BE NULL-TERMINATED
 */
double from_si(char* str) {
	if (str=="0") {
		return 0.;
	}
	char *num_str;
	char *unit;
	double ret=strtod(str,&unit);
	printf("Number is %lf, unit is %d\n",ret,(int)unit[0]);
	if ((int)unit[0]==0) {
		return ret;
	}
	int idx=index_of_char(units, unit[0]);
	if (idx==-1) {
		return ret;
	}
	int p=(idx-10)*3;
	printf("%c -> 1e%d\n",unit[0],p);
	return ret*pow(10.,p);
}

char * to_si(double in) {
	int sign=(in<0?-1:1);
	in=fabs(in);
	float inter=log10(in);
	int k;
	if (inter > 0) {
		k=3*floor(inter/3);
	} else {
		k=3*ceil(inter/3);
	}
	float new_val=in/pow(10,k);
	char *out=malloc(sizeof(char)*8);
	sprintf(out,"%.2f%c",new_val,units[(k/3)+10]);
	return out;
}

int main() {
	from_si("10K");
	from_si("100");
	from_si("20u");
	/* printf("%lf -> %s\n",1e4,to_si(1e4)); */
	/* printf("%lf -> %s\n",(double)100,to_si(100)); */
	printf("%lf -> %s\n",2e-5,to_si(2e-5));
}
