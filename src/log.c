#include "../include/log.h"

FILE * lf=NULL;
void log_init(char * log_name){
	if (lf==NULL) {
		lf=fopen(log_name,"a");
		fprintf(lf,"Log opened\n");
	}
}

void log_warn(char * msg) {
	fprintf(lf,"WARNING: %s\n",msg);
}
void log_err(char * msg) {
	fprintf(lf,"ERROR: %s\n",msg);
}
void log_msg(char * msg) {
	fprintf(lf,"MESSAGE: %s\n",msg);
}
void log_close() {
	if (lf != NULL) {
		fclose(lf);
		lf=NULL;
	}
}


