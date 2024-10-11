#include "../include/log.h"

void log_init(char * log_name){
	lf=fopen(log_name,"a");
	fprintf(lf,"Log opened\n");
}

void log_warn(char * msg) {
	fprintf(lf,"WARNING: %s",msg);
}
void log_err(char * msg) {
	fprintf(lf,"ERROR: %s",msg);
}
void log_msg(char * msg) {
	fprintf(lf,"MESSAGE: %s",msg);
}
void log_close() {
	fclose(lf);
}


