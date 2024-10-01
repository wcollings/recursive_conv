#ifndef __LOG_H__
#define __LOG_H__
#include <stdio.h>
#include <stdlib.h>

FILE * lf;
void log_init(char * log_name);
void log_warn(char * msg);
void log_err(char * msg);
void log_msg(char * msg);
void log_close();

#endif
