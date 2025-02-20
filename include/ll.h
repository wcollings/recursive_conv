#ifndef __LINKED_LIST_H__
#define __LINKED_LIST_H__

#include <string.h>
#include <stdlib.h>
struct node {
	char * name;
	void * contents;
	struct node * next;
};

extern struct node * root;
int add(char * name, void * ptr);
void * find_obj(char * name);

#endif
