#ifndef __LINKED_LIST_H__
#define __LINKED_LIST_H__

#include <string.h>
#include <stdlib.h>

/*
 * Implements essentially a hash map in the form of a linked list. Each element stores the key to associate with it, and the actual underlying contents. During lookup, you pass the key, and it searches the list and returns the associated solver struct.
 */

struct node {
	char * name;
	void * contents;
	struct node * next;
};

extern struct node * root;

/*
 * When a new solver object is created, this will add add it to the linked list
 */
int add(char * name, void * ptr);

/*
 * Perform the lookup and return the underlying stored struct.
 */
void * find_obj(char * name);

#endif
