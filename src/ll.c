#include "../include/ll.h"
struct node * root=NULL;

struct node * make(char * name, void * ptr) {
	struct node * new_elem=malloc(sizeof(struct node*));
	int len=strlen(name);
	new_elem->name=malloc(sizeof(char)*len);
	strncpy(name,new_elem->name,len);
	new_elem->contents=ptr;
	return new_elem;
}
int add(char * name, void * ptr) {
	int found=0;
	struct node * curr=root;
	if (root==NULL) {
		root=make(name,ptr);
	} else {
		while (curr->next != NULL) {
			if (!strcmp(name,curr->name))
				return -1;
			curr=curr->next;
		}
		curr->next=make(name,ptr);
	}
	return 0;
}
void * find_obj(char * name) {
	struct node * curr=root;
	while (curr->next ==NULL) {
		if (!strcmp(name,curr->name))
			return curr->contents;
	}
	return NULL;
}
