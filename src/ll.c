#include "../include/ll.h"
struct node * root=NULL;

struct node * make(char * name, void * ptr) {
	struct node * new_elem=malloc(sizeof(struct node*));
	int len=strlen(name)+1;
	new_elem->name=malloc(sizeof(char)*len);
	strncpy(new_elem->name,name,len);
	new_elem->contents=ptr;
	new_elem->next=NULL;
	return new_elem;
}
int add(char * name, void * ptr) {
	struct node * curr=root;
	if (root==NULL) {
		root=make(name,ptr);
	} else {
		do {
			if (!strcmp(name,curr->name))
				return -1;
			curr=curr->next;
		}while (curr->next != NULL);
		curr->next=make(name,ptr);
	}
	return 0;
}
void * find_obj(char * name) {
	struct node * curr=root;
	if (!strcmp(name,root->name))
		return root->contents;
	do {
		if (!strcmp(name,curr->name))
			return curr->contents;
		curr=curr->next;
	} while (curr != NULL);
	return NULL;
}
