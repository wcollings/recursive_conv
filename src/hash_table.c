#include "../include/hash_table.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

struct node * root;
int max_id;
void add_node(struct node * curr,struct node * n) {
	int place = strcmp(curr->name,n->name);
	if (place > 0) {
		if (curr->left != NULL) {
			add_node(curr->left,n);
			return;
		} else {
			curr->left=n;
		}
	} else if (place < 0) {
		if (curr->right != NULL) {
			add_node(curr->right,n);
			return;
		} else {
			curr->right=n;
		}
	} else {
		fprintf(stderr,"An object of name %s has already been registered; continuing.\n",n->name);
		return;
	}
}
struct node * create_node(char * name) {
	struct node * temp = malloc(sizeof(struct node));
	temp->name=name;
	temp->id=max_id++;
	if (root==NULL) {
		root=temp;
		printf("Registering %s as root node.\n",name);
	}
	return temp;
}
int find_node(struct node * curr, char * name) {
	int place = strcmp(curr->name,name);
	if (place==0) {
		return curr->id;
	} else if (place > 0) {
		if (curr->left==NULL) {
			return -1;
		}
		return find_node(curr->left,name);
	} else {
		if (curr->right==NULL) {
			return -1;
		}
		return find_node(curr->right,name);
	}
}

/* int main() { */
/* 	char *name[]={"part1","part2","part3","part2"}; */
/* 	for (int i=0; i < 4; ++i) { */
/* 		struct node * curr=create_node(name[i]); */
/* 		add_node(root,curr); */
/* 	} */
/* 	for (int i=0; i < 4; ++i) { */
/* 		int id=find_node(root,name[i]); */
/* 		printf("Name = %s, id = %d\n",name[i],id); */
/* 	} */
/* } */
