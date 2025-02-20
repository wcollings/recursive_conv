struct node {
	char * name;
	int id;
	struct node * left;
	struct node * right;
};

extern struct node * root;
extern int max_id;
struct node * create_node(char * name);
void add_node(struct node * curr,struct node * n);
int find_node(struct node * curr,char * name);

