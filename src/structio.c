#include <stdlib.h>
#include <stdio.h>
#include "../include/csv.h"
#include "../include/sara.h"
#include "poly.h"

void write_solver(struct Solver_t * s,char* fname) {
	FILE * fp = fopen(fname,"wb");
	fwrite(s->order,sizeof(int),1,fp); 
	int num_terms=s->eqs->num->num_terms;
	fwrite(num_terms,sizeof(int),1,fp); 
	fwrite(&s->eqs->offset,sizeof(prec_c_t),1,fp);
	fwrite(s->eqs->num->terms,sizeof(prec_c_t),num_terms,fp);
	fwrite(s->eqs->denom->terms,sizeof(prec_c_t),num_terms,fp);
	fwrite(&s->curr_t,sizeof(prec_t),1,fp);
	fwrite(&s->curr_x,sizeof(prec_t),1,fp);
	fwrite(s->tt,sizeof(prec_t),num_terms,fp);
	fwrite(s->xx,sizeof(prec_t),num_terms,fp);
	fwrite(s->yy,sizeof(prec_c_t),num_terms,fp);
	fclose(fp);
}
struct Solver_t * read_solver(char * fname) {
	FILE *fp = fopen(fname,"rb");
	struct Solver_t * self=malloc(sizeof(struct Solver_t));
	fread(&self->order,sizeof(int),1,fp);
	prec_c_t K0;
	fread(&K0,sizeof(prec_c_t),1,fp);
	// Read in the equation:
	// Ki terms
	struct Polynomial_t * num = poly_init_bare(self->order);
	num->terms = malloc(sizeof(prec_c_t)*self->order);
	fread(&num->terms,sizeof(prec_c_t),self->order,fp);
	//sigma_i terms
	struct Polynomial_t * denom = poly_init_bare(self->order);
	denom->terms = malloc(sizeof(prec_c_t)*self->order);
	fread(&denom->terms,sizeof(prec_c_t),self->order,fp);
	// Form the object
	struct Pade_t * eq = pade_init_poly(num,denom);
	eq->offset = K0;
	eq->vals=Roots;
	self->eqs = eq;
	
	fread(&self->curr_t,sizeof(prec_t),1,fp);
	fread(&self->curr_x,sizeof(prec_t),1,fp);
	self->tt = (prec_t*)malloc(sizeof(prec_t)*self->order);
	fread(self->tt,sizeof(prec_t),self->order,fp);
	self->xx = (prec_t*)malloc(sizeof(prec_t)*self->order);
	fread(self->xx,sizeof(prec_t),self->order,fp);
	self->yy = (prec_c_t*)malloc(sizeof(prec_t)*self->order);
	fread(self->yy,sizeof(prec_c_t),self->order,fp);
	switch (self->order) {
		case 1: self->qq=q1;
				  break;
		case 2: self->qq=q2;
				  break;
		case 3: self->qq=q3;
				  break;
		case 4: self->qq=q4;
				  break;
		default: self->qq=q2;
	}
	return self;
}
