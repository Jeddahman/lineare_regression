#ifndef _MYMATRICES_H_
#define _MYMATRICES_H_

#include <stdio.h>
#include <stdlib.h>

struct MATRIX{
	double **cnt;
	int dim_z;
	int dim_s;
	int sq;
};

typedef struct MATRIX mtx;

mtx *matrix_alloc(int n_z,int n_s);
void matrix_free(mtx *A);
void matrix_fprint(char *file, mtx *A);
void matrix_print(mtx *A);
void matrix_load(char *file,mtx *matrix);
void matrix_id(mtx *A);
void matrix_transpo(mtx *A);
void matrix_initup(mtx *A);
mtx *matrix_mult(mtx *A, mtx *B);

#endif
