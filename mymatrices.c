#include "mymatrices.h"

mtx *matrix_alloc(int n_z, int n_s) {
	mtx *a;
	double *tmp;
	int k, l;
	a = malloc(sizeof(mtx));
	a->cnt = malloc(n_z * sizeof(void *));
	if (a->cnt == NULL)
		return NULL;
	for (k = 0; k < n_z; ++k) {
		tmp = malloc(n_s * sizeof(double));
		if (tmp == NULL)
			break;
		else
			a->cnt[k] = tmp;
	}
	if (k < n_z) {
		for (l = 0; l < k; ++l) {
			free(a->cnt[l]);
			free(a);
		}
		return NULL;
	}
	a->dim_z = n_z;
	a->dim_s = n_s;
	if (n_z == n_s)
		a->sq = 1;
	else
		a->sq = 0;

	return a;
}

void matrix_free(mtx *A) {
	free(A->cnt);
	free(A);
}

void matrix_fprint(char *file, mtx *A) {
	int k, l;
	FILE *fp = fopen(file, "w");
	if (fp == NULL)
		return;
	for (k = 0; k < A->dim_z; ++k) {
		for (l = 0; l < A->dim_s; ++l) {
			fprintf(fp, "%05.3f", A->cnt[k][l]);
			if (l < A->dim_s - 1)
				fprintf(fp, " ");
		}
		if (k < A->dim_z - 1)
			fprintf(fp, "\n");
	}
	fclose(fp);
}

void matrix_print(mtx *A) {
	int k, l;
	for (k = 0; k < A->dim_z; ++k) {
		for (l = 0; l < A->dim_s; ++l)
			printf("%05.3f ", A->cnt[k][l]);
		printf("\n");
	}
}

void matrix_load(char *file, mtx *matrix) {
	int k, l, n_z = 1, n_s = 1;
	char ch;
	FILE *fp = fopen(file, "r");
	if (fp == NULL)
		return;

	while (!feof(fp)) {
		fscanf(fp, "%c", &ch);
		if (ch == '\n')
			++n_z;
		else if ((ch == ' ' || ch == '\t') && n_z == 1)
			++n_s;
	}
	fseek(fp, 0, SEEK_SET);
	for (k = 0; k < n_z; ++k) {
		for (l = 0; l < n_s; ++l) {
			fscanf(fp, "%lf", &(matrix->cnt[k][l]));
		}
	}
}

void matrix_id(mtx *A) {
	int k;
	if (!A->sq)
		return;
	for (k = 0; k < A->dim_s; ++k)
		A->cnt[k][k] = 1.;
}

void matrix_transpo(mtx *A) {
	double tmp;
	int k, l;
	if (!A->sq)
		return;
	for (k = 0; k < A->dim_z; ++k) {
		for (l = k; l < A->dim_z; ++l) {
			tmp = A->cnt[k][l];
			A->cnt[k][l] = A->cnt[l][k];
			A->cnt[l][k] = tmp;
		}
	}
}

void matrix_initup(mtx *A) {
	int k, l;
	for (k = 0; k < A->dim_z; ++k) {
		for (l = 0; l < A->dim_s; ++l)
			A->cnt[k][l] = l;
	}
}

mtx *matrix_mult(mtx *A, mtx *B) {
	int k, l, s;
	if (A->dim_z != A->dim_z || !A->sq || !A->sq)
		return NULL;
	mtx *C = matrix_alloc(A->dim_z, A->dim_s);
	if (C == NULL)
		return NULL;
	for (k = 0; k < C->dim_z; ++k) {
		for (l = k; l < C->dim_z; ++l) {
			C->cnt[k][l] = 0;
			C->cnt[l][k] = 0;
			for (s = 0; s < C->dim_z; ++s) {
				C->cnt[k][l] += A->cnt[k][s] * B->cnt[s][l];
				if (l != k)
					C->cnt[l][k] += A->cnt[l][s] * B->cnt[s][k];
			}
		}
	}
	return C;
}
