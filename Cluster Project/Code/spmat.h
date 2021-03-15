#ifndef _SPMAT_H
#define _SPMAT_H
#include "group.h"
#include <stdio.h>

typedef struct cell {
	int col;
	double value;
	struct cell *nextCell;
} cell;

typedef struct _spmat {
	/* Matrix size (n*n) */
	int n;

	/*M in the algorithm*/
	int M;

	/*Matrix 1-norma, calculate once at start*/
	double norma;

	/*Degree of all nodes.*/
	int *KList;

	/*Ki/M*/
	double *KFactorList;

	/*Number of division*/
	int *numberOfGroups;

	/*Group P*/
	group **P;

	/*Group O*/
	group **O;

	/* Adds row i the matrix. Called before any other call,
	 * exactly n times in order (i = 0 to n-1) */
	void (*add_row)(struct _spmat *A, int* indices, int neighbourNumber,
			int row);

	/* Frees all resources used by A */
	void (*free)(struct _spmat *A);

	/* Multiplies matrix A by vector v, into result (result is pre-allocated) */
	void (*mult)(const struct _spmat *A, const double *v, double *result);

	/* Multiplies matrix A by int vector v, into result (result is pre-allocated) */
	void (*multInt)(const struct _spmat *A, int *v, double *result);

	/*copy source matrix to destination matrix*/
	void (*copyMatrix)(const struct _spmat *source, struct _spmat *destination);

	/*calculates the 1-norm of the matrix A, putting in norma*/
	void (*calculateNorma)(struct _spmat *A);

	/*return the value of the matrix in cell [row,column]*/
	double (*actualValue)(const struct _spmat *A, int row, int column);

	/*return the sum of row in the sub matrix according to indices*/
	double (*sumOfRowB)(const struct _spmat *A, int *indices, int indices_size,
				int row);

	/* Private field for inner implementation.
	 * Should not be read or modified externally */
	cell** private_property;

} spmat;

/* Allocates a new linked-lists sparse matrix of size n */
spmat* allocateMatrix(int n);

void generateGroups(int size, group **P, group **O);

#endif
