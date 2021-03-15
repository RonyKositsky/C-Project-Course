#include <math.h>
#include "spmat.h"
#include <stdio.h>
#include <stdlib.h>
#include "eHandler.h"

/* load info provided into the cell: column & data*/
void createCell(cell *cell, int column, double data) {
	cell->col = column;
	cell->value = data;
	cell->nextCell = NULL;
}

/* read a Row from original matrix and store all non-zero elements as linked-list in Matrix*/
void add_row_list(struct _spmat *A, int* indices, int neighbourNumber,
		int row) {
	int firstElement = 1;
	int index, j;
	cell* current, *newCell;
	(A->private_property)[row] = NULL;
	for (j = 0; j < neighbourNumber; j++) {
		newCell = (cell*) calloc(1, sizeof(cell));
		allocCheck(newCell, __FILE__, __LINE__);
		index = indices[j];
		createCell(newCell, index, 1);
		if (firstElement) {
			current = newCell;
			(A->private_property)[row] = current;
			firstElement = 0;
		} else {
			current->nextCell = newCell;
			current = current->nextCell;
		}
	}
}

/*frees all Matrix memory*/
void free_list(struct _spmat *A) {
	int j = 0, n = A->n, i;
	cell* current, *previous;
	for (; j < n; ++j) {
		current = A->private_property[j];
		while (current != NULL) {
			previous = current;
			current = current->nextCell;
			free(previous);
		}
		free(current);
	}

	for (i = 0; i < *(A->numberOfGroups); i++) {
			free(((A->O)[i]->indices));
		}


	for (i = 0; i < n; i++) {
		free((A->P)[i]);
		free((A->O)[i]);
	}

	free(A->P);
	free(A->O);
	free(A->private_property);
	free(A->KFactorList);
	free(A->KList);
	free(A->numberOfGroups);

}

/*multiply matrix A with vector v and saves the vector in result */
void mult_list(const struct _spmat *A, const double *v, double *result) {
	int j = 0, n = A->n, column;
	double data, dotproduct;
	cell* current;
	for (; j < n; ++j) {
		dotproduct = 0.0;
		current = A->private_property[j];
		while (current != NULL) {
			column = current->col;
			data = current->value;
			dotproduct += data * v[column];
			current = current->nextCell;
		}
		result[j] = dotproduct;
	}
}

/*return the value of the cell in a specific place.*/
double actualValue_list(const struct _spmat *A, int row, int column) {
	int size = A->n, i;
	double result;
	cell *cell;
	cell = A->private_property[row];
	if (cell == NULL) {
		return 0;
	}
	for (i = 0; i < size; i++) {
		if (cell) {
			if (cell->col < column) {
				cell = cell->nextCell;
			}

			else if (cell->col == column) {
				result = cell->value;
				return result;
			}

			else {
				break;
			}
		}

		else {
			break;
		}

	}
	return 0;
}

double sumOfRowB_List(const struct _spmat *A, int *indices, int indices_size,
		int row) {
	int i;
	double result = 0;
	for (i = 0; i < indices_size; i++) {
		result += (A->actualValue(A,row, indices[i])
				- (A->KFactorList[row]* A->KList[indices[i]]));
	}
	return result;
}

/*return the 1-norm of the matrix - max sum of absolutes by columns, for matrix shifting*/
void calculateNorma_list(struct _spmat *A) {
	int size = A->n, row, column;
	double max = 0, sum = 0;
	for (column = 0; column < size; column++) {
		for (row = 0; row < size; row++) {
			if ((A->private_property)[row]) {
				sum += fabs(actualValue_list(A, row, column)-(((A->KFactorList)[row])*((A->KList[column]))));
			}
		}
		if (sum > max) {
			max = sum;
		}
		sum = 0;
	}
	A->norma = max;
}

/*hard copy of source matrix to destination matrix*/
void copyMatrix_list(const struct _spmat *source, struct _spmat *destination) {
	int size = source->n, i;
	cell *sourceCell, *destCell;
	for (i = 0; i < size; i++) {
		if (source->private_property[i]) {
			sourceCell = source->private_property[i];
			destination->private_property[i] = destCell; /* Rony - double check*/
			while (sourceCell != NULL) {
				destCell = (cell*) calloc(1, sizeof(cell));
				allocCheck(destCell, __FILE__, __LINE__);
				createCell(destCell, sourceCell->col, sourceCell->value);
				sourceCell = sourceCell->nextCell;
				destCell = destCell->nextCell;
			}
		}
	}
}

/*Generating P and O group.*/
void generateGroups(int size, group** P, group** O) {
	int j;
	allocateGroupList(O, size);
	allocateGroupList(P,size);
	P[0]->indices = (int*) calloc(size, sizeof(int)); /*All the indices starts at P.*/
	allocCheck(P[0]->indices, __FILE__, __LINE__);
	for (j = 0; j < size; j++) {
		(P[0]->indices)[j] = j;
	}

	P[0]->size = size;
}

/* Allocates a new linked-lists sparse matrix of size n */
spmat* allocateMatrix(int n) {
	spmat *matrix = (spmat*) (calloc(1, sizeof(spmat)));
	allocCheck(matrix, __FILE__, __LINE__);
	matrix->free = &free_list;
	matrix->mult = &mult_list;
	matrix->private_property = (cell**) (calloc(n, sizeof(cell*)));
	allocCheck(matrix->private_property, __FILE__, __LINE__);
	matrix->calculateNorma = &calculateNorma_list;
	matrix->copyMatrix = &copyMatrix_list;
	matrix->add_row = add_row_list;
	matrix->numberOfGroups = malloc(sizeof(int));
	allocCheck(matrix->numberOfGroups, __FILE__, __LINE__);
	matrix->actualValue = &actualValue_list;
	matrix->sumOfRowB=&sumOfRowB_List;
	matrix->n = n;
	matrix->M = 0;
	*(matrix->numberOfGroups) = 0;
	matrix->KFactorList=(double*)calloc(n,sizeof(double));
	allocCheck(matrix->KFactorList, __FILE__, __LINE__);
	matrix->KList=(int*)calloc(n,sizeof(int));
	allocCheck(matrix->KList, __FILE__, __LINE__);
	matrix->P = (group**) calloc(n, sizeof(group*));
	allocCheck(matrix->P, __FILE__, __LINE__);
	matrix->O = (group**) calloc(n, sizeof(group*));
	allocCheck(matrix->O, __FILE__, __LINE__);

	return matrix;
}
