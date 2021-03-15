#include <math.h>
#include "spmat.h"
#include "algorithms.h"
#include "eHandler.h"
#include <stdio.h>
#include <stdlib.h>

/*Writing the results.*/
void writeResult(FILE* output, spmat* Matrix) {
	int i, j;
	int size = *(Matrix->numberOfGroups);

	fwrite(&size, sizeof(int), 1, output);
	for (i = 0; i < size; i++) {
		fwrite(&((Matrix->O)[i]->size), sizeof(int), 1, output);
		for (j = 0; j < ((Matrix->O)[i]->size); j++) {
			fwrite(&((Matrix->O)[i]->indices)[j], sizeof(int), 1, output);
		}
	}
}

/*Generate the matrix data from output file.*/
void generateMatrix(spmat* Matrix, FILE* input, int size){
	int j, i;
	int non; /*NumberOfNeighbors*/
	int *ind;
	for (j = 0; j < size; j++) {

		fread(&non, sizeof(int), 1, input);
		if (non != 0) {
			(Matrix->KList)[j] = non;
			ind = (int*) calloc(non, sizeof(int));
			allocCheck(ind, __FILE__, __LINE__);
			for (i = 0; i < non; i++) {
				fread(&ind[i], sizeof(int), 1, input);
			}
			Matrix->add_row(Matrix, ind, non, j);
			free(ind);
		}
	}

	for (j = 0; j < size; j++) {
		(Matrix->M) += (Matrix->KList)[j]; /* Calculate factor list.*/
	}
	checkM(Matrix->M);
	for (j = 0; j < size; j++) {
		if(Matrix->M==0){
			(Matrix->KFactorList)[j]=0;
		}else{
			(Matrix->KFactorList)[j] = (double) ((Matrix->KList)[j]) / (Matrix->M); /* Calculate factor list.*/
		}
	}

	generateGroups(Matrix->n,Matrix->P,Matrix->O);
	Matrix->calculateNorma(Matrix);
}

int main(int argc, char *argv[]) {
	FILE *input, *output;
	int size;
	spmat *Matrix;
	double *modularity = malloc(sizeof(double));
	allocCheck(modularity, __FILE__, __LINE__);

	input = fopen(argv[1], "rb"); /* Add error handling.*/
	fread(&size, sizeof(int), 1, input); /*Retrieve number of nodes.*/

	/*Create matrix*/
	Matrix = allocateMatrix(size);
	generateMatrix(Matrix, input, size);
	fclose(input);

	algorithmThree(Matrix, Matrix->P, Matrix->O, Matrix->numberOfGroups,modularity);

	output = fopen(argv[argc - 1], "wb"); /*Add error handling.*/
	writeResult(output, Matrix);
	fclose(output);

	Matrix->free(Matrix);
	free(modularity);
	free(Matrix);
	return 0;
}
