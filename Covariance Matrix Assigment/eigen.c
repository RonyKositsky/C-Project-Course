#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <stdio.h>
#include "SPBufferset.h"


int rows;
int colums;
double epsilon = 0.00001;
double **Matrix;
double *firstVector;
double *secondVector;
double *MTV;
FILE *input_file;
FILE *output_file;

void vectorPrint(double *vector) {
	int i;
	printf("\n");
	for (i = 0; i < rows; i++) {
		printf("\n---[i]: %f seconds ---\n", vector[i]);
	}
}
void generateRandomVector(int size) {
	int i;
	firstVector = (double*) calloc(size, sizeof(double));
	assert(firstVector !=NULL);
	for (i = 0; i < size; i++) {
		firstVector[i] = rand();
	}
}
void generateVectors(int size) {
	secondVector = (double*) calloc(size, sizeof(double));
	assert(secondVector !=NULL);
	MTV = (double*) calloc(size, sizeof(double));
	assert(MTV !=NULL);
}

int vectorDifference(double *newVector, double *oldVector) {
	int i;
	for (i = 0; i < rows; i++) {
		if (fabs(newVector[i] - oldVector[i]) >= epsilon) {
			return 1;
		}
	}
	return 0;
}

void normalize() {
	double size;
	double sum = 0;
	int i;
	for (i = 0; i < rows; i++) {
		sum += pow(secondVector[i], 2);
	}
	size = sqrt(sum);
	for (i = 0; i < rows; i++) {
		secondVector[i] = (secondVector[i] / size);
	}
}

void copyVectors(double *from, double *to) {
	int i;
	for (i = 0; i < rows; i++) {
		to[i] = from[i];
	}
}

void matrixTimesVector() {
	int i;
	int j;
	double sum;
	for (i = 0; i < rows; i++) {
		for (j = 0; j < rows; j++) {
			sum += Matrix[i][j] * firstVector[j];
		}
		MTV[i] = sum;
		sum = 0;
	}
	copyVectors(MTV, secondVector);
}

void nextInSeries() {
	matrixTimesVector();
	normalize();
}

void powerIteration() {
	nextInSeries();
	while (vectorDifference(secondVector, firstVector)) {
		copyVectors(secondVector, firstVector);
		nextInSeries();
	}

}

void CopyValues(FILE *input_file) {
	int i;

	for (i = 0; i < rows; i++) {
		fread(Matrix[i], sizeof(double), colums, input_file);
	}
}

void GenerateMatrix(FILE *file, int rows) {
	int i;
	Matrix = (double**) calloc(rows, sizeof(double*));
	for (i = 0; i < rows; i++) {
		Matrix[i] = (double*) calloc(rows, sizeof(double));
	}
	CopyValues(file);
	fclose(file);
}

void freeMatrix() {
	int i;
	for (i = 0; i < rows; i++) {
		free(Matrix[i]);
	}
}

void WriteToFile(char *fileName) {
	int col = 1;
	int row;
	int check;
	output_file = fopen(fileName, "w");
	rewind(output_file);
	col = fwrite(&col, sizeof(int), 1, output_file);
	row = fwrite(&rows, sizeof(int), 1, output_file);
	check = fwrite(secondVector, sizeof(double), rows, output_file);
	assert(row == 1);
	assert(col == 1);
	assert(check == rows);
	fclose(output_file);
	
}

int main(int argc, char *argv[]) {
	clock_t startTime, finishTime;
	double time1;
	assert(argc != 1);
	startTime = clock();
	input_file = fopen(argv[1], "r");
	srand(time(NULL));
	fread(&colums, sizeof(int), 1, input_file);
	fread(&rows, sizeof(int), 1, input_file);
	GenerateMatrix(input_file, rows);
	generateRandomVector(rows);
	generateVectors(rows);
	powerIteration();
	WriteToFile(argv[2]);
	free(firstVector);
	free(secondVector);
	free(MTV);
	freeMatrix();
	free(Matrix);
	finishTime = clock();
	time1 = ((double) (finishTime - startTime) / CLOCKS_PER_SEC);
	printf("\n---running time: %f seconds ---\n", time1);
	return 0;

}
