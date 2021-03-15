#include<stdlib.h>
#include <assert.h>
#include <time.h>
#include <stdio.h>

FILE* input_file;
FILE* output_file;
double** Matrix;
int rows = 0;
int colums = 0;
void CopyValues(FILE* input_file);

double mean(double* arr, int n)
{
	int i;
	double sum = 0;
	for (i = 0; i < n; i++)
		sum = sum + arr[i];
	return sum / n;
}

void StandrizeMatrix() {
	int i, j;
	for (i = 0; i < rows; i++) {
		double rowMean = mean(Matrix[i], colums);
		for (j = 0; j < colums; j++) {
			Matrix[i][j] = Matrix[i][j] - rowMean;
		}
	}
}

void GenerateMatrix(FILE* file,int rows) {
	int i;
	Matrix = (double**)calloc(rows, sizeof(double*));
	for (i = 0; i < rows; i++) {
		Matrix[i] = (double*)calloc(rows, sizeof(double));
	}

	CopyValues(file);
	fclose(file);

	StandrizeMatrix();
}

void CopyValues(FILE* input_file) {
	int i;

	for (i = 0; i < rows; i++) {
		fread(Matrix[i], sizeof(double), colums, input_file);
	}
}

double covariance(double *arr1, double *arr2, int n)
{
	int i; double sum = 0;
	for (i = 0; i < n; i++)
		sum = sum + (arr1[i]) *(arr2[i]);
	return sum;
}

void WriteToFile(char* fileName) {
	int i, j;
	double* CovRoW;
	
	output_file = fopen(fileName, "w");
	rewind(output_file);

	fwrite(&rows, sizeof(int), 1, output_file);
	fwrite(&rows, sizeof(int), 1, output_file);

	CovRoW = calloc(sizeof(int), rows);

	for (i = 0; i < rows; i++) {
		for (j = 0; j < rows; j++) {
			CovRoW[j] = covariance(Matrix[i], Matrix[j], colums);
		}

		fwrite(CovRoW, sizeof(double), rows, output_file);
	}

	fclose(output_file);
}

int main(int argc, char* argv[]) {

	clock_t startTime, finishTime;
	double time;
	assert(argc != -1);
	startTime = clock();
	input_file= fopen(argv[1], "r");

	fread(&colums, sizeof(int), 1, input_file);
	fread(&rows, sizeof(int), 1, input_file);

	GenerateMatrix(input_file,rows);
	WriteToFile(argv[2]);
	finishTime = clock();
	time = ((double)(finishTime - startTime) / CLOCKS_PER_SEC);
	printf("\n---running time: %f seconds ---\n", time);
	return 0;
}
