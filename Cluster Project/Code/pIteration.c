#include <math.h>
#include "spmat.h"
#include <stdio.h>
#include <stdlib.h>
#include "pIteration.h"
#include "eHandler.h"
/********* Section 3 - Leading Eigenpair ********/

/*  Performing power iterations to obtain vector bk*/

#define epsilon 0.00001

void generateRandomVector(int size, double *vector) {
	int i;
	for (i = 0; i < size; i++) {
		vector[i] = rand();
	}
}

void copyVectors(double *from, double *to, int size) {
	int i;
	for (i = 0; i < size; i++) {
		to[i] = from[i];
	}
}

int vectorDifference(double *newVector, double *oldVector, int size) {
	int i;
	for (i = 0; i < size; i++) {
		if (fabs(newVector[i] - oldVector[i]) >= epsilon) {
			return 1;
		}
	}
	return 0;
}

void normalize(double *secondVector, int size) {
	double norma;
	double sum = 0;
	int i;
	for (i = 0; i < size; i++) {
		sum += pow(secondVector[i], 2);
	}
	norma = sqrt(sum);
	for (i = 0; i < size; i++) {
		secondVector[i] = (secondVector[i] / norma);
	}
}

/*multiplying submatrix (according to arr) of cTag by firstVector. Result in subMTV1 */
void mult1(spmat *cTag, int *arr, int arrSize, double *firstVector,
		double *subMTV1) {
	int i, j;
	double result = 0;
	for (i = 0; i < arrSize; i++) {
		for (j = 0; j < arrSize; j++) {
			result +=  (cTag->actualValue(cTag, arr[i], arr[j]))*firstVector[j];
		}
		subMTV1[i] = result;
		result = 0;
	}
}

/*multiplying submatrix (according to arr) of cTag by firstVector. Result in subMTV1 */
void mult1SVector(spmat *matrix, int *arr, int arrSize, int *sVector,
		double *subVector1) {
	int i, j;
	double result = 0;
	for (i = 0; i < arrSize; i++) {
		for (j = 0; j < arrSize; j++) {
			result += matrix->actualValue(matrix, arr[i], arr[j])
					* (double) sVector[j];
		}
		subVector1[i] = result;
		result = 0;
	}
}

void mult2(spmat *cTag, int *arr, int arrSize, double *firstVector,
		double *subMTV2) {
	int i, j;
	double result = 0;
	for (i = 0; i < arrSize; i++) {
		for (j = 0; j < arrSize; j++) {
			result += cTag->KFactorList[arr[i]] * cTag->KList[arr[j]]
					* firstVector[j];
			if (i == j) {
				result += (cTag->sumOfRowB(cTag, arr, arrSize, i))
						* firstVector[i];
			}
		}
		subMTV2[i] = result;
		result = 0;
	}
}

/*multiplying the K normalized submatrix (according to arr) of cTag by firstVector. Result in subMTV2*/
void mult2SVector(spmat *matrix, int *arr, int arrSize, int *sVector,
		double *subVector2) {
	int i, j;
	double result = 0;
	for (i = 0; i < arrSize; i++) {
		for (j = 0; j < arrSize; j++) {
			result += matrix->KFactorList[arr[i]] * matrix->KList[arr[j]]
					* (double) sVector[j];
		}
		subVector2[i] = result;
		result = 0;
	}
}

/*multiplying firstVector by factor of cTag's 1-norm. Result in subMTV3*/
void mult3(spmat *cTag, int arrSize, double *firstVector, double *subMTV3) {
	int i;
	double factor = cTag->norma;
	for (i = 0; i < arrSize; i++) {
		subMTV3[i] = firstVector[i] * factor;
	}
}

/*multiplying firstVector by factor of cTag's 1-norm. Result in subMTV3*/
void mult3SVector(spmat *matrix,int* indices, int arrSize, int *sVector, double *subVector3) {
	int i;
	for (i = 0; i < arrSize; i++) {
		subVector3[i] = -(double) sVector[i] * (matrix->sumOfRowB(matrix,indices,arrSize,i));
	}
}

/*linear combination of the 3 subVectors creating the multiplication result of sub-matrix by firstVector
 * i.e completing one "round" of power iteration (before normalizing)*/

void addVectors(int arrSize, double *subMTV1, double *subMTV2, double *subMTV3,
		double *MTV) {
	int i;
	for (i = 0; i < arrSize; i++) {
		MTV[i] = subMTV1[i] - subMTV2[i] + subMTV3[i];
	}
}

void matrixTimesVector(double *firstVector, double *MTV, double *secondVector,
		spmat *cTag, int *arr, int arrSize) {
	double *subMTV1, *subMTV2, *subMTV3;
	subMTV1 = (double*) calloc(arrSize, sizeof(double));
	allocCheck(subMTV1, __FILE__, __LINE__);
	subMTV2 = (double*) calloc(arrSize, sizeof(double));
	allocCheck(subMTV2, __FILE__, __LINE__);
	subMTV3 = (double*) calloc(arrSize, sizeof(double));
	allocCheck(subMTV3, __FILE__, __LINE__);
	mult1(cTag, arr, arrSize, firstVector, subMTV1);/*first part - A times vector v, put result in subMTV1*/
	mult2(cTag, arr, arrSize, firstVector, subMTV2);/*second part -subtraction of K times v, put result in subMTV2*/
	mult3(cTag, arrSize, firstVector, subMTV3);/*third part - v times the 1-norm, put result in subMTV3*/
	addVectors(arrSize, subMTV1, subMTV2, subMTV3, MTV);
	copyVectors(MTV, secondVector, arrSize);
	free(subMTV1);
	free(subMTV2);
	free(subMTV3);
}

void nextInSeries(double *firstVector, double *MTV, double *secondVector,
		spmat *cTag, int *arr, int arrSize) {
	matrixTimesVector(firstVector, MTV, secondVector, cTag, arr, arrSize);
	normalize(secondVector, arrSize);
}

/*performing power iteration on sub-matrix (according to arr) of cTag, Result in secondVector */
void powerIteration(double *firstVector, double *MTV, double *secondVector,
		spmat *cTag, int *arr, int arrSize) {
	int counter = 0;
	nextInSeries(firstVector, MTV, secondVector, cTag, arr, arrSize);
	while (vectorDifference(secondVector, firstVector, arrSize)) {
		if(counter > 1000* arrSize){
			infiniteLoop();
		}
		copyVectors(secondVector, firstVector, arrSize);
		nextInSeries(firstVector, MTV, secondVector, cTag, arr, arrSize);
		counter+=1;
	}

}

double dotProduct(double *vec1, double *vec2, int size) {
	int i;
	double sum = 0;
	for (i = 0; i < size; i++) {
		sum += vec1[i] * vec2[i];
	}
	return sum;
}

double dotProductInt(int *vec1, double *vec2, int size) {
	int i;
	double sum = 0;
	for (i = 0; i < size; i++) {
		sum += ((double) vec1[i]) * vec2[i];
	}
	return sum;
}

/*finding leading eigenpair of sub-matrix according to arr, putting the eigenvector in vector, and the eigenvalue in value*/
void leadEigenPair(spmat *matrix, int *arr, int arrSize, double *vector,
		double *value) {
	double norma = matrix->norma;
	double eigenValue;
	double *firstVector, *MTV, *lastVector1, *lastVector2, *lastVector3, *last,
			*secondVector;

	/* *** Alocations *** */
	firstVector = (double*) calloc(arrSize, sizeof(double));
	allocCheck(firstVector, __FILE__, __LINE__);
	MTV = (double*) calloc(arrSize, sizeof(double));
	allocCheck(MTV, __FILE__, __LINE__);
	lastVector1 = (double*) calloc(arrSize, sizeof(double));
	allocCheck(lastVector1, __FILE__, __LINE__);
	lastVector2 = (double*) calloc(arrSize, sizeof(double));
	allocCheck(lastVector2, __FILE__, __LINE__);
	lastVector3 = (double*) calloc(arrSize, sizeof(double));
	allocCheck(lastVector3, __FILE__, __LINE__);
	last = (double*) calloc(arrSize, sizeof(double));
	allocCheck(last, __FILE__, __LINE__);
	secondVector = (double*) calloc(arrSize, sizeof(double));
	allocCheck(secondVector, __FILE__, __LINE__);

	/* *** PI & leading eigen pair computation *** */
	generateRandomVector(arrSize, firstVector);
	powerIteration(firstVector, MTV, secondVector, matrix, arr, arrSize);
	mult1(matrix, arr, arrSize, secondVector, lastVector1);
	mult2(matrix, arr, arrSize, secondVector, lastVector2);
	mult3(matrix, arrSize, secondVector, lastVector3);
	addVectors(arrSize, lastVector1, lastVector2, lastVector3, last);
	eigenValue = (dotProduct(secondVector, last, arrSize)
			/ dotProduct(secondVector, secondVector, arrSize));
	eigenValue -= norma;
	copyVectors(secondVector, vector, arrSize);
	*value = eigenValue;
	free(firstVector);
	free(secondVector);
	free(MTV);
	free(lastVector1);
	free(lastVector2);
	free(lastVector3);
	free(last);

}
