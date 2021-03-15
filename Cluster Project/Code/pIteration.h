#ifndef PITERATION_H_
#define PITERATION_H_

#include <stdio.h>
#include "spmat.h"

void leadEigenPair(spmat *A, int *arr, int size, double *destVec,
		double *destValue);

void mult1SVector(spmat *matrix, int *arr, int arrSize, int *sVector,
		double *subVector1);

void mult2SVector(spmat *matrix, int *arr, int arrSize, int *sVector,
		double *subVector2);

void mult3SVector(spmat *matrix,int* indices, int arrSize, int *sVector, double *subVector3);

void addVectors(int arrSize, double *subMTV1, double *subMTV2, double *subMTV3,
		double *MTV);

double dotProductInt(int *vec1, double *vec2, int size);

#endif
