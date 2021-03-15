#include <math.h>
#include <stdlib.h>
#include "spmat.h"
#include "algorithms.h"
#include "pIteration.h"
#include "eHandler.h"
#define ISPOSITIVE(X) ((X) > 0.00001)

/*Calculate modularity give a matrix and partition s.*/
double modularityCalculation(spmat *matrix, int *s, int *indices, int size) {
	int i;
	double value = 0;
	double *subVector1, *subVector2, *subVector3, *sumVector;
	subVector1 = (double*) calloc(size, sizeof(double));
	allocCheck(subVector1, __FILE__, __LINE__);
	subVector2 = (double*) calloc(size, sizeof(double));
	allocCheck(subVector2, __FILE__, __LINE__);
	subVector3 = (double*) calloc(size, sizeof(double));
	allocCheck(subVector3, __FILE__, __LINE__);
	sumVector = (double*) calloc(size, sizeof(double));
	allocCheck(sumVector, __FILE__, __LINE__);
	mult1SVector(matrix, indices, size, s, subVector1);
	mult2SVector(matrix, indices, size, s, subVector2);
	mult3SVector(matrix,indices, size, s, subVector3);
	addVectors(size, subVector1, subVector2, subVector3, sumVector);

	for (i = 0; i < size; i++) {
		value += (double) s[i] * sumVector[i];
	}

	free(subVector1);
	free(subVector2);
	free(subVector3);
	free(sumVector);

	return value/2;
}

/*Remove element from list.*/
void remove_element(int *array, int index, int array_length) {
	int i;
	for (i = index; i < array_length - 1; i++)
		array[i] = array[i + 1];
}

/*Removing group from group list..*/
void remove_group(group **group, int index, int array_length, int maxSize) {
	int i = index;
	int *temp;
	if (array_length == maxSize) {
		group[i]->indices = NULL;
		group[i]->size = 0;
	} else {
		temp = (group[index])->indices;
		for (; i < array_length; i++) {
			group[i]->indices = group[i + 1]->indices;
			group[i]->size = group[i + 1]->size;
		}
		free(temp);
	}
}

void buildMatrixRowVector(spmat *matrix, int row, int size, int *indices, double *result) {

	int i;
	double aValue, kValue;
	for (i = 0; i < size; i++) {
		aValue = matrix->actualValue(matrix, indices[row], indices[i]);
		kValue = (matrix->KList)[indices[i]]
				* (matrix->KFactorList)[indices[row]];
		result[i] = aValue - kValue;
	}
}

/*Find the best modularity change with movement of 1 vertex. Algorithm 4 from the forum.*/
void findModularityPartitionChange(spmat *matrix, int *s, int size,
		int *indices, double *improve, int *unmoved,
		int *improvementIndices) {
	int i, k, j, indexToRemove,first;
	double max; double *score;
	double *rowMatrix;

	for (i = 0; i < size; i++) {
		score = (double*) calloc(size, sizeof(double));
		allocCheck(score, __FILE__, __LINE__);
		for (k = 0; k < size; k++) {
			if(!unmoved[k]){
				double q0 = ((matrix->KFactorList)[indices[k]]
								* (matrix->KList)[indices[k]]);
				s[k] = -s[k];
				rowMatrix = (double*) calloc(size, sizeof(double));
				buildMatrixRowVector(matrix, k, size, indices,rowMatrix);
				score[k] = 4 * ((s[k] * dotProductInt(s, rowMatrix, size)) + q0);
				s[k] = -s[k];
				free(rowMatrix);
			}
		}

		first = 1;
		for (j = 0; j < size; j++) {
			if(!unmoved[j]){
				if(first){
					max = score[j];
					indexToRemove = j;
					first=0;
				}else{
					if (score[j] > max) {
						max = score[j];
						indexToRemove = j;
					}
				}
			}
		}

		s[indexToRemove] = -s[indexToRemove];
		improvementIndices[i] = indexToRemove;
		if (i == 0)
			improve[i] = max;
		else
			improve[i] = improve[i - 1] + max;
		unmoved[indexToRemove]=1;

		free(score);
	}
}

/*Improving initial division. Algorithm 4 from the forum.*/
double improvingDivision(spmat *Matrix, int *s, int *indices, int size) {
	int i, j, maxImproveIndex,counter=0;
	double maxImprove;
	int *unmoved, *improvementIndices;
	double *improve;
	double deltaQ;
	deltaQ = (double) 1;

	while (ISPOSITIVE(deltaQ)) {
		if(counter >  1000 * size){
			infiniteLoop();
		}
		unmoved = (int*) calloc(size, sizeof(int));
		allocCheck(unmoved, __FILE__, __LINE__);
		improve = (double*) calloc(size, sizeof(double));
		allocCheck(improve, __FILE__, __LINE__);
		improvementIndices = (int*) calloc(size, sizeof(int));
		allocCheck(improvementIndices, __FILE__, __LINE__);

		for (i = 0; i < size; i++) {
			unmoved[i] = 0;
		}

		/*Create improve list which detail how each vertex move change the modularity.*/
		findModularityPartitionChange(Matrix, s, size, indices, improve,
				unmoved, improvementIndices);

		maxImprove = improve[0];
		maxImproveIndex = 0;

		for (j = 0; j < size; j++) {
			if (improve[j] > maxImprove) {
				maxImprove = improve[j];
				maxImproveIndex = j;
			}
		}

		for (i = size - 1; i > maxImproveIndex; i--) {
			j = improvementIndices[i];
			s[j] = -s[j];
		}

		if (maxImproveIndex == size - 1) {
			deltaQ = 0;
		} else {
			deltaQ = maxImprove;
		}

		counter+=1;
	}

	/*Free all.*/

	free(unmoved);
	free(improve);
	free(improvementIndices);
	return deltaQ;
}

/*Generate s vector according to leading eigenvector.*/
void changeSVector(int *s, double *v, int size) {
	int i;
	for (i = 0; i < size; i++) {
		if (v[i] > 0) {
			s[i] = 1;
		} else {
			s[i] = -1;
		}
	}
}

/*Algorithm two.*/
void algorithmTwo(spmat *matrix, int *indices, int size, int *s,
		double *modularity) {
	double *leadingEigenValue, *leadingEigenVector;
	double mod;
	leadingEigenValue = malloc(sizeof(double));
	allocCheck(leadingEigenValue, __FILE__, __LINE__);
	leadingEigenVector = (double*) calloc(size, sizeof(double));
	allocCheck(leadingEigenVector, __FILE__, __LINE__);
	leadEigenPair(matrix, indices, size, leadingEigenVector, leadingEigenValue);
	if ((*leadingEigenValue) <= 0) {
		free(leadingEigenVector);
		free(leadingEigenValue);
		return;
	}

	changeSVector(s, leadingEigenVector, size);
	mod = modularityCalculation(matrix, s, indices, size);
	mod += improvingDivision(matrix, s, indices, size);
	if (mod < 0) {
		free(leadingEigenVector);
		free(leadingEigenValue);
		return;
	}

	*modularity += mod;
	free(leadingEigenVector);
	free(leadingEigenValue);
}

/*Adding the indices group to the relevant spot in the array.*/
void addIndicesToGroup(int *indices, group **groupList, int *groupSize,
		int size) {
	int i;
	(groupList[*groupSize])->indices = (int*)calloc(size,sizeof(int));
	allocCheck((groupList[*groupSize])->indices, __FILE__, __LINE__);
	for (i = 0; i < size; i++) {
		((groupList[*groupSize])->indices)[i] = indices[i];
	}
	(groupList[*groupSize])->size = size;
}

/*Sorting the division of the 2 groups, and inserting each one of them to the relevant
 group according to algorithm 3.*/
void sortingDividedGroups(int *g1, int *g2, int *sizeG1, int *sizeG2, group **P,
		group **O, int *sizeP, int *sizeO) {
	if (*sizeG1 > 1 && *sizeG2 > 1) {
		addIndicesToGroup(g1, P, sizeP, *sizeG1);
		(*sizeP) += 1;
		addIndicesToGroup(g2, P, sizeP, *sizeG2);
		(*sizeP) += 1;
	} else if (*sizeG1 == 1 && *sizeG2 == 1) {
		addIndicesToGroup(g2, O, sizeO, 1);
		(*sizeO) += 1;
		addIndicesToGroup(g2, O, sizeO, 1);
		(*sizeO) += 1;
	} else {
		if (*sizeG1 == 1) {
			addIndicesToGroup(g1, O, sizeO, 1);
			addIndicesToGroup(g2, P, sizeP, *sizeG2);
		} else {
			addIndicesToGroup(g2, O, sizeO, 1);
			addIndicesToGroup(g1, P, sizeP, *sizeG1);
		}
		(*sizeO) += 1;
		(*sizeP) +=1;
	}
}

/*Checking the division of the s vector to see if algorithm 2 split our group to just 1 group or 2.
 Returns true if the split is only to 1 group, false if the split is to 2 different group*/
int checkDivision(int *s, int *indices, int *g1, int *g2, int size,
		int *sizeFirst, int *sizeSecond) {
	int i;
	*sizeFirst = 0;
	*sizeSecond = 0;
	for (i = 0; i < size; i++) {
		if (s[i] == 1) {
			g1[*sizeFirst] = indices[i];
			(*sizeFirst)++;
		} else {
			g2[*sizeSecond] = indices[i];
			(*sizeSecond)++;
		}
	}

	if (*sizeFirst == 0 || *sizeFirst == size)
		return 1;

	return 0;
}

/*Algorithm three*/
void algorithmThree(spmat *matrix, group **P, group **O, int *sizeO,double* modularity) {
	int *sizeP, *g1, *g2, *indices, *s, *sizeG1, *sizeG2;
	int i, pGroupSize,counter=0;
	sizeP = malloc(sizeof(int));
	allocCheck(sizeP, __FILE__, __LINE__);
	*sizeP = 1;
	*modularity = (double) 0;

	while ((*sizeP) > 0) {
		if(counter > 1000 * matrix->n){
			infiniteLoop();
		}
		pGroupSize = P[0]->size;

		g1 = (int*) calloc(pGroupSize, sizeof(int));
		allocCheck(g1, __FILE__, __LINE__);
		g2 = (int*) calloc(pGroupSize, sizeof(int));
		allocCheck(g2, __FILE__, __LINE__);
		indices = (int*) calloc(pGroupSize, sizeof(int));
		allocCheck(indices, __FILE__, __LINE__);
		s = (int*) calloc(pGroupSize, sizeof(int));
		allocCheck(s, __FILE__, __LINE__);
		for (i = 0; i < pGroupSize; i++) {
			indices[i] = (P[0]->indices)[i];
			s[i] = 1;
		}

		remove_group(P, 0, *sizeP, matrix->n);
		(*sizeP) -=1;

		sizeG1 = malloc(sizeof(int));
		allocCheck(sizeG1, __FILE__, __LINE__);
		sizeG2 = malloc(sizeof(int));
		allocCheck(sizeG2, __FILE__, __LINE__);
		algorithmTwo(matrix, indices, pGroupSize, s, modularity);
		if (checkDivision(s, indices, g1, g2, pGroupSize, sizeG1, sizeG2)) {
			addIndicesToGroup(indices, O, sizeO, pGroupSize);
			(*sizeO) +=1;
		} else {
			sortingDividedGroups(g1, g2, sizeG1, sizeG2, P, O, sizeP, sizeO);
		}

		counter += 1;

		free(indices);
		free(s);
		free(sizeG1);
		free(sizeG2);
		free(g1);
		free(g2);
	}

	free(sizeP);
}
