#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eHandler.h"

void allocCheck(void *element, const char *fileName, int line) {
	if (element == NULL) {
		printf("\n Error allocating memory in file %s, line number %d",
				fileName, --line);
		exit(1);
	}
}

void checkM(int M) {
	if (M == 0) {
		printf("\n Error calculating M value of the matrix - division by zero");
		exit(2);
	}
}

void infiniteLoop(){
	printf("\n Infinite loop detected.");
	exit(2);
}
