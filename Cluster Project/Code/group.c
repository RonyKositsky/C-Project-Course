#include "group.h"
#include "eHandler.h"
#include <stdlib.h>

/*Allocating group, depending on it's size.*/
void allocateGroupList(group** list, int size){
	int i;
	for (i = 0; i < size; i++) {
		group* instance = (group*)calloc(1,sizeof(group));
		allocCheck(instance, __FILE__, __LINE__);
		list[i] = instance;
	}
}


