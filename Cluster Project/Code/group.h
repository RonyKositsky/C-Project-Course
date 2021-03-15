typedef struct _group {
	int *indices;
	int size;
} group;

void allocateGroupList(group** list, int size);

