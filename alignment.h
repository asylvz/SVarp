#ifndef __ALIGNMENT
#define __ALIGNMENT

#include "common.h"
#include "reference.h"

typedef struct _alignment
{
	char* read_name;
	int start;
	int end;
	char strand;
	char* node;
	char* path;	
	struct _alignment *next;	
} alignment;

int read_alignments(parameters*, gfa**);

#endif
