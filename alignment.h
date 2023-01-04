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
} alignment;


typedef struct _gaf { /* table entry: */
	char* contig;
    alignment *node; /* replacement text */
    struct _gaf *next; /* next entry in chain */
} gaf;


int read_alignments(parameters *params, gfa* gfa_table[], gaf* gaf_table[]);
void gaf_traverse(gaf** hashtab);

#endif
