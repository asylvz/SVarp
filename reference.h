#ifndef __REFERENCE
#define __REFERENCE
#include "common.h"


typedef struct _refs
{
	char* name; /* name of the node */
	char* contig; /* name of the contig that the node belongs to, e.g., SN:Z:chr1 */
	char* sequence; /* sequence of the node - second column of the GFA */
	int len; /*node length, e.g., LN:i:74 */
	int offset; /*node offset within the contig e.g., SO:i:10746 */
} reference;

typedef struct _gfa { /* table entry: */ 		
    reference *node; /* replacement text */
	struct _gfa *next; /* next entry in chain */
} gfa;


gfa *gfa_insert(gfa **hashtab, char *name, reference *node);
gfa *gfa_lookup(gfa **hashtab, char *s);
void gfa_traverse(gfa** hashtab);

int free_gfa(gfa** hash_table);
int read_reference_graph(parameters* params, gfa**);
void init_reference(reference** ref);

#endif
