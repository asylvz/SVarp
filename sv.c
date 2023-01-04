#include <stdio.h>
#include <stdlib.h>
#include "sv.h"
#include "reference.h"

sv* generate_sv_node(gfa* gfa_table, int ref_pos, char strand, char* path, int path_start, int path_end, int sv_size)
{
	sv* svar = (sv*) getMem(sizeof(variant));
	svar->var->sv_size = sv_size;
}
