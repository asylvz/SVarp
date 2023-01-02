#include <stdio.h>
#include "common.h"
#include "cmdline.h"
#include "alignment.h"
#include "reference.h"


int main(int argc, char** argv)
{
	printf("Loading...\n");
	parameters* params;
	
	gfa *gfa_table[HASHSIZE];
	for(int i = 0; i < HASHSIZE; i++)
		gfa_table[i] = NULL;
	init_params(&params);
	
	int return_value = parse_command_line( argc, argv, params);			
	printf("The input files are: %s\n%s\n", params->gaf, params->ref_graph);
	
	// Read the S lines in the GFA file
	read_reference_graph(params, gfa_table);
	gfa_traverse(gfa_table);
	
	
	
	
	free_gfa(gfa_table);
	return 0;
}
