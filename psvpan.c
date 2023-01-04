#include <stdio.h>
#include "common.h"
#include "cmdline.h"
#include "reference.h"
#include "alignment.h"
#include "free.h"

//ctags *.c

int main(int argc, char** argv)
{
	parameters* params;

	gfa *gfa_table[HASHSIZE];
	for(int i = 0; i < HASHSIZE; i++)
		gfa_table[i] = NULL;
	
	gaf *gaf_table[HASHSIZE];
	for(int i = 0; i < HASHSIZE; i++)
		gaf_table[i] = NULL;
		
	init_params(&params);
	
	int return_value = parse_command_line(argc, argv, params);			
	printf("\nThe input files are: \n\t%s\n\t%s\n", params->gaf, params->ref_graph);
	
	printf("\nReading the reference graph (GFA)\n");
	// Read the S lines in the GFA file
	
	read_reference_graph(params, gfa_table);
	//gfa_traverse(gfa_table);
	
	printf("Reading the alignments (GAF)\n");
	read_alignments(params, gfa_table, gaf_table);	
	
	//printf("Traversing\n");
	//gaf_traverse(gaf_table);
	
	printf("Freeing GFA table\n");
	free_gfa(gfa_table);
	

	printf("Freeing GAF table\n");
	free_gaf(gaf_table);
	return 0;
}
