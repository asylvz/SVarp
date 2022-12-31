#include <stdio.h>
#include "common.h"
#include "cmdline.h"
#include "alignment.h"
#include "reference.h"


int main(int argc, char** argv)
{
	printf("Loading...\n");
	parameters* params;
	
	init_params(&params);
	
	int return_value = parse_command_line( argc, argv, params);			
	printf("The input files are: %s\n%s\n", params->gaf, params->ref_graph);
	
	// Read the S lines in the GFA file
	read_reference_graph(params);
	
	/*rdict* dnm;
	install("Arda", "birinci");
	install("Ferda", "ikinci");
	dnm = lookup("Aasds");
	printf("%s - %s\n", dnm->name, dnm->defn);	
	*/

	return 0;
}
