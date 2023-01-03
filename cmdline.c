#include <stdio.h>
#include <getopt.h>
#include "cmdline.h"


int parse_command_line(int argc, char** argv, parameters* params)
{
	int index, o;
	
	static struct option long_options[] = 
	{
		{"gaf" , required_argument, NULL, 'i'},
		{"graph" , required_argument, NULL, 'g'},
		{NULL, 0, NULL, 0}
	};
	
	while((o = getopt_long( argc, argv, "i:g:", long_options, &index)) != -1)
	{
		switch(o)
		{
			case 'i':
				set_str(&(params->gaf), optarg);
				break;	
			case 'g':
				set_str(&(params->ref_graph), optarg);
				break;
		}
	}	
	/* check if --ref   is invoked */
	if( params->ref_graph == NULL)
	{
		fprintf( stderr, "[PSVPAN CMDLINE ERROR] Please enter reference graph genome file (GFA) using the --graph option.\n");
		return -1;
	}

	if( params->gaf == NULL)
	{
		fprintf( stderr, "[PSVPAN CMDLINE ERROR] Please enter alignment file (GAF) using the --gaf option.\n");
		return -1;
	}	

	return 1;
}

