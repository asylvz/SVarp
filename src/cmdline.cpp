#include <iostream>
#include <getopt.h>
#include "cmdline.h"


int parse_command_line(int argc, char** argv, parameters* params)
{
	int index, o;
	
	static struct option long_options[] = 
	{
		{"gaf" , required_argument, NULL, 'a'},
		{"graph" , required_argument, NULL, 'g'},
		{"fastq" , required_argument, NULL, 'f'},
		{NULL, 0, NULL, 0}
	};
	
	while((o = getopt_long( argc, argv, "a:g:f:", long_options, &index)) != -1)
	{
		switch(o)
		{
			case 'a':
				params->gaf = optarg;
				break;	
			case 'g':
				params->ref_graph = optarg;
				break;
			case 'f':
				params->fastq = optarg;
				break;
		}
	}	
	/* check if --ref   is invoked */
	if((params->ref_graph).empty())
	{
		fprintf( stderr, "[PSVPAN CMDLINE ERROR] Please enter reference graph genome file (GFA) using the --graph (-g) option.\n");
		return RETURN_ERROR;
	}

	if((params->gaf).empty())
	{
		fprintf( stderr, "[PSVPAN CMDLINE ERROR] Please enter alignment file (GAF) using the --gaf (-a) option.\n");
		return RETURN_ERROR;
	}	

	return RETURN_SUCCESS;
}

