#include <iostream>
#include <set>
#include <map>
#include "cmdline.h"
#include "gfa.h"
#include "alignment.h"

//ctags *.c

int main(int argc, char** argv)
{
	std::multimap <std::string, variant*> insertions;	
	std::set<std::string> contigs;	
	parameters* params;
	init_params(&params);
	
	int return_value = parse_command_line(argc, argv, params);			
	printf("\nThe input files are: \n\t%s\n\t%s\n", params->gaf, params->ref_graph);
	
	std::map<std::string, gfaNode*> ref = read_gfa(params, contigs);
	std::multimap<std::string, alignment*> alignments = read_alignments(params, ref, insertions);	
	find_supporting_reads(ref, alignments, contigs, insertions);
		
	return 0;
}
