#include <iostream>
#include <set>
#include <map>
#include "cmdline.h"
#include "gfa.h"
#include "alignment.h"
#include "assembly.h"

//ctags *.c

int main(int argc, char** argv)
{
	std::multimap <std::string, variant*> insertions;	
	std::multimap<std::string, alignment*> alignments;
	std::set<std::string> contigs;	
	std::map<std::string, gfaNode*> ref;

	parameters* params = new parameters;	
	if (parse_command_line(argc, argv, params) != RETURN_SUCCESS)
		return RETURN_ERROR;

	
	std::cout<<"\nThe input files are:\n\t"<<params->gaf<<"\n\t"<< params->ref_graph<<"\n\t"<<params->fasta<<std::endl;
	
	ref = read_gfa(params, contigs);
	alignments = read_alignments(params, ref, insertions);	

	find_supporting_reads(ref, alignments, contigs, insertions);
	
	run_assembly(params, insertions);	

	return RETURN_SUCCESS;
}


