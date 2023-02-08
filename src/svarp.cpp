#include <iostream>
#include <map>
#include <unistd.h>
#include <algorithm>
#include "cmdline.h"
#include "reference.h"
#include "alignment.h"
#include "assembly.h"
#include "phasing.h"
#include "remap.h"

//ctags *.c


int main(int argc, char** argv)
{
	std::map <std::string, variant*> variations_tmp;	
	std::map <std::string, gfaNode*> gfa;
	std::map <std::string, phase*> phased_reads;
	std::map <std::string, std::vector<svtig*>> insertions, deletions;
	std::map <std::string, Contig*> ref;
	
	parameters* params = new parameters;	
	if (parse_command_line(argc, argv, params) != RETURN_SUCCESS)
		return RETURN_ERROR;
	
	init_logs(params);	

	if (read_gfa(params, ref, gfa) != RETURN_SUCCESS)
		return RETURN_ERROR;

	if (read_alignments(params, ref, gfa, variations_tmp) != RETURN_SUCCESS)
		return RETURN_ERROR;

	refine_svs(variations_tmp, insertions, deletions);

	//Read the TSV file and phase the reads
	std::cout<<"Phasing"<<std::endl;	
	if (read_phase_file(params, phased_reads) == RETURN_SUCCESS)
		phase_svs(phased_reads, insertions, deletions);
	
	find_deletions(params, deletions);

	run_assembly(params, ref, insertions);	

	//merge assembly output in a single fasta file
	merge_assemblies();		
	
	//Run minigraph for remapping insertion assemblies
	remap_assemblies(params);	
	
	//Read remapped svtigs
	read_remappings(gfa);

	char *username= new char[50];
	getlogin_r(username, 50);
	
	std::cout<<"\nThank you for using SVarp "<<username<< "... Tschüs, güle güle, adios, bye...\n" <<std::endl;
	
	logFile.close();
	//delete[] username;

	return RETURN_SUCCESS;
}
