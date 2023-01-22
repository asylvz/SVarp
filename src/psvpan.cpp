#include <iostream>
#include <set>
#include <map>
#include <unistd.h>
#include <filesystem>
#include "cmdline.h"
#include "gfa.h"
#include "alignment.h"
#include "assembly.h"
#include "phasing.h"

//ctags *.c

int main(int argc, char** argv)
{

	std::map <std::string, variant*> insertions;	
	//std::multimap<std::string, alignment*> alignments;
	std::set<std::string> contigs;	
	std::map<std::string, gfaNode*> ref;
	std::map <std::string, phase*> phased_reads;

	parameters* params = new parameters;	
	if (parse_command_line(argc, argv, params) != RETURN_SUCCESS)
		return RETURN_ERROR;


	std::string cwd = std::filesystem::current_path().string();
	std::string log_path = cwd + "/log/";	
	std::cout<<"\nLogs will be written to: "<<log_path<<std::endl;
	if(std::filesystem::exists(log_path))
		std::filesystem::remove_all(log_path);

	std::cout<<"\nInput files are:\n\t"<<params->gaf<<"\n\t"<< params->ref_graph<<"\n\t"<<params->fasta<<std::endl;
	
	ref = read_gfa(params, contigs);
	if (read_alignments(params, ref, insertions) != RETURN_SUCCESS)
		std::cout<<"Alignments could not be read\n";

	//find_supporting_reads(ref, alignments, contigs, insertions);
	//Read the TSV file and phase the reads
	if (read_phase_file(params, phased_reads) == RETURN_SUCCESS)
		phase_svs(phased_reads, insertions);

	run_assembly(params, insertions);	

	char *username= new char[50];
	getlogin_r( username, 50);
	
	std::cout<<"\nThank you for using pSVpan "<<username<< "... Hope to see you again.\n" <<std::endl;
	
	delete[] username;
	
	return RETURN_SUCCESS;
}


void check_size(parameters* params, std::map<std::string, variant*> insertions, std::map<std::string, unsigned long> fasta_index, std::map <std::string, phase*> phased_reads)
{
	std::map<std::string, phase*>::iterator itr;
	for (itr=phased_reads.begin(); itr != phased_reads.end(); ++itr)
	{
		
		if(itr->second->haplotype == "none" || itr->second->phase_set == "none")
		{
				continue;
		}

	
	}
}

