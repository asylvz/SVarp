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

#include <fstream>
//ctags *.c


void check_size(parameters* params, std::map <std::string, phase*> phased_reads)
{
	std::cout<<"Check read size\n";
	std::map<std::string, unsigned long> fasta_index;
	index_fasta(params, fasta_index);	
	std::map<std::string, phase*>::iterator itr;
	
	std::ifstream fp_read(params->fasta);
	
	unsigned long none_total = 0, total = 0;
	int total_line = 0, none_line = 0;
	for (itr=phased_reads.begin(); itr != phased_reads.end(); ++itr)
	{	
		unsigned long read_size = 0;	
		unsigned long char_pos;	
		std::string line;	
		string read = itr->first;
		
		if (fasta_index.find(read)!=fasta_index.end())
		{
			char_pos = fasta_index[read];

			fp_read.seekg(char_pos, std::ios::beg);
			
			getline(fp_read, line);
			//cout<<read<<" - "<<line<<std::endl;
			line.clear();
			while(getline(fp_read, line))
			{
				read_size += line.length();
				if(!fp_read || line[0] == '>')
					break;
				line.clear();
			}
		}
		else
			std::cout<<"Not found"<<read<<std::endl;	
	

		if(itr->second->haplotype == "none" || itr->second->phase_set == "none")
		{
			none_line++;
			none_total += read_size;
		}
		else
		{
			total_line++;
			total += read_size;
		}
		
		//cout<<none_line<<" - "<<total_line<<std::endl;
	}
	std::cout<<"None = "<<(double) none_total / none_line <<"\nPhased = "<<total/total_line<<std::endl;
}


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

	//Read the TSV file and phase the reads
	if (read_phase_file(params, phased_reads) == RETURN_SUCCESS)
		phase_svs(phased_reads, insertions);
	
	//check_size(params, phased_reads);
	run_assembly(params, insertions);	

	char *username= new char[50];
	getlogin_r( username, 50);
	
	std::cout<<"\nThank you for using pSVpan "<<username<< "... Hope to see you again.\n" <<std::endl;
	
	delete[] username;
	
	return RETURN_SUCCESS;
}



