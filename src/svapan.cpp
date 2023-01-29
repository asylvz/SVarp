#include <iostream>
#include <set>
#include <map>
#include <unistd.h>
#include <filesystem>
#include "cmdline.h"
#include "reference.h"
#include "alignment.h"
#include "assembly.h"
#include "phasing.h"

#include <fstream>
//ctags *.c

std::ofstream logFile;

int main(int argc, char** argv)
{
	std::map <std::string, variant*> insertions_tmp;	
	std::map <std::string, gfaNode*> gfa;
	std::map <std::string, phase*> phased_reads;
	std::map <std::string, std::vector<svtig*>> insertions;
	std::map <std::string, Contig*> ref;

	parameters* params = new parameters;	
	if (parse_command_line(argc, argv, params) != RETURN_SUCCESS)
		return RETURN_ERROR;

	std::string cwd = std::filesystem::current_path().string();
	std::string log_path = cwd + "/log/";	
	std::cout<<"\nLogs will be written to: "<<log_path<<std::endl;
	if(std::filesystem::exists(log_path))
		std::filesystem::remove_all(log_path);
	
  	if (!std::filesystem::create_directory(log_path))
	{
		std::cerr << "Error creating log folder" << std::endl;
		exit(-1);
  	}
  	if (!std::filesystem::create_directory(log_path + "in/"))
	{
		std::cerr << "Error creating log/in/" << std::endl;
		exit(-1);
  	}
  	if (!std::filesystem::create_directory(log_path + "out/"))
	{
		std::cerr << "Error creating log/out/" << std::endl;
		exit(-1);
  	}

	logFile.open(log_path + "psvpan.log");
	
	std::cout<<"\nInput files are:\n\t"<<params->gaf<<"\n\t"<< params->ref_graph<<"\n\t"<<params->fasta<<std::endl;
	
	read_gfa(params, ref, gfa);
	if (read_alignments(params, ref, gfa, insertions_tmp) != RETURN_SUCCESS)
		std::cout<<"Alignments could not be read\n";
			
	std::map<std::string, Contig*>::iterator it;
	for (it=ref.begin(); it != ref.end(); ++it)
		logFile<< it->first<<"---> LEN= "<<it->second->contig_length<<" - Mapped bases= "<< it->second->mapped_bases<<" - Cov= "<<it->second->coverage <<std::endl;	
	
	refine_svs(insertions_tmp, insertions);

	std::cout<<"Phasing"<<std::endl;	
	//Read the TSV file and phase the reads
	if (read_phase_file(params, phased_reads) == RETURN_SUCCESS)
		phase_svs(phased_reads, insertions);
	
	//check_size(params, phased_reads);
	run_assembly(params, insertions);	

	char *username= new char[50];
	getlogin_r( username, 50);
	
	std::cout<<"\nThank you for using pSVpan "<<username<< "... Hope to see you again.\n" <<std::endl;
	
	logFile.close();
	//delete[] username;

	return RETURN_SUCCESS;
}



/*void check_size(parameters* params, std::map <std::string, phase*> phased_reads)
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
}*/
