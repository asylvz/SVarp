#include <iostream>
#include <set>
#include <map>
#include <unistd.h>
#include "cmdline.h"
#include "reference.h"
#include "alignment.h"
#include "assembly.h"
#include "phasing.h"

//ctags *.c


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

	init_logs(params);	
	
	read_gfa(params, ref, gfa);
	if (read_alignments(params, ref, gfa, insertions_tmp) != RETURN_SUCCESS)
		std::cout<<"Alignments could not be read\n";
				
	refine_svs(insertions_tmp, insertions);

	//Read the TSV file and phase the reads
	std::cout<<"Phasing"<<std::endl;	
	if (read_phase_file(params, phased_reads) == RETURN_SUCCESS)
		phase_svs(phased_reads, insertions);
	
	//check_size(params, phased_reads);
	run_assembly(params, insertions);	

	char *username= new char[50];
	getlogin_r(username, 50);
	
	std::cout<<"\nThank you for using svapan "<<username<< "... Hope to see you again.\n" <<std::endl;
	
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
