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


void calculate_n50(parameters* params, std::map <std::string, phase*> phased_reads);

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

	//if (read_gfa(params, ref, gfa) != RETURN_SUCCESS)
	//	return RETURN_ERROR;

	//if (read_alignments(params, ref, gfa, insertions_tmp) != RETURN_SUCCESS)
	//	return RETURN_ERROR;

	//refine_svs(insertions_tmp, insertions);

	//Read the TSV file and phase the reads
	//std::cout<<"Phasing"<<std::endl;	
	//if (read_phase_file(params, phased_reads) == RETURN_SUCCESS)
	//	phase_svs(phased_reads, insertions);
	
	//calculate_n50(params, phased_reads);
	//run_assembly(params, insertions);	

	//merge assembly output in a single fasta file
	//merge_assemblies();		
	
	//Run minigraph for remapping insertion assemblies
	//remap_assemblies(params);	
	
	//Read remapped svtigs
	read_remappings(gfa);

	char *username= new char[50];
	getlogin_r(username, 50);
	
	std::cout<<"\nThank you for using svapan "<<username<< "... Tschüs, güle güle, adios, bye...\n" <<std::endl;
	
	logFile.close();
	//delete[] username;

	return RETURN_SUCCESS;
}



void calculate_n50(parameters* params, std::map <std::string, phase*> phased_reads)
{
	std::cout<<"\nCalculate N50\n";
	std::map<std::string, unsigned long> fasta_index;
	index_fasta(params, fasta_index);	
	std::map<std::string, phase*>::iterator itr;
	
	std::vector <int> phased_sizes;
	std::vector <int> unphased_sizes;

	std::ifstream fp_read(params->fasta);
	
	long unphased_total = 0, phased_total = 0, all_total = 0, total_line = 0, none_line = 0;
	for (itr=phased_reads.begin(); itr != phased_reads.end(); ++itr)
	{	
			
		long read_size = 0;	
		long char_pos;	
		std::string line;	
		std::string read = itr->first;
		

		if (fasta_index.find(read)!=fasta_index.end())
		{
			char_pos = fasta_index[read];

			fp_read.seekg(char_pos, std::ios::beg);
			
			//std::cout<<read<<" "<<char_pos<<std::endl;	
			getline(fp_read, line);
			//cout<<read<<" - "<<line<<std::endl;
			line.clear();
			while(getline(fp_read, line))
			{
				//std::cout<<line<<std::endl;
				if(line[0] == '>')
					break;
			
				//std::cout<<line.length()<<std::endl;
				read_size += line.length();
				line.clear();
				//std::cout<<read_size<<std::endl;
			}
			//std::cout<<read<<" "<<read_size<<std::endl;
		}
		else
		{
			std::cout<<"Not found "<<read<<std::endl;	
			continue;
		}
	

		if(itr->second->haplotype == "none" || itr->second->phase_set == "none")
		{
			none_line++;
			unphased_total += read_size;
			unphased_sizes.push_back(read_size);
			all_total += read_size;
		}
		else
		{
			total_line++;
			phased_total += read_size;
			phased_sizes.push_back(read_size);
			all_total += read_size;
		}

		//std::cout<<phased_total <<" "<<unphased_total<<" "<<read_size<<" " <<all_total <<std::endl;	
		//cout<<none_line<<" - "<<total_line<<std::endl;
	}

	std::cout<<"Unphased = "<<(double) unphased_total / none_line <<"\nPhased = "<<phased_total/total_line<<std::endl;

	sort(phased_sizes.begin(), phased_sizes.end(), std::greater<>());
	sort(unphased_sizes.begin(), unphased_sizes.end(), std::greater<>());
	
	long tmp_total = 0;
	bool test1 = false, test2 = false;
	for(size_t i = 0; i< phased_sizes.size(); i++)
	{
		//std::cout<<"PHASED "<<phased_sizes[i]<<std::endl;
		tmp_total += phased_sizes[i];
		if (tmp_total >= (all_total / 2) && !test1)
		{
			std::cout<<"(all total) N50 of phased is "<<phased_sizes[i]<<std::endl;
			test1 = true;
		}

		if (tmp_total >= (phased_total / 2) && !test2)
		{
			std::cout<<"(phased total) N50 of phased is "<<phased_sizes[i]<<std::endl;
			test2 = true;
		}
	}

	tmp_total = 0;
	test1 = false, test2 = false;
	for(size_t i = 0; i< unphased_sizes.size(); i++)
	{
		//std::cout<<"UNPHASED "<<unphased_sizes[i]<<std::endl;
		tmp_total += unphased_sizes[i];
		if (tmp_total >= (all_total / 2) and !test1)
		{
			std::cout<<"(all total) N50 of unphased is "<<unphased_sizes[i]<<std::endl;
			test1 = true;
		}

		if (tmp_total >= (unphased_total / 2) && !test2)
		{
			std::cout<<"(unphased total) N50 of unphased is "<<unphased_sizes[i]<<std::endl;
			test2 = true;
		}
	}			
}
