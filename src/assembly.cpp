#include <iostream>
#include <fstream>
#include <filesystem>
#include "assembly.h"


void generate_fastq_file(parameters* params, std::map<std::string, unsigned long>& fasta_index, std::set <std::string>& reads, std::string file_path)
{
	std::ifstream fp_read(params->fasta);
	std::ofstream fp_write(file_path);
	
	unsigned long char_pos;	
	std::string line;	
	
	for (auto &read: reads)
	{
		if (fasta_index.find(read)!=fasta_index.end())
		{
			char_pos = fasta_index[read];

			fp_read.seekg(char_pos, std::ios::beg);
			
			getline(fp_read, line);
			fp_write<< line << std::endl;
			line.clear();
			while(getline(fp_read, line))
			{
				if(!fp_read || line[0] == '>')
					break;
				fp_write<< line << std::endl;
				line.clear();
			}
		}
		else
			std::cout<<"Not found"<<read<<std::endl;	
	}
}

int index_fasta(parameters* params, std::map<std::string, unsigned long>& fasta_index)
{
	std::cout<<"--->indexing the fasta file"<<std::endl;
	std::ifstream fp_read(params->fasta);
	std::string line, read_name;
	unsigned long char_count = 0;

	if(!fp_read.good())
	{
        std::cerr << "Error opening '"<<params->fasta<< std::endl;
        return RETURN_ERROR;
    }
	while(std::getline(fp_read, line))
	{
        if(line[0] == '>')
		{
			int pos = line.find(" ");
        	read_name = line.substr(1, pos - 1);

			fasta_index.insert(std::pair<std::string, unsigned long>(read_name, char_count));
			read_name.clear();
		}
		// +1 is for the newline character
		char_count += line.length() + 1;
    }
	return RETURN_SUCCESS;
}


void run_assembly(parameters* params, std::map<std::string, std::vector<svtig*>>& insertions)
{
	std::map<std::string, unsigned long> fasta_index;
	std::map<std::string, std::vector<svtig*>>::iterator itr;

	std::string cwd = std::filesystem::current_path().string();
		
	std::string log_path = cwd + "/log/";

	std::cout<<"\nAssembly..."<<std::endl;		
	
	index_fasta(params, fasta_index);	

	std::cout<<"--->assembling the reads using wtdbg2"<<std::endl;
	int svtig_cnt = 0;
	int h1 = 0, h2 = 0;
	for (itr=insertions.begin(); itr != insertions.end(); ++itr)
	{	
		//Generate fastq files
		for (auto &sv : itr->second) 
		{
			if ((sv->reads_h1).size() < MIN_READ_SUPPORT)
				continue;

			std::string filename = sv->contig + "_svtig" + std::to_string(svtig_cnt++) + "_H1";
				
			//logFile<<"H1 "<<sv->reads_h1.size()<<" "<<itr->first <<std::endl;
			//for (auto &a: sv->reads_h1)
			//	logFile<<"\t"<<a<<std::endl;

			std::string file_path = log_path + "in/" + filename + ".fasta";	
			std::string output_path = log_path + "out/" + filename;
			
			generate_fastq_file(params, fasta_index, sv->reads_h1, file_path);	
			
  			if (! std::filesystem::create_directory(output_path))
        		std::cerr << "Error creating the folder "<<output_path << std::endl;

			int var_size = (sv->sv_size * 2) / 1000;
			if(var_size == 0)
				var_size = 1;
			
			h1++;
			std::string wtdbg2_cmd = "wtdbg2.pl -t 16 -x ont -g " + std::to_string(var_size) + "m -o " + output_path + " " + file_path + " &>"+ output_path+".log";
			
			//std::cout<<"\n"<<wtdbg2_cmd<<"\nSV size=" << sv->sv_size<< std::endl;
			//system(wtdbg2_cmd.c_str());
		}

		for (auto &sv : itr->second) 
		{

			if ((sv->reads_h2).size() < MIN_READ_SUPPORT)
				continue;
			
			std::string filename = sv->contig + "_svtig" + std::to_string(svtig_cnt++) + "_H2";
				
			//logFile<<"H2"<<sv->reads_h1.size()<<" "<<itr->first <<std::endl;
			//for (auto &a: sv->reads_h1)
			//	logFile<<"\t"<<a<<std::endl;

			std::string file_path = log_path + "in/" + filename + ".fasta";	
			std::string output_path = log_path + "out/" + filename;
			
			generate_fastq_file(params, fasta_index, sv->reads_h2, file_path);	
			
  			if (! std::filesystem::create_directory(output_path))
        		std::cerr << "Error creating the folder "<<output_path << std::endl;

			int var_size = (sv->sv_size * 2) / 1000;
			if(var_size == 0)
				var_size = 1;
			
			h2++;
			std::string wtdbg2_cmd = "wtdbg2.pl -t 16 -x ont -g " + std::to_string(var_size) + "m -o " + output_path + " " + file_path + " &>"+ output_path+".log";	
			
			//std::cout<<"\n"<<wtdbg2_cmd<<"\nSV size=" << sv->sv_size<< std::endl;
			//system(wtdbg2_cmd.c_str());
		}
	}
	std::cout<<"\nThere are "<<h1+h2<<" SVs that have >"<<MIN_READ_SUPPORT<<" minimum read support ("<<h1<<" H1 - "<<h2<<" H2)" <<std::endl;
}

