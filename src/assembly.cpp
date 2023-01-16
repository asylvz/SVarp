#include <iostream>
#include <fstream>
#include <filesystem>
#include "assembly.h"


void generate_fastq_file(parameters* params, std::map<std::string, unsigned long> fasta_index, std::set <std::string> reads, std::string file_path)
{
	//std::cout<<"Generating fastq file"<<std::endl;
	
	std::ifstream fp_read(params->fasta);
	std::ofstream fp_write(file_path);
	
	unsigned long char_pos;	
	std::string line;	

	/*for(auto &a : fasta_index)
	{
		std::cout << a.first << " " << a.second<<"\n";
		if (fp_write.is_open())
  		{	
			fp_write << a.first << " " << a.second<<"\n";
		}
		else
			std::cout << "Unable to open file";
	}
	fp_write.close();*/
	
	for (auto &read: reads)
	{
		if (fasta_index.find(read)!=fasta_index.end())
		{
			char_pos = fasta_index[read];

			//std::cout<<"READ = \t"<<read<<" - " <<char_pos <<std::endl;
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
	std::cout<<"Indexing the fasta file"<<std::endl;
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


void run_assembly(parameters* params, std::multimap<std::string, variant*>& insertions)
{
	std::map<std::string, unsigned long> fasta_index;
	std::multimap<std::string, variant*>::iterator itr;

	std::string cwd = std::filesystem::current_path().string();
		
	std::string log_path = cwd + "/log/";
	if(std::filesystem::exists(log_path))
		std::filesystem::remove_all(log_path);

  	if (!std::filesystem::create_directory(log_path))
	{
		std::cerr << "Error creating log folder" << std::endl;
		exit(-1);
  	}
  	if (!std::filesystem::create_directory(log_path + "assembly_input/"))
	{
		std::cerr << "Error creating log/assembly_input/" << std::endl;
		exit(-1);
  	}
  	if (!std::filesystem::create_directory(log_path + "assembly_output/"))
	{
		std::cerr << "Error creating log/assembly_output/" << std::endl;
		exit(-1);
  	}

	std::cout<<"Assembly using wtdbg2..."<<std::endl;		

	index_fasta(params, fasta_index);	

	for (itr=insertions.begin(); itr != insertions.end(); ++itr)
	{
		if (itr->second->reads.size() == 0)
			std::cout<< "ALARM - there is no read supporting the SV ("<<itr->second->ref_start<<" - "<<itr->second->ref_end<<") in "<< itr->second->contig <<std::endl;
		
		//Generate fastq files
		if (itr->second->reads.size() > 2)
		{
			std::string filename = itr->second->contig + "_" + std::to_string(itr->second->ref_start) + "_" + std::to_string(itr->second->ref_end);
				
			std::string file_path = log_path + "assembly_input/" + filename + ".fasta";
			generate_fastq_file(params, fasta_index, itr->second->reads, file_path);	
			std::string output_path = log_path + "assembly_output/" + filename + "/";

  			if (! std::filesystem::create_directory(output_path)) {
        		std::cerr << "Error creating the folder "<<output_path << std::endl;
				exit(-1);
  			}

			int var_size = itr->second->sv_size / 1000;
			if(var_size == 0)
				var_size = 1;

			std::string wtdbg2_cmd1 = "wtdbg2 -t 16 -x ont -g " + std::to_string(var_size) + "k -o " + output_path;			
			
			std::string wtdbg2_cmd2 = "wtpoa-cns -t 16 -i " + output_path + ".ctg.lay.gz -fo " + output_path + filename + ".fa";
			
			//std::cout<< wtdbg2_cmd1<<"\n"<<wtdbg2_cmd2<<"\nSV size=" << itr->second->sv_size<< std::endl;
		system(wtdbg2_cmd1.c_str());
		system(wtdbg2_cmd2.c_str());
		}
	}



    /*for (const auto & entry : std::filesystem::directory_iterator(path + "assembly_input/"))
	{
		if (entry.path().extension() != ".fasta")
			continue;
		
		std::string rawname = entry.path().string().substr(entry.path().string().find_last_of("/\\") + 1);	
		std::string output_path = path + "assembly_out/" + rawname.substr(0, rawname.find_last_of(".")); 

		std::string flye_cmd = "flye --nano-raw " + entry.path().string() + " --out-dir " + output_path;
		std::cout<<"Running flye with " <<flye_cmd<<std::endl;
		system(flye_cmd.c_str());

	}*/

}
