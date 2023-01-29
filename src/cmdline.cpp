#include <iostream>
#include <filesystem>
#include <fstream>
#include <getopt.h>
#include "cmdline.h"


std::ofstream logFile;

int parse_command_line(int argc, char** argv, parameters* params)
{
	int index, o;
	
	static struct option long_options[] = 
	{
		{"gaf" , required_argument, NULL, 'a'},
		{"fasta" , required_argument, NULL, 'f'},
		{"graph" , required_argument, NULL, 'g'},
		{"help"   , no_argument,         0, 'h'},
		{"phase" , required_argument, NULL, 'p'},
		{NULL, 0, NULL, 0}
	};
	
	while((o = getopt_long( argc, argv, "a:f:g:p:h:", long_options, &index)) != -1)
	{
		switch(o)
		{
			case 'a':
				params->gaf = optarg;
				break;	
			case 'f':
				params->fasta = optarg;
				break;
			case 'g':
				params->ref_graph = optarg;
				break;
			case 'p':
				params->phase_tags = optarg;
				break;
			case 'h':
				print_help();
				exit(0);
		}
	}	
	/* check if --ref   is invoked */
	if((params->ref_graph).empty())
	{
		std::cerr<<"[SVAPAN CMDLINE ERROR] Please enter reference graph genome file (GFA) using the --graph (-g) option.\n";
		return RETURN_ERROR;
	}

	if((params->gaf).empty())
	{
		std::cerr<<"[SVAPAN CMDLINE ERROR] Please enter alignment file (GAF) using the --gaf (-a) option.\n";
		return RETURN_ERROR;
	}	
	
	if((params->phase_tags).empty())
	{
		std::cerr<<"\nNo phase information provided (--phase). SVs will not be phased...\n";
		return RETURN_SUCCESS;
	}	

	return RETURN_SUCCESS;
}


void init_logs(parameters* params)
{	
	std::string cwd = std::filesystem::current_path().string();
	std::string log_path = cwd + "/log/";	
	std::cout<<"\n...Hello, svapan is running...\n";
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

	logFile.open(log_path + "svapan.log");
	
	std::cout<<"\nInput files:\n\t"<<params->gaf<<"\n\t"<< params->ref_graph<<"\n\t"<<params->fasta<<std::endl;
	std::cout<<"\nLog folder:\n\t"<<log_path<<std::endl;
}


void print_help() 
{
	std::cerr << std::endl;
	std::cout << "svapan: Phased structural variation discovery in pangenomes" << std::endl;
	
	std::cerr << std::endl;
	std::cerr << "Required arguments"<<std::endl;
	std::cerr << "\t--gaf (-a)          : GAF alignment file"<<std::endl;
	std::cerr << "\t--graph (-g)        : GFA pangenome file"<<std::endl;
	std::cerr << "\t--fasta (-f)        : Fasta sequence file"<<std::endl;
	std::cerr << "\t--phase (-p)        : WhatsHap haplotag output file (in .tsv)"<<std::endl;
	std::cerr << std::endl;
	std::cerr << "Optional arguments"<<std::endl;
	std::cerr << "\t--help              : Print this help menu"<<std::endl;
	std::cerr << std::endl;
}

