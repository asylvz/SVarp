#include <iostream>
#include <filesystem>
#include <fstream>
#include <getopt.h>
#include "cmdline.h"


int parse_command_line(int argc, char** argv, parameters& params)
{
	int index, o;
	std::string support, dist_threshold, threads, as, pc, map_ratio;
	
	static struct option long_options[] = 
	{	
		{"gaf" , required_argument, NULL, 'a'},
		{"assembler" , required_argument, NULL, 'b'},
		{"pc" , required_argument, NULL, 'c'},
		{"dist-threshold" , required_argument, NULL, 'd'},
		{"as" , required_argument, NULL, 'e'},
		{"fasta" , required_argument, NULL, 'f'},
		{"graph" , required_argument, NULL, 'g'},
		{"help"   , no_argument, 0, 'h'},
		{"skip-untagged" , no_argument, 0, 'j'},
		{"asm" , no_argument, NULL, 'm'},
		{"map-ratio" , required_argument, NULL, 'n'},
		{"out" , required_argument, NULL, 'o'},
		{"phase" , required_argument, NULL, 'p'},
		{"no-remap" , no_argument, 0, 'r'},
		{"support" , required_argument, NULL, 's'},
		{"sample" , required_argument, NULL, 'i'},
		{"threads" , required_argument, NULL, 't'},
		{"debug" , no_argument, NULL, 'u'},
		{NULL, 0, NULL, 0}
	};
	
	while((o = getopt_long( argc, argv, "a:b:c:d:e:f:g:h:i:j:m:n:o:p:r:s:t:u:", long_options, &index)) != -1)
	{
		switch(o)
		{
			case 'a':
				params.gaf = optarg;
				break;
			case 'b':
				params.assembler = optarg;
				break;
			case 'c':
				pc = optarg;
				break;	
			case 'd':
				dist_threshold = optarg;
				break;	
			case 'e':
				as = optarg;
				break;	
			case 'f':
				params.fasta = optarg;
				break;
			case 'g':
				params.ref_graph = optarg;
				break;
			case 'i':
				params.sample_name = optarg;
				break;
			case 'j':
				params.skip_untagged = true;
				break;
			case 'm':
				params.asm_mode = true;
				break;
			case 'n':
				map_ratio = optarg;
				break;	
			case 'o':
				params.output_path = optarg;
				break;	
			case 'p':
				params.phase_tags = optarg;
				break;
			case 'r':
				params.no_remap = true;
				break;
			case 's':
				support = optarg;
				break;
			case 't':
				threads = optarg;
				break;
			case 'u':
				params.debug = true;
				break;
			case 'h':
				print_help();
				exit(0);
		}
	}	
	
	std::cout<<"\n";
	/* check if --ref   is invoked */
	if((params.ref_graph).empty())
	{
		std::cerr<<"[SVARP CMDLINE ERROR] Please enter reference graph genome file (GFA) using the --graph (-g) option.\n";
		return RETURN_ERROR;
	}

	if((params.gaf).empty())
	{
		std::cerr<<"[SVARP CMDLINE ERROR] Please enter alignment file (GAF) using the --gaf (-a) option.\n";
		return RETURN_ERROR;
	}
	if((params.fasta).empty())
	{
		std::cerr<<"[SVARP CMDLINE ERROR] Please enter the bgzipped fasta file path using the --fasta (-f) option. Or do NOT use exhaustive parameter\n";
		return RETURN_ERROR;
	}
	
	if((params.assembler).empty() || params.assembler == "wtdbg" || params.assembler == "wtdbg2")
	{
		params.assembler = "wtdbg2";
		std::cerr<<"Using wtdbg2 for assembly...\n";
	}
	else if (params.assembler == "abpoa" || params.assembler == "poa")
	{
		params.assembler = "abpoa";
		std::cerr<<"Using abPOA for assembly...\n";
	}
	else if(params.assembler == "shasta" || params.assembler == "Shasta")
	{
		params.assembler = "Shasta";
		std::cerr<<"Using Shasta for assembly...\n";
	}
	else
	{
		std::cerr<<"SVarp does not use "<<params.assembler<< " assembler. Switching to wtdbg2\n";
		params.assembler = "wtdbg2";
	}
			
	if(params.no_remap)
		std::cerr<<"Skipping realignment (not suggested, you may have too many false positives and duplicates)...\n";
		
	if((params.output_path).empty())
	{
		std::string cwd = std::filesystem::current_path().string();
		params.log_path = cwd + "/svarp_log/";	
		//std::cerr<<"No output path provided (--out).\n";
	}
	else
	{
		if(!std::filesystem::exists(params.output_path))
		{
			std::cerr<<params.output_path<<" does not exist, generating...\n";
			std::filesystem::create_directory(params.output_path);
		}

		if (params.output_path.at(params.output_path.size()-1) == '/')
			params.log_path = params.output_path;
		else
			params.log_path = params.output_path + "/";
	}

	if((params.sample_name).empty())
		params.sample_name = "sample";
	
	if(dist_threshold.empty())
		params.dist_threshold = 100;
	else
	{
		try {
			params.dist_threshold = stoi(dist_threshold);
		} catch (const std::invalid_argument&) {
			std::cerr << "[SVARP CMDLINE ERROR] dist_threshold must be an integer: " << dist_threshold << std::endl;
			return RETURN_ERROR;
		}
	}
	
	if(as.empty())
		params.min_alignment_score = 5000;
	else
	{
		try {
			params.min_alignment_score = stoi(as);
		} catch (const std::invalid_argument&) {
			std::cerr << "[SVARP CMDLINE ERROR] alignment_score must be an integer: " << as << std::endl;
			return RETURN_ERROR;
		}
	}

	if(pc.empty())
		params.min_precise_clipping = 0.97;
	else
	{
		try {
			params.min_precise_clipping = stod(pc);
		} catch (const std::invalid_argument&) {
			std::cerr << "[SVARP CMDLINE ERROR] precise_clipping must be a float: " << pc << std::endl;
			return RETURN_ERROR;
		}
	}

	if(map_ratio.empty())
		params.min_map_ratio = 0.90;
	else
	{
		try {
			params.min_map_ratio = stod(map_ratio);
		} catch (const std::invalid_argument&) {
			std::cerr << "[SVARP CMDLINE ERROR] map_ratio must be a float: " << map_ratio << std::endl;
			return RETURN_ERROR;
		}
	}

	if((params.phase_tags).empty())
		std::cerr<<"No phase information provided (--phase). SVs will not be phased...\n";

	if(support.empty())
		params.support = 5;
	else
	{
		try {
			params.support = stoi(support);
		} catch (const std::invalid_argument&) {
			std::cerr << "[SVARP CMDLINE ERROR] support must be an integer: " << support << std::endl;
			return RETURN_ERROR;
		}
	}
	
	if(params.asm_mode)
	{
		std::cerr<<"Using assembly mode...\n";
		params.support = 0;
		params.dist_threshold = 0;
		params.vcf_path = params.log_path + params.sample_name + ".vcf";
	}

	if(params.debug)
	{
		std::cerr<<"Running in debug mode...\n";
	}

	if(threads.empty())
		params.threads = 16;
	else
	{
		try {
			params.threads = stoi(threads);
		} catch (const std::invalid_argument&) {
			std::cerr << "[SVARP CMDLINE ERROR] threads must be an integer: " << threads << std::endl;
			return RETURN_ERROR;
		}
	}

	// Validate input files exist
	if (!std::filesystem::exists(params.gaf))
	{
		std::cerr << "[SVARP CMDLINE ERROR] GAF file not found: " << params.gaf << std::endl;
		return RETURN_ERROR;
	}

	if (!std::filesystem::exists(params.ref_graph))
	{
		std::cerr << "[SVARP CMDLINE ERROR] GFA file not found: " << params.ref_graph << std::endl;
		return RETURN_ERROR;
	}

	if (!std::filesystem::exists(params.fasta))
	{
		std::cerr << "[SVARP CMDLINE ERROR] FASTA file not found: " << params.fasta << std::endl;
		return RETURN_ERROR;
	}

	return RETURN_SUCCESS;
}


void init_logs(parameters& params)
{	
	
	std::cout<<"\n...hallo, merhaba, ola, salaam, hello!!! SVarp is running...\n";
	if(std::filesystem::exists(params.log_path))
		std::filesystem::remove_all(params.log_path);
  	
	try {
		if (!std::filesystem::create_directory(params.log_path))
			error("Error creating log folder");
  	
		if (params.debug)
		{
  			if (!std::filesystem::create_directory(params.log_path + "in/"))
				error("Error creating log/in/ directory");
  		
  			if (!std::filesystem::create_directory(params.log_path + "out/"))
				error("Error creating log/out/ directory");
		}
	
  		if (!std::filesystem::create_directory(params.log_path + "tmp/"))
			error("Error creating log/tmp/ directory");
	} catch (const std::filesystem::filesystem_error& e) {
		std::string msg = std::string("Error creating log directories: ") + e.what();
		error(msg.c_str());
	}
	params.fp_logs.open(params.log_path + params.sample_name + ".log");
		
	params.remap_gaf_path = params.log_path + params.sample_name + "_remap.gaf";


	std::cout<<"\nParameters:\n\tMinimum read support: "<<params.support<<"\n\tMinimum distance threshold: "<< params.dist_threshold<<"\n"<<std::endl;
	std::cout<<"\tMinimum map ratio: "<<params.min_map_ratio<<"\n\tPrecise clipping (Graphaligner): "<< params.min_precise_clipping<<"\n\tAlignment score (Graphaligner): "<< params.min_alignment_score<<"\n"<<std::endl;
	std::cout<<"\nInput files:\n\t"<<params.gaf<<"\n\t"<< params.ref_graph<<"\n\t"<<params.fasta<<std::endl;
	std::cout<<"\nLog folder:\n\t"<<params.log_path<<std::endl;

	
	if (params.fp_logs.is_open()) {
		params.fp_logs << "\nParameters:\n\tMinimum read support: " << params.support << "\n\tMinimum distance threshold: " << params.dist_threshold << "\n\n";
		params.fp_logs << "\tMinimum map ratio: " << params.min_map_ratio << "\n\tPrecise clipping (Graphaligner): " << params.min_precise_clipping << "\n\tAlignment score (Graphaligner): " << params.min_alignment_score << "\n\n";
		params.fp_logs << "\nInput files:\n\t" << params.gaf << "\n\t" << params.ref_graph << "\n\t" << params.fasta << "\n";
		params.fp_logs << "\nLog folder:\n\t" << params.log_path << "\n\n";
	}
}


void print_help() 
{
	std::cerr << std::endl;
	std::cout << "SVarp: pangenome-based structural variation discovery" << std::endl;
#ifndef SVARP_VERSION
#define SVARP_VERSION "unknown"
#endif
#ifndef SVARP_UPDATE
#define SVARP_UPDATE "unknown"
#endif
	std::cout<< "\tVersion "<<SVARP_VERSION<<", Last update: "<<SVARP_UPDATE<<"\n";
	std::cerr << std::endl;
	std::cerr << "Required arguments"<<std::endl;
	std::cerr << "\t--gaf (-a)                  : GAF alignment file"<<std::endl;
	std::cerr << "\t--graph (-g)                : GFA pangenome file"<<std::endl;
	std::cerr << "\t--fasta (-f)                : Fasta sequence file"<<std::endl;
	std::cerr << std::endl;
	std::cerr << "Optional arguments"<<std::endl;
	std::cerr << "\t--sample (-i)               : Sample name."<<std::endl;
	std::cerr << "\t--out (-o)                  : Output folder."<<std::endl;
	std::cerr << "\t--debug                     : Output multiple log files for debugging purpose."<<std::endl;
	std::cerr << "\t--skip-untagged             : Output only phased variants (~30\% faster)."<<std::endl;
	std::cerr << "\t--dist_threshold (-d)       : Distance threshold to merge SV breakpoints (default=100)"<<std::endl;
	std::cerr << "\t--out (-o)                  : Output folder path"<<std::endl;
	std::cerr << "\t--phase (-p)                : WhatsHap haplotag file in .tsv (https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-haplotag)"<<std::endl;
	std::cerr << "\t--reads(-r)                 : Bgzipped FASTA file of reads for extensive mode (needed for WFA realignment)"<<std::endl;
	std::cerr << "\t--sample (-i)               : Sample (Individual) name"<<std::endl;
	std::cerr << "\t--support (-s)              : Minimum support for a cluster to be assembled (default=5 for diploid samples)"<<std::endl;
	std::cerr << "\t--threads(-t)               : Number of threads for assembly and realignment (default:32)"<<std::endl;
	//std::cerr << "\t--assembler                 : (Experimental) Assembler can be either \"Shasta\", \"wtdbg2\" or \"bsalign\" (default:wtdbg2) "<<std::endl;
	//std::cerr << "\t--asm                       : (Experimental) Runs in assembly mode. You can provide assembly to find the variations. This outputs exact breakpoints instead of SVtigs"<<std::endl;
	std::cerr << "\t--help                      : Print this help menu"<<std::endl;

	std::cerr << "\nCommand to run\n\tbuild/svarp -a xxx.gaf -g xxx.gfa --fasta xxx.fasta.gz --phase read_tags.tsv -i SAMPLE_NAME -o OUTPUT_FOLDER\n\n";
	std::cerr << std::endl;
}

