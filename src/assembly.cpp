#include <iostream>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <zlib.h>
#include "assembly.h"
//#include "msa.h"



void Assembly::generate_fasta_file(parameters& params, faidx_t*& fasta_index, std::set <std::string>& reads, std::string file_path)
{
	std::ofstream fp_write(file_path);	
	
	int loc_length;
	std::string line;	
	const size_t line_len = 60;
	std::set <std::string> read_seqs;

	for (auto &read: reads)
	{
		char* tmp = faidx_fetch_seq(fasta_index, read.c_str(), 0, MAX_FETCH_LEN, &loc_length);
		if (tmp == NULL)
			continue;

		std::string seq(tmp);

		fp_write<<">"<<read << std::endl;

		for (unsigned int i = 0; i < seq.size(); i += line_len)
		{
			std::string tmp = (seq.substr(i, line_len));
			fp_write<< tmp << std::endl;
    	}
		read_seqs.insert(seq);
	}

	fp_write.close();
}


int Assembly::write_svtigs(std::string& f_path, const std::string& f_name, int pos, std::string& contig, int coverage, std::ofstream& fp_write)
{
	std::ifstream fp_read(f_path);	
	std::string line;
	int contig_cnt = 0;
	std::string svtig_name = f_name;

	if (!fp_read)
	{
     	std::cerr << "Error opening "<<f_path<< std::endl;
		exit(1);
	}
	
	while(getline(fp_read, line))
	{
		if (line[0] == '>')
		{	
			if (contig_cnt > 0)
			{
				//continue;
				svtig_name = f_name + "_" + std::to_string(contig_cnt + 1);
			}
			
			contig_cnt++;
			if (contig == "" && pos == 0 && coverage == 0)
				fp_write<< ">"<<svtig_name<< std::endl;
			else
				fp_write<< ">"<<svtig_name<<" contig="<<contig<<" pos="<<pos<<" support="<<coverage<< std::endl;
		}
		else
			fp_write<< line << std::endl;
				
		line.clear();
	}
	return contig_cnt;
}


int Assembly::merge_svtigs(parameters& params)
{
	int cnt = 0;
	std::string asm_folder_path = params.log_path + "tmp/";
	
	//Find the files ending with ".cns.fa"
	for (const auto& entry : std::filesystem::directory_iterator(asm_folder_path))
	{
		if (entry.path().extension() == ".fa" && entry.path().stem().extension() == ".cns")
		{
			if (entry.file_size() == 0)
			{
				unassembled_cnt++;
				continue;
			}
			std::string file_path = entry.path();
			std::string file_name = entry.path().stem().stem();
			std::string tmp = "";
			cnt += write_svtigs(file_path, file_name, 0, tmp, 0,  params.fp_svtigs);

			if(params.debug)
			{
				std::string tmp_cmd = "cp " + file_path + " " + params.log_path + "out/";
				system(tmp_cmd.c_str());
			}
		}
	}
	return cnt;
}


//Assemble SV clusters
int Assembly::final_assembly(parameters& params, faidx_t*& fasta_index, std::set <std::string>& read_set, std::string& svtig_name, double& contig_depth, SVCluster*& sv, std::map <std::string, SVtig*>& final_svtigs)
{
	unsigned int support_threshold = static_cast<unsigned int>(params.support)/2;
	if (support_threshold < 3) support_threshold = 3;
	int svtig_tmp_cnt = 0;
	
	if (contig_depth * 2 < read_set.size())
	{
		this->filter_hicov++;
		return 0;
	}
	else if (contig_depth > 5 * read_set.size())
	{
		this->filter_lowcov++;
		return 0;
	}
	else if (read_set.size() < support_threshold)
	{
		this->filter_support++;
		return 0;
	}
	else if (contig_depth > MAX_CONTIG_DEPTH)
	{
		this->filter_hicov++;
		return 0;
	}

	//std::cout<<svtig_name<<" Read size: "<<read_set.size()<< " Contig depth: "<< contig_depth<<"\n";
	
	if (sv != nullptr)
	{
		SVtig *tmp = new SVtig;
		tmp->name = svtig_name;
		tmp->pos = sv->ref_pos;
		(tmp->reads).insert(read_set.begin(), read_set.end());
		tmp->contig = sv->contig;
		final_svtigs.insert(std::pair<std::string, SVtig*>(tmp->name, tmp));
	}

	std::string file_path = params.log_path + "tmp/" + svtig_name + ".fasta";	
	std::string output_path = params.log_path + "tmp/" + svtig_name;

	generate_fasta_file(params, fasta_index, read_set, file_path);
	int var_size = 4;
	//--ctg-min-nodes 2 -p 0 -k 15 -AS 2 --edge-min 1
	std::string wtdbg2_cmd = "wtdbg2.pl -t " + std::to_string(params.threads) + " -x ont -g " + std::to_string(var_size) + "m -o " + output_path + " " + file_path + " > /dev/null 2>&1";
	system(wtdbg2_cmd.c_str());
	
	//Write the assembly to the tmp fasta file
	svtig_tmp_cnt = merge_svtigs(params);

	//Remove the files in the folder
	for (const auto& entry : std::filesystem::directory_iterator(params.log_path + "tmp/")) 
	{
		if(params.debug)
		{
			if (entry.path() == file_path)
			{
				std::string tmp_cmd = "mv " + file_path + " " + params.log_path + "in/";
				system(tmp_cmd.c_str());
			}
		}
		std::filesystem::remove_all(entry.path());
	}	
	return svtig_tmp_cnt;
}


int Assembly::assemble_clusters(parameters& params, faidx_t*& fasta_index, std::vector<SVCluster*>& sv_cluster, std::map <std::string, Contig*>& depth, std::map <std::string, SVtig*>& final_svtigs)
{
	int initial_svtigs_cnt = 0;
	//Iterating over SV clusters of a contig
	for (auto &sv : sv_cluster) 
	{
		double contig_depth = depth[sv->contig]->coverage;
		if (contig_depth < 5)
			contig_depth = 5;
		
		std::string svtig_name;
		if ((params.phase_tags).empty())
		{
			//Create Svtigs for untagged reads
			svtig_name = sv->node + "_" + std::to_string(sv->start_pos);
			initial_svtigs_cnt += final_assembly(params, fasta_index, sv->reads_untagged, svtig_name, contig_depth, sv, final_svtigs);
		}
		else
		{
			svtig_name = "H1-" +sv->node + "_" + std::to_string(sv->start_pos);	
			initial_svtigs_cnt += final_assembly(params, fasta_index, sv->reads_h1, svtig_name, contig_depth, sv, final_svtigs);

			svtig_name = "H2-" +sv->node + "_" + std::to_string(sv->start_pos);	
			initial_svtigs_cnt += final_assembly(params, fasta_index, sv->reads_h2, svtig_name, contig_depth, sv, final_svtigs);
			
			if(!params.skip_untagged)
			{
				svtig_name = "None-" +sv->node + "_" + std::to_string(sv->start_pos);	
				initial_svtigs_cnt += final_assembly(params, fasta_index, sv->reads_untagged, svtig_name, contig_depth, sv, final_svtigs);	
			}
		}
	}

	return initial_svtigs_cnt;
}


void Assembly::run_assembly(parameters& params, std::map <std::string, Contig*>& depth, std::map<std::string, std::vector<SVCluster*>>& vars, std::set <std::string>& unmapped, std::map <std::string, SVtig*>& final_svtigs)
{
	int initial_svtigs_cnt = 0;
	std::map<std::string, std::vector<SVCluster*>>::iterator itr;

	auto t1 = std::chrono::steady_clock::now();
	std::cout<<"\nAssembly..."<<std::endl;		
	std::cout<<"--->assembling reads using "<<params.assembler<< std::endl;
		
	std::string svtigs_tmp_path = params.log_path + params.sample_name + "_svtigs_tmp.fa";
	params.fp_svtigs.open(svtigs_tmp_path);
	
	faidx_t* fasta_index = fai_load((params.fasta).c_str());

	//Iterating over contigs; i.e., check SV clusters of each contig
	for (itr=vars.begin(); itr != vars.end(); ++itr)
		initial_svtigs_cnt += assemble_clusters(params, fasta_index, itr->second, depth, final_svtigs);

	// Assemble unmapped reads
	double unmapped_count = static_cast<double>(unmapped.size());
	//std::cout <<"Size of unmapped is "<<unmapped_count<<"\n";
	if(unmapped_count > 0)
	{
		SVCluster* tmp = nullptr;
		std::string svtig_name = "None-unmapped_" + std::to_string(0);	
		initial_svtigs_cnt += final_assembly(params, fasta_index, unmapped, svtig_name, unmapped_count, tmp, final_svtigs);
	}
		
	if(std::filesystem::exists(params.log_path + "tmp/"))
		std::filesystem::remove_all(params.log_path + "tmp/");	

	params.fp_svtigs.close();
	fai_destroy(fasta_index);
	
	std::cout<<"--->"<<(filter_hicov)+(this->filter_lowcov)+(this->filter_support)<<" filtered (" <<this->filter_hicov<<" high, "<<this->filter_lowcov<<" low coverage read clusters and "<<this->filter_support<<" low read support)\n";
	
	params.fp_logs<<"--->"<<(this->filter_hicov)+(this->filter_lowcov)+(this->filter_support)<<" filtered (" <<this->filter_hicov<<" high, "<<this->filter_lowcov<<" low coverage read clusters and "<<this->filter_support<<" low read support)\n";
	
	std::cout<<"--->"<<unassembled_cnt<<" clusters cannot be assembled\n";	
	params.fp_logs<<"--->"<<unassembled_cnt<<" clusters cannot be assembled\n";

	std::cout<<"--->there are "<< initial_svtigs_cnt <<" svtigs before final filtering \n";
	params.fp_logs<<"--->there are "<< initial_svtigs_cnt <<" svtigs before final filtering \n";
	
	auto t2 = std::chrono::steady_clock::now();

	std::cout<<"--->assembly execution time: "<<std::chrono::duration<double> (t2 - t1).count()<<" sec.\n";
}
