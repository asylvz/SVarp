#include <iostream>
#include <filesystem>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <iterator>
#include <chrono>
#include <algorithm>
#include <htslib/faidx.h>
#include <sys/wait.h>
#include <unistd.h>   // access, X_OK
#include <cstdlib>
#include "common.h"
#include "remap.h"
#include "alignment.h"
#include "variant.h"
#include "reference.h"
#include "bindings/cpp/WFAligner.hpp"



inline bool cmp(Read* i1, Read* i2)
{
	if (i1->node != i2->node)
		return (i1->node < i2->node);
	else
	{
		if (i1->start != i2->start)
    		return (i1->start < i2->start);
		else
			if(i1->end != i2->end)
				return (i1->end < i2->end);
			else
				return (i1->end > i2->end);
	}
}


int get_middle_string(std::string& s)
{
	int y = (s.length() - (s.length()/10)) / 2;
	s.erase(0, y);
	s.erase(s.length() - y, s.length());
	return y;
}


int write_final_svtigs_fasta(faidx_t*& fasta_index, std::string& svtig_name, int pos, std::string& contig, int coverage, std::ostream& fp_write)
{
	
	int loc_length;
	std::string line;	
	const size_t line_len = 60;

	std::string seq = faidx_fetch_seq(fasta_index, svtig_name.c_str(), 0, MAX_FETCH_LEN, &loc_length);
	//int new_pos = get_middle_string(seq);
	fp_write<< ">"<<svtig_name<<" contig="<<contig<<" pos="<<pos<<" support="<<coverage<< std::endl;

	for (unsigned int i = 0; i < seq.size(); i += line_len)
	{
		std::string tmp = (seq.substr(i, line_len));
		fp_write<< tmp << std::endl;
    }
	return RETURN_SUCCESS;
}


std::pair<int, int> remove_duplicates(std::vector <Read*>& tmp_svtig, std::map <std::string, SVtig*>& final_svtigs, int& extra_added)
{
	std::pair<int, int> dup_legit(0,0);
	std::map<std::string, SVtig*>::iterator it_svtigs;
	std::sort(tmp_svtig.begin(), tmp_svtig.end(), cmp);
	
	for(unsigned int i = 0; i < tmp_svtig.size(); i++)
	{
		bool toAdd = true;
		Read *r = tmp_svtig[i];
		if(r->duplicate)
			continue;

		for(unsigned int j = 0; j < tmp_svtig.size(); j++)
		{
			if(tmp_svtig[i]->node != tmp_svtig[j]->node)
				continue;

			if(tmp_svtig[j]->duplicate || i == j)
				continue;
			
			double overlap = overlap_ratio(tmp_svtig[i]->start, tmp_svtig[i]->end, tmp_svtig[j]->start, tmp_svtig[j]->end);
			
			if (overlap > 0.9)
			{
				if(tmp_svtig[i]->svtig_size >= tmp_svtig[j]->svtig_size)
				{
					tmp_svtig[j]->duplicate = true;
					continue;	
				}
				else
				{
					//std::cout<<"\t duplicate and longer, so do not add the inn\n";
					r->duplicate = true;
					dup_legit.first++;
					toAdd = false;
					break;
				}
			}
		}
		//Add if not a duplicate
		if (toAdd)
		{
			//If the assembly is fragmented (i.e., multiple contigs), then we need to add the additional contig because it is not in final_svtigs DS (only None-s237302_123 is in final_svtigs). e.g., None-s237302_123 and None-s237302_123_2
			auto npos = r->rname.rfind("_");
			auto npos2 = r->rname.find("_");
			std::string tmp_str = "";
			if (npos != npos2)
			{
				tmp_str = r->rname.substr(0, npos);
				it_svtigs = final_svtigs.find(tmp_str);
				
				//Add it to the final_svtigs
				if (it_svtigs != final_svtigs.end())
				{
					SVtig *tmp = new SVtig;
					tmp->name = r->rname;
					tmp->pos = it_svtigs->second->pos;
					(tmp->reads).insert(it_svtigs->second->reads.begin(), it_svtigs->second->reads.end());
					tmp->contig = it_svtigs->second->contig;
					tmp->output = true;
					final_svtigs.insert(std::pair<std::string, SVtig*>(r->rname, tmp));
					extra_added++;
				}
				else
					std::cout<<"Error - SVtig= "<<r->rname<<" not found (multiple contig)...\n";
			}
			else
			{
				std::map<std::string, SVtig*>::iterator it_dup = final_svtigs.find(r->rname);
				if (it_dup != final_svtigs.end())
					it_svtigs->second->output = true;
				else
					std::cout<<"Error - SVtig= "<<r->rname<<" not found...\n";
			}
			dup_legit.second++;
			//std::cout<<"Legit, added\n";
		}
	}
	//Add the unmapped reads as well
	
	/*bool is_aligned = false;
	for (it_svtigs=final_svtigs.begin(); it_svtigs != final_svtigs.end(); ++it_svtigs)
	{
		std::cout<<"Checking "<<it_svtigs->first<<"\t";
		is_aligned = false;
		for (auto & element : tmp_svtig)
		{
			if (it_svtigs->first == element->rname)
			{
				std::cout<<"Aligned\n";
				is_aligned = true;
				break;
			}	
		}
		if (!is_aligned)
		{
			std::cout<<"Unmapped\n";
			it_svtigs->second->output = true;
		}
	}*/

	return dup_legit;
}


void wfa_align(std::map<std::string, gfaNode*>& gfa, std::string& cigar, std::string &query_name, int query_start, int query_end, std::string path, int gfa_start, int gfa_end, wfa::WFAlignerGapAffine &aligner, faidx_t*& fasta_index)
{
	std::string ref_tmp = "";
	
	size_t p = 0;
	while (p < path.size())
	{
		char strand = path[p];
		++p;
		size_t q = p;
		while (q < path.size() && path[q] != '>' && path[q] != '<') ++q;
		std::string node_name = path.substr(p, q - p);
		p = q;
		
		if (strand == '>')
			ref_tmp += gfa[node_name]->sequence;
		else if (strand == '<')
			ref_tmp += reverse_complement(gfa[node_name]->sequence);
		else
			std::cout<<"Strand resolution issue in wfa_align()\n";
	}
	
	std::string ref = ref_tmp.substr(gfa_start, gfa_end - gfa_start + 1);

	int loc_length;	
	std::string query = faidx_fetch_seq(fasta_index, query_name.c_str(), query_start, query_end, &loc_length);
	
	aligner.alignEnd2End(query, ref); // Align
	cigar = aligner.getCIGAR(true);
}


int read_remappings(parameters& params, std::map<std::string, gfaNode*>& gfa, std::map <std::string, SVtig*>& final_svtigs, faidx_t*& fasta_index)
{
	int secondary = 0, primary = 0, lowmq = 0, extra_added = 0;

	std::string line;
	std::vector <std::string> tokens;
	std::ifstream fp;
	
	std::cout<<"--->reading remappings from "<< params.remap_gaf_path <<std::endl;	
	fp.open(params.remap_gaf_path);

	if(!fp.good())
	{
		std::cerr << "Error opening '"<<params.remap_gaf_path << std::endl;
        return RETURN_ERROR;
    }

	//Read the gaf once to find the largest mappings of the svtigs. If that's more than half of the svtig, get it	
	std::map<std::string, Read*>::iterator it;
	std::map <std::string, Read*> reads;
	std::map <std::string, std::string> non_dup_by_pos, non_dup_by_rname;
	std::vector<int> cigarLen;
	std::vector<char> cigarOp;
	
	wfa::WFAlignerGapAffine aligner(4, 6, 2, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryMed);
	
	double map_ratio = 0.0;
	while(getline(fp, line))
	{
		Gaf g;
		parse_gaf_line(line, g);
		
		bool hasSV = false, LowMQ = false;	
		
		if(g.mapping_quality < MINMAPQREMAP)
		{
			lowmq++;
			LowMQ = true;
		}

		if (!g.is_primary)
			secondary++;
		else if (g.is_primary && !LowMQ)
			primary++;
		

		cigarLen.clear();
		cigarOp.clear();

		std::string cigar;
				
		//Use the cigar of WFA realignment
		wfa_align(gfa, cigar, g.query_name, g.query_start, g.query_end, g.path, g.path_start, g.path_end, aligner, fasta_index);

		int cigar_cnt = decompose_cigars(cigar, cigarLen, cigarOp);
		for (int c = 0; c < cigar_cnt; c++)
		{
			if ((cigarOp[c] == INSERTION || cigarOp[c] == DELETION || cigarOp[c] == MISMATCH) && cigarLen[c] > MINSVSIZE)
			{
				hasSV = true;
				break;
			}
		}

		map_ratio = static_cast<double> ((double) g.query_end - g.query_start) / g.query_length;

		it = reads.find(g.query_name);
		if (it != reads.end())
		{
			if(LowMQ || !g.is_primary)
			{
				if((it->second)->sv_in_cigar == false && hasSV)
				{
					(it->second)->sv_in_cigar = true;
				
					if(map_ratio > (it->second)->highest_map_ratio)
						(it->second)->highest_map_ratio = map_ratio;

					(it->second)->start = g.path_start;
					(it->second)->end = g.path_end;
					(it->second)->node = g.path;
				}
			}
			else
			{
				if(map_ratio > (it->second)->highest_map_ratio)
					(it->second)->highest_map_ratio = map_ratio;

				if(((it->second)->end - (it->second)->start) < (g.path_end - g.path_start))
				{
					(it->second)->start = g.path_start;
					(it->second)->end = g.path_end;
					(it->second)->node = g.path;
				}

				(it->second)->freq++;
			}
		}
		else
		{	
			Read *r = new Read();
			r->highest_map_ratio = map_ratio;
			r->freq = 1;
			r->rname = g.query_name;
			r->node = g.path;
			r->start = g.path_start;
			r->end = g.path_end;
			r->sv_in_cigar = hasSV;
			r->svtig_size = g.query_length;
			reads.insert(std::pair<std::string, Read*>(g.query_name, r));
		}
	}

	std::vector <Read*> tmp_svtig;
	int filtered = 0;
	
	for (auto &t: reads)
	{
		if(t.second->highest_map_ratio > 0)
		{	
			if(t.second->freq > 1)
				tmp_svtig.push_back(t.second);
			else if(t.second->highest_map_ratio < params.min_map_ratio)
				tmp_svtig.push_back(t.second);
			else if (t.second->sv_in_cigar)
				tmp_svtig.push_back(t.second);
			else
				filtered++;
		}
		else
			filtered++;
	}
	//If an svtig is in tmp_svtig, check if it is a duplicate. If not set "->output=true"
	//Check the unmapped svtigs and add them tooo
	std::pair<int, int> dup_legit = remove_duplicates(tmp_svtig, final_svtigs, extra_added);

	std::cout<<"--->"<< filtered << " filtered - " << dup_legit.first<<" duplicate\n";
	std::cout<<"--->"<< primary<<" primary, "<<secondary<<" secondary mappings, "<<lowmq<<" low MAPQ(<"<<MINMAPQREMAP<<")" <<"; svtigs from multiple contig assemblies = "<<extra_added<<"\n";

	return RETURN_SUCCESS;
}


int write_final_svtigs(faidx_t*& fasta_index, std::map <std::string, SVtig*>& final_svtigs, std::string& out_file, std::string haplotype)
{
	int cnt = 0;
	std::map<std::string, SVtig*>::iterator itr;
	std::ofstream fp_write(out_file);
	
	//Find the files ending with ".cns.fa"
	std::string file_name = "";

	//params.fp_logs << "\n\n------->Reads contributing to each SVtig\n\n";
	
	for (itr=final_svtigs.begin(); itr != final_svtigs.end(); ++itr)
	{
		if (itr->second->output == true)
		{
			file_name = itr->second->name;
			if (haplotype != "None" && (file_name.find(haplotype) != std::string::npos))
			{
				write_final_svtigs_fasta(fasta_index, file_name, itr->second->pos, itr->second->contig, itr->second->reads.size(), fp_write);
				cnt++;
			}
			else if(haplotype == "None" && (file_name.find("H1") == std::string::npos) && (file_name.find("H2") == std::string::npos))
			{
				write_final_svtigs_fasta(fasta_index, file_name, itr->second->pos, itr->second->contig, itr->second->reads.size(), fp_write);
				cnt++;
			}

    		//params.fp_logs << file_name <<" contig="<<itr->second->contig<<" pos="<<itr->second->pos<<" support="<<itr->second->reads.size() << "\n";
			//for (auto r : itr->second->reads)
			//	params.fp_logs << r << std::endl;
  			//params.fp_logs << "\n";
		}
	}

	fp_write.close();
	return cnt;
}


int filter_svtigs(parameters& params, std::map<std::string, gfaNode*>& gfa, std::map <std::string, SVtig*>& final_svtigs)
{
	std::cout<<"\nFiltering svtigs"<<std::endl;	
	std::string svtigs_tmp_path = params.log_path + params.sample_name + "_svtigs_tmp.fa";

	std::cout<<"--->remapping svtigs onto the graph using Graphaligner"<<std::endl;			

	// GraphAligner PATH'te mi?
	std::string graphaligner_bin = find_executable("GraphAligner");
	if (graphaligner_bin.empty())
	{
		std::string msg = "[filter_svtigs] GraphAligner not found in PATH. "
		                  "Please install GraphAligner or add it to your PATH.";
		if (params.fp_logs.is_open())
			params.fp_logs << msg << std::endl;
		error(msg.c_str());
	}

	// Komutu GraphAligner binary'si ile kur, stdout/stderr sessiz
	std::string graphaligner_cmd =
		graphaligner_bin +
		" -g " + params.ref_graph +
		" -f " + svtigs_tmp_path +
		" -a " + params.remap_gaf_path +
		" -t " + std::to_string(params.threads) +
		" -x vg"
		" --precise-clipping " + std::to_string(params.min_precise_clipping) +
		" --min-alignment-score " + std::to_string(params.min_alignment_score) +
		" --multimap-score-fraction 0.9"
		" >/dev/null 2>&1";
		
	run_and_log(graphaligner_cmd, params, "GraphAligner", 2, 2, true);
		
	if (std::filesystem::is_empty(params.remap_gaf_path))
	{
		std::cout<< "Error: Graphaligner did not run successfully..."<< "\n";
		std::cout<< "--->Command: " << graphaligner_cmd << "\n";
		exit(0);
	}
	
	std::string file_path = params.log_path + params.sample_name + "_svtigs_tmp.fa";
	faidx_t* fasta_index = fai_load((file_path).c_str());
	if (!fasta_index)
		error("Error loading FASTA index for remapping: file not found or corrupted");

	//Now the ones that we want to output have final_svtigs->output = true	
	read_remappings(params, gfa, final_svtigs, fasta_index);	
	
	std::string svtigs_path;
	if ((params.phase_tags).empty())
	{
		svtigs_path = params.log_path + params.sample_name + "_svtigs.fa";
		int h1 = write_final_svtigs(fasta_index, final_svtigs, svtigs_path, "None");

		std::cout<<"--->there are "<<h1<<" svtigs"<<"\n";
	}
	else
	{
		svtigs_path = params.log_path + params.sample_name + "_svtigs_H1.fa";
		int h1 = write_final_svtigs(fasta_index, final_svtigs, svtigs_path, "H1");

		svtigs_path = params.log_path + params.sample_name + "_svtigs_H2.fa";
		int h2 = write_final_svtigs(fasta_index, final_svtigs, svtigs_path, "H2");
		
		if (!params.skip_untagged)	
		{
			svtigs_path = params.log_path + params.sample_name + "_svtigs_untagged.fa";
			int untagged = write_final_svtigs(fasta_index, final_svtigs, svtigs_path, "None");
			
			std::cout<<"--->there are "<<h1<<" haplotype 1, " <<h2<<" haplotype 2 and "<<untagged<<" untagged svtigs"<<"\n";
		}
		else
			std::cout<<"--->there are "<<h1<<" haplotype 1, " <<h2<<" haplotype 2 svtigs"<<"\n";
	}
	fai_destroy(fasta_index);

	if (!params.debug)
	{
		if(std::filesystem::exists(svtigs_tmp_path))
			std::filesystem::remove_all(svtigs_tmp_path);
		if(std::filesystem::exists(svtigs_tmp_path + ".fai"))
			std::filesystem::remove_all(svtigs_tmp_path + ".fai");
		if(std::filesystem::exists(params.remap_gaf_path))
			std::filesystem::remove_all(params.remap_gaf_path);
		if(std::filesystem::exists(params.log_path + "in/"))
			std::filesystem::remove_all(params.log_path + "in/");
		if(std::filesystem::exists(params.log_path + "out/"))
			std::filesystem::remove_all(params.log_path + "out/");
	}

	return RETURN_SUCCESS;
}
