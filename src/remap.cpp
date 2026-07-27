#include <iostream>
#include <filesystem>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
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


//Header of an output svtig. path/map_ratio come from remapping, so they are
//omitted when the svtig was never remapped (--no-remap).
std::string svtig_header(const SVtig* svtig)
{
	std::ostringstream out;
	out << ">" << svtig->name << " contig=" << svtig->contig << " pos=" << svtig->pos
		<< " support=" << svtig->reads.size();

	if (svtig->map_ratio >= 0)
		out << " path=" << svtig->remap_path << " map_ratio=" << svtig->map_ratio;

	return out.str();
}


int write_final_svtigs_fasta(faidx_t*& fasta_index, SVtig* svtig, std::ostream& fp_write)
{

	int loc_length;
	const size_t line_len = 60;

	char *tmp_seq = faidx_fetch_seq(fasta_index, (svtig->name).c_str(), 0, MAX_FETCH_LEN, &loc_length);
	if (tmp_seq == nullptr)
		return RETURN_ERROR;
	std::string seq(tmp_seq);
	free(tmp_seq);
	fp_write<< svtig_header(svtig) << std::endl;

	for (unsigned int i = 0; i < seq.size(); i += line_len)
	{
		std::string tmp = (seq.substr(i, line_len));
		fp_write<< tmp << std::endl;
    }
	return RETURN_SUCCESS;
}


// An svtig still carries an SV when its realignment to the graph diverges by at
// least MINSVSIZE, the same bound the CIGAR walk applies when the signal is
// first picked up from the reads.
bool cigar_has_sv(const std::string& cigar)
{
	std::vector<int> cigar_len;
	std::vector<char> cigar_op;

	int cigar_cnt = decompose_cigars(cigar, cigar_len, cigar_op);
	for (int c = 0; c < cigar_cnt; c++)
	{
		if ((cigar_op[c] == INSERTION || cigar_op[c] == DELETION || cigar_op[c] == MISMATCH)
		    && cigar_len[c] >= MINSVSIZE)
			return true;
	}
	return false;
}


// Svtigs assembled from different haplotypes describe different alleles of the
// same locus, so they are never duplicates of each other. Without a phase file
// the names carry no prefix and all svtigs fall into one class.
static std::string haplotype_of(const std::string& svtig_name)
{
	if (svtig_name.rfind("H1-", 0) == 0)
		return "H1";
	if (svtig_name.rfind("H2-", 0) == 0)
		return "H2";
	if (svtig_name.rfind("None-", 0) == 0)
		return "None";
	return "";
}


std::pair<int, int> remove_duplicates(std::vector <Read*>& tmp_svtig, std::map <std::string, SVtig*>& final_svtigs, int& extra_added)
{
	std::pair<int, int> dup_legit(0,0);
	std::map<std::string, SVtig*>::iterator it_svtigs;

	// Sort by haplotype and node first, then by svtig_size descending (greedy: largest first).
	// Size alone decides the representative; map ratio and sv_in_cigar are on the
	// Read as well and would rank candidates better.
	std::sort(tmp_svtig.begin(), tmp_svtig.end(), [](const Read* a, const Read* b) {
		const std::string hap_a = haplotype_of(a->rname), hap_b = haplotype_of(b->rname);
		if (hap_a != hap_b)
			return hap_a < hap_b;
		if (a->node != b->node)
			return a->node < b->node;
		return a->svtig_size > b->svtig_size; // descending by size
	});

	// Process each haplotype/node group independently
	unsigned int i = 0;
	while (i < tmp_svtig.size())
	{
		// Find the range [i, group_end) for the current haplotype and node.
		// Paths have to match exactly, so svtigs of one locus reached over
		// slightly different paths never compare.
		unsigned int group_end = i + 1;
		while (group_end < tmp_svtig.size() && tmp_svtig[group_end]->node == tmp_svtig[i]->node
		       && haplotype_of(tmp_svtig[group_end]->rname) == haplotype_of(tmp_svtig[i]->rname))
			group_end++;

		// Greedy: iterate from largest svtig_size to smallest
		// Each candidate is compared only against previously kept entries
		std::vector<unsigned int> kept_indices;

		for (unsigned int k = i; k < group_end; k++)
		{
			Read *r = tmp_svtig[k];
			if (r->duplicate)
				continue;

			bool is_dup = false;
			for (unsigned int ki : kept_indices)
			{
				double overlap = overlap_ratio(r->start, r->end, tmp_svtig[ki]->start, tmp_svtig[ki]->end);
				if (overlap > MIN_DUP_OVERLAP)
				{
					r->duplicate = true;
					dup_legit.first++;
					is_dup = true;
					break;
				}
			}

			if (!is_dup)
			{
				kept_indices.push_back(k);

				// Handle fragmented assemblies (multiple contigs from same svtig)
				auto npos = r->rname.rfind("_");
				auto npos2 = r->rname.find("_");
				if (npos != npos2)
				{
					std::string tmp_str = r->rname.substr(0, npos);
					it_svtigs = final_svtigs.find(tmp_str);

					if (it_svtigs != final_svtigs.end())
					{
						SVtig *tmp = new SVtig;
						tmp->name = r->rname;
						tmp->pos = it_svtigs->second->pos;
						(tmp->reads).insert(it_svtigs->second->reads.begin(), it_svtigs->second->reads.end());
						tmp->contig = it_svtigs->second->contig;
						tmp->output = true;
						tmp->remap_path = r->node;
						tmp->map_ratio = r->highest_map_ratio;
						final_svtigs.insert(std::pair<std::string, SVtig*>(r->rname, tmp));
						extra_added++;
					}
					else
						std::cerr<<"Error - SVtig= "<<r->rname<<" not found (multiple contig)...\n";
				}
				else
				{
					std::map<std::string, SVtig*>::iterator it_dup = final_svtigs.find(r->rname);
					if (it_dup != final_svtigs.end())
					{
						it_dup->second->output = true;
						it_dup->second->remap_path = r->node;
						it_dup->second->map_ratio = r->highest_map_ratio;
					}
					else
						std::cerr<<"Error - SVtig= "<<r->rname<<" not found...\n";
				}
				dup_legit.second++;
			}
		}

		i = group_end;
	}

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

		if (gfa.count(node_name) == 0)
			continue;

		if (strand == '>')
			ref_tmp += gfa[node_name]->sequence;
		else if (strand == '<')
		{
			// reverse_complement mutates in place; work on a copy so the
			// stored node sequence is not corrupted.
			std::string rc = gfa[node_name]->sequence;
			ref_tmp += reverse_complement(rc);
		}
		else
			std::cerr<<"Strand resolution issue in wfa_align()\n";
	}
	
	// A dropped path node can leave ref_tmp shorter than gfa_start; skip
	// rather than substr past the end (which would throw out_of_range).
	if (gfa_start < 0 || static_cast<size_t>(gfa_start) > ref_tmp.size())
		return;
	std::string ref = ref_tmp.substr(gfa_start, gfa_end - gfa_start + 1);

	int loc_length;
	char *tmp_query = faidx_fetch_seq(fasta_index, query_name.c_str(), query_start, query_end, &loc_length);
	if (tmp_query == nullptr)
		return;
	std::string query(tmp_query);
	free(tmp_query);

	aligner.alignEnd2End(query, ref); // Align
	cigar = aligner.getCIGAR(true);
}


int read_remappings(parameters& params, std::map<std::string, gfaNode*>& gfa, std::map <std::string, SVtig*>& final_svtigs, faidx_t*& fasta_index)
{
	int secondary = 0, primary = 0, lowmq = 0, extra_added = 0;

	std::string line;
	std::vector <std::string> tokens;
	std::ifstream fp;
	
	std::cout<<"--> reading remappings from "<< params.remap_gaf_path <<std::endl;
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
		

		std::string cigar;

		//Use the cigar of WFA realignment
		wfa_align(gfa, cigar, g.query_name, g.query_start, g.query_end, g.path, g.path_start, g.path_end, aligner, fasta_index);

		hasSV = cigar_has_sv(cigar);

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
			// freq counts the good alignments of this svtig, and later records
			// only add to it when they are primary and pass MINMAPQREMAP, so the
			// first record has to be weighed the same way.
			r->freq = (LowMQ || !g.is_primary) ? 0 : 1;
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
			else {
				filtered++;
				if (params.fp_remap_log.is_open())
					params.fp_remap_log << t.first << "\tFILTERED\treason=high_map_ratio\tmap_ratio=" << t.second->highest_map_ratio << "\tno_sv_in_cigar\n";
			}
		}
		else {
			filtered++;
			if (params.fp_remap_log.is_open())
				params.fp_remap_log << t.first << "\tFILTERED\treason=no_mapping\n";
		}
	}
	//If an svtig is in tmp_svtig, check if it is a duplicate. If not set "->output=true"
	//Check the unmapped svtigs and add them tooo
	std::pair<int, int> dup_legit = remove_duplicates(tmp_svtig, final_svtigs, extra_added);

	if (params.fp_remap_log.is_open()) {
		for (auto &r : tmp_svtig) {
			if (r->duplicate)
				params.fp_remap_log << r->rname << "\tDUPLICATE\tnode=" << r->node << "\tstart=" << r->start << "\tend=" << r->end << "\tsize=" << r->svtig_size << "\n";
			else
				params.fp_remap_log << r->rname << "\tKEPT\tnode=" << r->node << "\tmap_ratio=" << r->highest_map_ratio << "\tsv_in_cigar=" << (r->sv_in_cigar ? "yes" : "no") << "\tsize=" << r->svtig_size << "\n";
		}
	}

	std::cout<<"--> "<< filtered << " filtered - " << dup_legit.first<<" duplicate\n";
	// Detailed mapping stats only in log file

	if (params.fp_logs.is_open()) {
		params.fp_logs << "--> " << filtered << " filtered - " << dup_legit.first << " duplicate\n";
		params.fp_logs << "--> " << primary << " primary, " << secondary << " secondary mappings, " << lowmq << " low MAPQ(<" << MINMAPQREMAP << "); svtigs from multiple contig assemblies = " << extra_added << "\n";
	}

	for (auto& p : reads)
		delete p.second;

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
				write_final_svtigs_fasta(fasta_index, itr->second, fp_write);
				cnt++;
			}
			else if(haplotype == "None" && (file_name.find("H1") == std::string::npos) && (file_name.find("H2") == std::string::npos))
			{
				write_final_svtigs_fasta(fasta_index, itr->second, fp_write);
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


static void remap_and_flag(parameters& params, std::map<std::string, gfaNode*>& gfa,
                           std::map <std::string, SVtig*>& final_svtigs, faidx_t*& fasta_index,
                           std::string& svtigs_tmp_path)
{
	std::cout<<"--> remapping svtigs onto the graph using GraphAligner"<<std::endl;

	// Is GraphAligner on PATH?
	std::string graphaligner_bin = find_executable("GraphAligner");
	if (graphaligner_bin.empty())
	{
		std::string msg = "[filter_svtigs] GraphAligner not found in PATH. "
		                  "Please install GraphAligner or add it to your PATH.";
		if (params.fp_logs.is_open())
			params.fp_logs << msg << std::endl;
		error(msg.c_str());
	}

	std::string ga_log = params.log_path + params.sample_name + "_graphaligner.log";
	std::string ga_redir = params.debug ? (" >" + ga_log + " 2>&1") : " >/dev/null 2>&1";

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
		+ ga_redir;
		
	run_and_log(graphaligner_cmd, params, "GraphAligner", 2, 2, true);
		
	if (std::filesystem::is_empty(params.remap_gaf_path))
	{
		std::cerr<< "Error: Graphaligner did not run successfully..."<< "\n";
		std::cerr<< "--->Command: " << graphaligner_cmd << "\n";
		exit(0);
	}

	//Now the ones that we want to output have final_svtigs->output = true	
	read_remappings(params, gfa, final_svtigs, fasta_index);
}


int filter_svtigs(parameters& params, std::map<std::string, gfaNode*>& gfa, std::map <std::string, SVtig*>& final_svtigs)
{
	std::cout<<"\nFiltering svtigs"<<std::endl;
	std::string svtigs_tmp_path = params.log_path + params.sample_name + "_svtigs_tmp.fa";

	faidx_t* fasta_index = fai_load(svtigs_tmp_path.c_str());
	if (!fasta_index)
		error("Error loading FASTA index for remapping: file not found or corrupted");

	if (params.no_remap)
	{
		// write_final_svtigs only writes svtigs flagged by read_remappings,
		// so skipping the remap means flagging them all here.
		std::cout<<"--> skipping GraphAligner remapping (--no-remap)"<<std::endl;
		if (params.fp_logs.is_open())
			params.fp_logs << "--> skipping GraphAligner remapping (--no-remap)\n";

		std::map<std::string, SVtig*>::iterator it;
		for (it = final_svtigs.begin(); it != final_svtigs.end(); ++it)
			it->second->output = true;
	}
	else
		remap_and_flag(params, gfa, final_svtigs, fasta_index, svtigs_tmp_path);

	std::string svtigs_path;
	if ((params.phase_tags).empty())
	{
		svtigs_path = params.log_path + params.sample_name + "_svtigs.fa";
		int h1 = write_final_svtigs(fasta_index, final_svtigs, svtigs_path, "None");

		std::cout<<"--> "<<h1<<" svtigs after filtering\n";
		if (params.fp_logs.is_open())
			params.fp_logs << "--> " << h1 << " svtigs after filtering\n";
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
			
			std::cout<<"--> "<<h1<<" haplotype 1, " <<h2<<" haplotype 2 and "<<untagged<<" untagged svtigs after filtering\n";
			if (params.fp_logs.is_open())
				params.fp_logs << "--> " << h1 << " haplotype 1, " << h2 << " haplotype 2 and " << untagged << " untagged svtigs after filtering\n";
		}
		else
		{
			std::cout<<"--> "<<h1<<" haplotype 1, " <<h2<<" haplotype 2 svtigs after filtering\n";
			if (params.fp_logs.is_open())
				params.fp_logs << "--> " << h1 << " haplotype 1, " << h2 << " haplotype 2 svtigs after filtering\n";
		}
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
