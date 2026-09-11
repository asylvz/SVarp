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
#include "bindings/cpp/WFAligner.hpp"
#include "alignment.h"
#include "variant.h"
#include "reference.h"



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


//Non-reference nodes of the remap path as "first-last:contig" blocks; "" if a node is unknown
//True when the path walks reference nodes (SR 0) only, in reference order, without a skipped
//node or a change of direction: such a path spells the reference sequence itself. Nodes without
//an SR tag cannot be judged and give false.
bool reference_colinear(const std::string& path, std::map<std::string, gfaNode*>& gfa)
{
	if (path.empty() || (path[0] != '>' && path[0] != '<'))
		return false;
	const gfaNode* prev = nullptr;
	char dir = path[0];
	size_t p = 0;
	while (p < path.size())
	{
		char strand = path[p++];
		size_t q = p;
		while (q < path.size() && path[q] != '>' && path[q] != '<') ++q;
		auto it = gfa.find(path.substr(p, q - p));
		p = q;
		if (it == gfa.end() || it->second->rank != 0 || strand != dir)
			return false;
		const gfaNode* n = it->second;
		if (prev)
		{
			if (prev->contig != n->contig)
				return false;
			bool adjacent = (dir == '>') ? prev->offset + prev->len == n->offset : n->offset + n->len == prev->offset;
			if (!adjacent)
				return false;
		}
		prev = n;
	}
	return prev != nullptr;
}


std::string collect_alt_nodes(const std::string& path, std::map<std::string, gfaNode*>& gfa)
{
	//GAF also allows a stable sequence name ("chr21:100-900"), which has no nodes
	if (path.empty() || (path[0] != '>' && path[0] != '<'))
		return "";

	std::ostringstream out;
	std::string first, last, contig;
	bool have_block = false;

	auto flush = [&]() {
		if (!have_block)
			return;
		if (out.tellp() > 0)
			out << ';';
		out << first;
		if (last != first)
			out << '-' << last;
		out << ':' << contig;
		have_block = false;
	};

	size_t p = 0;
	while (p < path.size())
	{
		++p;   //skip '>' or '<'
		size_t q = p;
		while (q < path.size() && path[q] != '>' && path[q] != '<') ++q;
		std::string node_name = path.substr(p, q - p);
		p = q;

		if (node_name.empty())
			continue;

		auto it = gfa.find(node_name);
		if (it == gfa.end())
			return "";   //path and graph disagree

		const gfaNode* n = it->second;

		//rank 0 is the reference, -1 means no SR tag to judge by
		if (n->rank <= 0)
		{
			flush();
			continue;
		}

		if (have_block && n->contig == contig)
			last = node_name;
		else
		{
			flush();
			first = last = node_name;
			contig = n->contig;
			have_block = true;
		}
	}
	flush();

	return out.str();
}


//Fill alt_nodes for every svtig that survived remapping.
void fill_alt_nodes(std::map<std::string, SVtig*>& final_svtigs, std::map<std::string, gfaNode*>& gfa)
{
	for (auto& s : final_svtigs)
		if (s.second->output && !s.second->remap_path.empty())
			s.second->alt_nodes = collect_alt_nodes(s.second->remap_path, gfa);
}


//Output svtig header; remap fields are omitted with --no-remap
std::string svtig_header(const SVtig* svtig)
{
	std::ostringstream out;
	out << ">" << svtig->name << " contig=" << svtig->contig << " pos=" << svtig->pos
		<< " support=" << svtig->reads.size();

	if (svtig->map_ratio >= 0)
	{
		out << " path=" << svtig->remap_path << " graph_cov=" << svtig->map_ratio;
		if (svtig->graph_identity >= 0)
			out << " graph_identity=" << svtig->graph_identity;
		out << " max_gap=" << svtig->max_gap << " max_indel=" << svtig->max_indel
		    << " graph_explained=" << (svtig->graph_explained ? "yes" : "no");
		if (!svtig->alt_nodes.empty())
			out << " alt_nodes=" << svtig->alt_nodes;
		if (svtig->trim_end > 0)
			out << " trim=" << svtig->trim_start << "-" << svtig->trim_end;
	}

	return out.str();
}


int write_final_svtigs_fasta(faidx_t*& fasta_index, SVtig* svtig, std::ostream& fp_write)
{

	hts_pos_t loc_length;
	const size_t line_len = 60;

	hts_pos_t n = faidx_seq_len64(fasta_index, (svtig->name).c_str());
	if (n <= 0)
		return RETURN_ERROR;
	hts_pos_t s = 0, e = n - 1;
	if (svtig->trim_end > 0 && svtig->trim_end <= n) { s = svtig->trim_start; e = svtig->trim_end - 1; }
	char *tmp_seq = faidx_fetch_seq64(fasta_index, (svtig->name).c_str(), s, e, &loc_length);
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


// Largest indel in a CIGAR, merging runs of the same operation that are separated
// by at most max_gap matched bases (an aligner splits a long indel in repeats).
int merged_indel(const std::string& cigar, int max_gap)
{
	std::vector<int> len;
	std::vector<char> op;
	int n = decompose_cigars(cigar, len, op);
	int best = 0, cur = 0, gap = 0;
	char cur_op = 0;
	for (int i = 0; i < n; i++)
	{
		if (op[i] == INSERTION || op[i] == DELETION)
		{
			if (op[i] == cur_op && gap <= max_gap)
				cur += len[i];
			else
			{
				cur_op = op[i];
				cur = len[i];
			}
			gap = 0;
			if (cur > best) best = cur;
		}
		else
			gap += len[i];
	}
	return best;
}


bool cigar_has_sv(const std::string& cigar)
{
	return merged_indel(cigar) >= MINSVSIZE;
}


// Same merging as merged_indel, keeping the query span of every run so indels
// can be attributed to the kept part of a trimmed svtig.
void indel_runs(const std::string& cigar, int query_start, int max_gap, std::vector<IndelRun>& runs)
{
	std::vector<int> len;
	std::vector<char> op;
	int n = decompose_cigars(cigar, len, op);
	int q = query_start, gap = 0;
	char cur_op = 0;
	IndelRun cur{0, 0, 0};
	for (int i = 0; i < n; i++)
	{
		if (op[i] == INSERTION || op[i] == DELETION)
		{
			if (op[i] == cur_op && gap <= max_gap)
				cur.len += len[i];
			else
			{
				if (cur.len > 0) runs.push_back(cur);
				cur_op = op[i];
				cur = {q, q, len[i]};
			}
			gap = 0;
			if (op[i] == INSERTION) q += len[i];
			cur.qe = q;
		}
		else
		{
			gap += len[i];
			q += len[i];
		}
	}
	if (cur.len > 0) runs.push_back(cur);
}


// Identity of the svtig against the graph per TRIMWINDOW bp of its sequence,
// summed over the counted records.
void window_identity(const std::string& cigar, int query_start, Read* r)
{
	if (r->svtig_size <= 0) return;
	size_t nw = (r->svtig_size + TRIMWINDOW - 1) / TRIMWINDOW;
	if (r->win_match.size() != nw)
	{
		r->win_match.assign(nw, 0);
		r->win_aligned.assign(nw, 0);
	}
	std::vector<int> len;
	std::vector<char> op;
	int n = decompose_cigars(cigar, len, op);
	int q = query_start;
	for (int i = 0; i < n; i++)
	{
		if (q < 0 || q >= r->svtig_size) break;
		if (op[i] == DELETION)
		{
			r->win_aligned[q / TRIMWINDOW] += len[i];
			continue;
		}
		int left = len[i];
		while (left > 0 && q < r->svtig_size)
		{
			int w = q / TRIMWINDOW, step = std::min(left, (w + 1) * TRIMWINDOW - q);
			r->win_aligned[w] += step;
			if (op[i] != INSERTION && op[i] != MISMATCH) r->win_match[w] += step;
			q += step;
			left -= step;
		}
	}
}


// Kept part of the svtig: windows are dropped from both ends while fewer than
// half their bases align or the aligned bases fall below min_identity.
std::pair<int, int> trim_bounds(const Read* r, double min_identity)
{
	int nw = r->win_match.size();
	if (nw == 0) return {0, r->svtig_size};
	auto good = [&](int w) {
		long need = std::min(TRIMWINDOW, r->svtig_size - w * TRIMWINDOW);
		return r->win_aligned[w] * 2 >= need && r->win_match[w] >= min_identity * r->win_aligned[w];
	};
	int lo = 0;
	while (lo < nw && !good(lo)) lo++;
	int hi = nw - 1;
	while (hi >= lo && !good(hi)) hi--;
	if (lo > hi) return {0, 0};
	return {lo * TRIMWINDOW, std::min(r->svtig_size, (hi + 1) * TRIMWINDOW)};
}


// Restrict the read to [lo, hi): coverage intervals, size and indels follow.
void apply_trim(Read* r, int lo, int hi)
{
	bool trimmed = (lo > 0 || hi < r->svtig_size);
	r->trimmed = trimmed;
	r->trim_start = lo;
	r->trim_end = trimmed ? hi : 0;
	if (trimmed)
	{
		std::vector<std::pair<int, int>> kept;
		for (auto& iv : r->ivals)
		{
			int s = std::max(iv.first, lo), e = std::min(iv.second, hi);
			if (e > s) kept.push_back({s - lo, e - lo});
		}
		r->ivals = kept;
		r->svtig_size = std::max(0, hi - lo);
	}
	long m = 0, al = 0;
	for (int w = lo / TRIMWINDOW; w < (int) r->win_match.size() && w * TRIMWINDOW < hi; w++)
	{
		m += r->win_match[w];
		al += r->win_aligned[w];
	}
	r->identity = al > 0 ? static_cast<double>(m) / al : -1;
	r->max_indel = 0;
	for (auto& run : r->runs)
		if (run.qs < hi && std::max(run.qe, run.qs + 1) > lo && run.len > r->max_indel)
			r->max_indel = run.len;
}


// Coverage of the svtig by its counted graph alignments: union of query
// intervals, largest uncovered stretch between them, largest indel inside them.
void graph_fit(Read* r)
{
	r->cov = 0; r->max_gap = 0;
	if (r->ivals.empty() || r->svtig_size <= 0) return;
	std::sort(r->ivals.begin(), r->ivals.end());
	long covered = 0; int cs = r->ivals[0].first, ce = r->ivals[0].second;
	for (size_t i = 1; i < r->ivals.size(); i++)
	{
		int s = r->ivals[i].first, e = r->ivals[i].second;
		if (s <= ce) { if (e > ce) ce = e; }
		else
		{
			covered += ce - cs;
			if (s - ce > r->max_gap) r->max_gap = s - ce;
			cs = s; ce = e;
		}
	}
	covered += ce - cs;
	r->cov = static_cast<double>(covered) / r->svtig_size;
}


// The graph explains a kept svtig when its alignments leave no SV-sized hole and
// open no SV-sized indel; coverage itself is the output gate (min_graph_cov).
bool explained_by_graph(const Read* r)
{
	return r->max_gap < MINSVSIZE && r->max_indel < MINSVSIZE;
}


// Svtigs assembled from different haplotypes describe different alleles of the
// same locus, so they are never duplicates of each other. Without a phase file
// the names carry no prefix and all svtigs fall into one class.
std::string haplotype_of(const std::string& svtig_name)
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
						tmp->map_ratio = r->cov;
						tmp->max_gap = r->max_gap;
						tmp->max_indel = r->max_indel;
						tmp->graph_explained = r->explained;
						tmp->graph_identity = r->identity;
						tmp->trim_start = r->trim_start;
						tmp->trim_end = r->trim_end;
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
						it_dup->second->map_ratio = r->cov;
						it_dup->second->max_gap = r->max_gap;
						it_dup->second->max_indel = r->max_indel;
						it_dup->second->graph_explained = r->explained;
						it_dup->second->graph_identity = r->identity;
						it_dup->second->trim_start = r->trim_start;
						it_dup->second->trim_end = r->trim_end;
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
		{
			cigar.clear(); //unknown node: no reference to align against
			return;
		}

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


// Primary records with MAPQ >= MINMAPQREMAP decide the representative alignment
// (longest query span), the map ratio and the SV evidence. Low-MAPQ or secondary
// records only stand in while no such record exists.
void update_read(Read* r, const Gaf& g, bool good, bool has_sv, double map_ratio)
{
	int span = g.query_end - g.query_start;
	if (good)
	{
		if (!r->has_good)
		{
			r->has_good = true;
			r->span = 0;
			r->highest_map_ratio = 0;
			r->sv_in_cigar = false;
		}
		r->freq++;
		r->sv_in_cigar = r->sv_in_cigar || has_sv;
	}
	else if (r->has_good)
		return;

	if (map_ratio > r->highest_map_ratio)
		r->highest_map_ratio = map_ratio;
	if (span > r->span)
	{
		r->span = span;
		r->node = g.path;
		r->start = g.path_start;
		r->end = g.path_end;
	}
}


int read_remappings(parameters& params, std::map<std::string, gfaNode*>& gfa, std::map <std::string, SVtig*>& final_svtigs, faidx_t*& fasta_index)
{
	int secondary = 0, primary = 0, lowmq = 0, extra_added = 0;

	std::string line;
	std::ifstream fp;

	std::cout<<"--> reading remappings from "<< params.remap_gaf_path <<std::endl;
	fp.open(params.remap_gaf_path);
	if(!fp.good())
	{
		std::cerr << "Error opening '"<<params.remap_gaf_path << std::endl;
		return RETURN_ERROR;
	}

	std::map<std::string, Read*>::iterator it;
	std::map <std::string, Read*> reads;

	wfa::WFAlignerGapAffine aligner(4, 6, 2, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryMed);

	while(getline(fp, line))
	{
		Gaf g;
		if (parse_gaf_line(line, g) != RETURN_SUCCESS)
			continue;

		bool LowMQ = g.mapping_quality < MINMAPQREMAP;
		if (LowMQ) lowmq++;
		if (!g.is_primary) secondary++;
		else if (!LowMQ) primary++;

		// largest indel from the aligner's CIGAR and from a gap-affine realignment
		// of the same segment; either one may keep a split indel in one piece
		std::string wfa_cigar;
		wfa_align(gfa, wfa_cigar, g.query_name, g.query_start, g.query_end, g.path, g.path_start, g.path_end, aligner, fasta_index);
		int indel = std::max(merged_indel(g.cigar), merged_indel(wfa_cigar));
		double map_ratio = static_cast<double> ((double) g.query_end - g.query_start) / g.query_length;

		it = reads.find(g.query_name);
		Read* r;
		if (it != reads.end())
			r = it->second;
		else
		{
			r = new Read();
			r->rname = g.query_name;
			r->svtig_size = g.query_length;
			r->freq = 0;
			r->highest_map_ratio = 0;
			reads.insert(std::pair<std::string, Read*>(g.query_name, r));
		}
		update_read(r, g, g.is_primary && !LowMQ, indel >= MINSVSIZE, map_ratio);

		// every record at or above the identity floor counts as coverage, low MAPQ
		// included: a multimapping svtig is present in the graph, not absent
		if (g.identity < 0 || g.identity >= params.min_identity)
		{
			r->ivals.push_back(std::make_pair(g.query_start, g.query_end));
			indel_runs(g.cigar, g.query_start, 20, r->runs);
			indel_runs(wfa_cigar, g.query_start, 20, r->runs);
			window_identity(g.cigar, g.query_start, r);
		}
	}

	// Output gate: long enough and anchored in the graph; whether the graph also
	// holds the SV allele is recorded, not filtered.
	std::vector <Read*> tmp_svtig;
	int unaligned = 0, too_short = 0, low_cov = 0, reference = 0;

	auto trim_tag = [](const Read* r) { return r->trimmed ? "\ttrim=" + std::to_string(r->trim_start) + "-" + std::to_string(r->trim_start + r->svtig_size) : std::string(); };
	for (auto &t: reads)
	{
		Read* r = t.second;
		// noisy ends: windows aligning below trim_identity are cut before the gates
		std::pair<int, int> keep = params.trim_identity > 0 ? trim_bounds(r, params.trim_identity) : std::make_pair(0, r->svtig_size);
		apply_trim(r, keep.first, keep.second);
		graph_fit(r);
		if (r->svtig_size <= 0 || r->svtig_size < params.min_svtig_len)
		{
			too_short++;
			if (params.fp_remap_log.is_open())
				params.fp_remap_log << t.first << "\tFILTERED\treason=short\tsize=" << r->svtig_size << trim_tag(r) << "\n";
			continue;
		}
		if (r->cov < params.min_graph_cov)
		{
			low_cov++;
			if (params.fp_remap_log.is_open())
				params.fp_remap_log << t.first << "\tFILTERED\treason=low_graph_cov\tcov=" << r->cov << "\tmax_gap=" << r->max_gap << "\tmax_indel=" << r->max_indel << "\tsize=" << r->svtig_size << trim_tag(r) << "\n";
			continue;
		}
		r->explained = explained_by_graph(r);
		if (r->explained && !params.keep_reference && reference_colinear(r->node, gfa))
		{
			reference++;
			if (params.fp_remap_log.is_open())
				params.fp_remap_log << t.first << "\tFILTERED\treason=reference\tcov=" << r->cov << "\tnode=" << r->node << "\tsize=" << r->svtig_size << trim_tag(r) << "\n";
			continue;
		}
		tmp_svtig.push_back(r);
	}

	// svtigs without any alignment are dropped, but counted
	int nseq = faidx_nseq(fasta_index);
	for (int i = 0; i < nseq; i++)
	{
		const char* nm = faidx_iseq(fasta_index, i);
		if (nm && reads.find(nm) == reads.end())
		{
			unaligned++;
			if (params.fp_remap_log.is_open())
				params.fp_remap_log << nm << "\tFILTERED\treason=no_alignment\n";
		}
	}

	std::pair<int, int> dup_legit = remove_duplicates(tmp_svtig, final_svtigs, extra_added);
	fill_alt_nodes(final_svtigs, gfa);

	if (params.fp_remap_log.is_open()) {
		for (auto &r : tmp_svtig) {
			if (r->duplicate)
				params.fp_remap_log << r->rname << "\tDUPLICATE\tnode=" << r->node << "\tstart=" << r->start << "\tend=" << r->end << "\tsize=" << r->svtig_size << trim_tag(r) << "\n";
			else
				params.fp_remap_log << r->rname << "\tKEPT\tgraph_explained=" << (r->explained ? "yes" : "no") << "\tcov=" << r->cov << "\tidentity=" << r->identity << "\tmax_gap=" << r->max_gap << "\tmax_indel=" << r->max_indel << "\tnode=" << r->node << "\tsize=" << r->svtig_size << trim_tag(r) << "\n";
		}
	}

	int kept = 0, in_graph = 0;
	for (auto &r : tmp_svtig)
		if (!r->duplicate) { kept++; in_graph += r->explained; }

	std::cout << "--> " << kept << " svtigs kept (" << in_graph << " present in the graph, " << kept - in_graph << " not); filtered: " << too_short << " short, " << low_cov << " low graph coverage, " << reference << " reference-identical, " << unaligned << " without alignment, " << dup_legit.first << " duplicate\n";
	if (params.fp_logs.is_open()) {
		params.fp_logs << "--> " << kept << " svtigs kept (" << in_graph << " present in the graph, " << kept - in_graph << " not); filtered: " << too_short << " short, " << low_cov << " low graph coverage, " << reference << " reference-identical, " << unaligned << " without alignment, " << dup_legit.first << " duplicate\n";
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
			std::string hap = haplotype_of(file_name);
			if ((haplotype != "None" && hap == haplotype) || (haplotype == "None" && hap != "H1" && hap != "H2"))
				if (write_final_svtigs_fasta(fasta_index, itr->second, fp_write) == RETURN_SUCCESS)
					cnt++;

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
		" --min-alignment-score " + std::to_string(params.min_alignment_score) +
		" --multimap-score-fraction 0.9" +
		(params.min_precise_clipping > 0 ? " --precise-clipping " + std::to_string(params.min_precise_clipping) : std::string("")) +
		ga_redir;
		
	if (params.fp_logs.is_open())
		params.fp_logs << "--> GraphAligner " << tool_version(graphaligner_bin) << "\n--> " << graphaligner_cmd << "\n";
	run_and_log(graphaligner_cmd, params, "GraphAligner", 2, 2, true);
		
	if (std::filesystem::is_empty(params.remap_gaf_path))
	{
		std::cout << "[warning] GraphAligner aligned no svtig; every svtig is reported without an alignment" << std::endl;
		if (params.fp_logs.is_open())
			params.fp_logs << "[warning] GraphAligner aligned no svtig: " << graphaligner_cmd << "\n";
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

		// Flag what the assembly produced: clusters without a contig never reached
		// the FASTA, extra contigs of a cluster (<name>_N) exist only there.
		int nseq = faidx_nseq(fasta_index);
		for (int i = 0; i < nseq; i++)
		{
			std::string nm = faidx_iseq(fasta_index, i);
			auto it = final_svtigs.find(nm);
			if (it == final_svtigs.end())
			{
				auto us = nm.rfind('_');
				auto base = (us == std::string::npos) ? final_svtigs.end() : final_svtigs.find(nm.substr(0, us));
				if (base == final_svtigs.end())
					continue;
				SVtig* tmp = new SVtig;
				tmp->name = nm;
				tmp->pos = base->second->pos;
				tmp->reads = base->second->reads;
				tmp->contig = base->second->contig;
				it = final_svtigs.insert(std::pair<std::string, SVtig*>(nm, tmp)).first;
			}
			it->second->output = true;
		}
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

	if (!params.debug && !params.keep_remap)
	{
		if(std::filesystem::exists(svtigs_tmp_path))
			std::filesystem::remove_all(svtigs_tmp_path);
		if(std::filesystem::exists(svtigs_tmp_path + ".fai"))
			std::filesystem::remove_all(svtigs_tmp_path + ".fai");
		if(std::filesystem::exists(params.remap_gaf_path))
			std::filesystem::remove_all(params.remap_gaf_path);
	}
	if (!params.debug)
	{
		if(std::filesystem::exists(params.log_path + "in/"))
			std::filesystem::remove_all(params.log_path + "in/");
		if(std::filesystem::exists(params.log_path + "out/"))
			std::filesystem::remove_all(params.log_path + "out/");
	}

	return RETURN_SUCCESS;
}
