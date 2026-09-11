#ifndef __REMAP
#define __REMAP

#include <map>
#include <vector>
#include <utility>
#include <iosfwd>
#include "common.h"
#include "reference.h"

class SVtig;
struct faidx_t;

#define TRIMWINDOW 1000 //svtig window (bp) for the end-trimming identity check

struct IndelRun { int qs; int qe; int len; }; //merged indel and the query span it occupies


typedef struct _read
{
	std::string rname;
	std::string node;
	int start; //query_start
	int end; //query_end
	double highest_map_ratio;
	int svtig_size;
	double highest_aln_identity;
	int freq;
	bool sv_in_cigar = false;
	bool duplicate = false;
	bool has_good = false; //seen a primary record with MAPQ >= MINMAPQREMAP
	int span = 0; //query span of the representative record
	std::vector<std::pair<int, int>> ivals; //query intervals of records at or above min identity
	int max_indel = 0; //largest merged indel in any counted record
	double cov = 0; //fraction of the svtig covered by ivals
	int max_gap = 0; //largest uncovered stretch inside the covered range
	bool explained = false; //no SV-sized gap or indel against the graph
	int trim_start = 0; //kept part of the assembled sequence, query coordinates
	int trim_end = 0; //0 when nothing was trimmed
	bool trimmed = false;
	std::vector<long> win_match, win_aligned; //per TRIMWINDOW bp, over counted records
	std::vector<IndelRun> runs; //merged indels of counted records

}Read;

int filter_svtigs(parameters& params, std::map<std::string, gfaNode*>& gfa, std::map <std::string, SVtig*>& final_svtigs);
std::pair<int, int> remove_duplicates(std::vector<Read*>& tmp_svtig, std::map<std::string, SVtig*>& final_svtigs, int& extra_added);
int merged_indel(const std::string& cigar, int max_gap = 20);
void indel_runs(const std::string& cigar, int query_start, int max_gap, std::vector<IndelRun>& runs);
void window_identity(const std::string& cigar, int query_start, Read* r);
std::pair<int, int> trim_bounds(const Read* r, double min_identity);
void apply_trim(Read* r, int lo, int hi);
bool cigar_has_sv(const std::string& cigar);
void graph_fit(Read* r);
bool explained_by_graph(const Read* r);
void update_read(Read* r, const Gaf& g, bool good, bool has_sv, double map_ratio);
std::string collect_alt_nodes(const std::string& path, std::map<std::string, gfaNode*>& gfa);
bool reference_colinear(const std::string& path, std::map<std::string, gfaNode*>& gfa);
void fill_alt_nodes(std::map<std::string, SVtig*>& final_svtigs, std::map<std::string, gfaNode*>& gfa);
std::string svtig_header(const SVtig* svtig);
int write_final_svtigs_fasta(faidx_t*& fasta_index, SVtig* svtig, std::ostream& fp_write);
std::string haplotype_of(const std::string& svtig_name);
int write_final_svtigs(faidx_t*& fasta_index, std::map <std::string, SVtig*>& final_svtigs, std::string& out_file, std::string haplotype);

#endif
