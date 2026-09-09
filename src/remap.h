#ifndef __REMAP
#define __REMAP

#include <map>
#include <vector>
#include <utility>
#include "common.h"
#include "reference.h"

class SVtig;
struct faidx_t;


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

}Read;

int filter_svtigs(parameters& params, std::map<std::string, gfaNode*>& gfa, std::map <std::string, SVtig*>& final_svtigs);
std::pair<int, int> remove_duplicates(std::vector<Read*>& tmp_svtig, std::map<std::string, SVtig*>& final_svtigs, int& extra_added);
int merged_indel(const std::string& cigar, int max_gap = 20);
bool cigar_has_sv(const std::string& cigar);
void graph_fit(Read* r);
bool explained_by_graph(const Read* r, double min_cov);
void update_read(Read* r, const Gaf& g, bool good, bool has_sv, double map_ratio);
std::string collect_alt_nodes(const std::string& path, std::map<std::string, gfaNode*>& gfa);
void fill_alt_nodes(std::map<std::string, SVtig*>& final_svtigs, std::map<std::string, gfaNode*>& gfa);
std::string svtig_header(const SVtig* svtig);
std::string haplotype_of(const std::string& svtig_name);
int write_final_svtigs(faidx_t*& fasta_index, std::map <std::string, SVtig*>& final_svtigs, std::string& out_file, std::string haplotype);

#endif
