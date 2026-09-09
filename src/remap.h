#ifndef __REMAP
#define __REMAP

#include <map>
#include "common.h"
#include "reference.h"

class SVtig;


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

}Read;

int filter_svtigs(parameters& params, std::map<std::string, gfaNode*>& gfa, std::map <std::string, SVtig*>& final_svtigs);
std::pair<int, int> remove_duplicates(std::vector<Read*>& tmp_svtig, std::map<std::string, SVtig*>& final_svtigs, int& extra_added);
bool cigar_has_sv(const std::string& cigar);
void update_read(Read* r, const Gaf& g, bool good, bool has_sv, double map_ratio);
std::string collect_alt_nodes(const std::string& path, std::map<std::string, gfaNode*>& gfa);
void fill_alt_nodes(std::map<std::string, SVtig*>& final_svtigs, std::map<std::string, gfaNode*>& gfa);
std::string svtig_header(const SVtig* svtig);

#endif
