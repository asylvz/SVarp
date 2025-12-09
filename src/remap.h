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

}Read;

int filter_svtigs(parameters& params, std::map<std::string, gfaNode*>& gfa, std::map <std::string, SVtig*>& final_svtigs);

#endif
