#ifndef __REMAP
#define __REMAP

#include <map>
#include "reference.h"
#include "common.h"
#include "variant.h"


typedef struct _read
{
	std::string rname;
	std::string node;
	int start;
	int end;
	double highest_map_ratio;
	int svtig_size;
	double highest_aln_identity;
	int freq;
	bool sv_in_cigar = false;
	bool duplicate = false;

}Read;

int filter_svtigs(parameters& params, std::map<std::string, gfaNode*>& gfa, std::map <std::string, FinalSvtig*>& final_svtigs);

#endif
