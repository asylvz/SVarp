#ifndef __PHASING
#define __PHASING

#include <map>
#include "common.h" 
#include "variant.h"

typedef struct _phase
{
	std::string read_name;
	std::string haplotype;
	std::string phase_set;
	std::string contig;

	_phase() {
    }
} phase;

void phase_svs(std::map<std::string, phase*> phased_reads, std::map<std::string, std::vector<SVCluster*>>& vars);
int read_phase_file(parameters& params, std::map<std::string, phase*>& phased_reads);

#endif
