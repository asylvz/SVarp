#ifndef __PHASING
#define __PHASING

#include <map>
#include "common.h" 
#include "sv.h"

typedef struct _phase
{
	std::string read_name;
	std::string haplotype;
	std::string phase_set;
	std::string contig;

	_phase() {
    }
} phase;

void phase_svs(std::map<std::string, phase*> phased_reads, std::map<std::string, variant*>& insertions);
int read_phase_file(parameters *params, std::map<std::string, phase*>& phased_reads);

#endif
