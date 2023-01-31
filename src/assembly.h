#ifndef __ASSEMBLY
#define __ASSEMBLY

#include <map>
#include "common.h"
#include "sv.h"


void run_assembly(parameters* params, std::map<std::string, std::vector<svtig*>>& insertions);
int index_fasta(parameters* params, std::map<std::string, unsigned long>& fasta_index);
int remap_assemblies(parameters* params);
int merge_assemblies();

#endif
