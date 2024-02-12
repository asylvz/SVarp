#ifndef __ASSEMBLY
#define __ASSEMBLY

#include <map>
#include "common.h"
#include "variant.h"


void run_assembly(parameters& params, std::map <std::string, Contig*>& depth, std::map<std::string, std::vector<Svtig*>>& vars, std::set <std::string>& unmapped, std::map <std::string, FinalSvtig*>& final_svtigs);


//int index_fasta(parameters* params, std::map<std::string, unsigned long>& fasta_index);
int remap_assemblies(parameters& params);

#endif
