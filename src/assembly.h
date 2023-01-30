#ifndef __ASSEMBLY
#define __ASSEMBLY

#include <map>
#include "common.h"
#include "sv.h"

const std::string REMAP_OUTPUT = "remap_output.gaf";
const std::string FASTA_OUTPUT = "merged_cns.fa";


void run_assembly(parameters* params, std::map<std::string, std::vector<svtig*>>& insertions);
int index_fasta(parameters* params, std::map<std::string, unsigned long>& fasta_index);
int remap_assemblies(parameters* params);
int merge_assemblies();

#endif
