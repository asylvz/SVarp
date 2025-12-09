#ifndef __ALIGNMENT
#define __ALIGNMENT

#include <vector>
#include <map>
#include "common.h"
#include "reference.h"
// Forward declare Variant only - Contig and gfaNode are defined in reference.h
class Variant;


int read_alignments(parameters& params, std::map <std::string, Contig*>& ref, std::map<std::string, gfaNode*>& gfa, std::map<std::string, Variant*>& vars, std::set <std::string>& unmapped);
int read_gz(parameters& params, std::map <std::string, Contig*>& ref, std::map<std::string, gfaNode*>& gfa, std::map<std::string, Variant*>& vars, std::set <std::string>& unmapped, std::map <std::string, int>& read_freq);

#endif
