#ifndef __ALIGNMENT
#define __ALIGNMENT

#include <vector>
#include <map>
#include "reference.h" 
#include "variant.h"
#include "common.h"


int read_alignments(parameters& params, std::map <std::string, Contig*>& ref, std::map<std::string, gfaNode*>& gfa, std::map<std::string, Variant*>& vars, std::set <std::string>& unmapped);

Gaf parse_gaf_line(std::string& line);

#endif
