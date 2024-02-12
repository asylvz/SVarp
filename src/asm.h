#ifndef __READ_ALIGNMENTS_ASM
#define __READ_ALIGNMENTS_ASM

#include <vector>
#include <map>
#include "reference.h" 
#include "variant.h"
#include "common.h"


int read_alignments_asm(parameters& params, std::map <std::string, Contig*>& ref, std::map<std::string, gfaNode*>& gfa, std::set <std::string>& unmapped);

#endif
