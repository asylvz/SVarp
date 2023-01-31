#ifndef __REMAP
#define __REMAP

#include <map>
#include "reference.h"
#include "common.h"

int remap_assemblies(parameters* params);
int read_remappings(std::map<std::string, gfaNode*> gfa);

#endif
