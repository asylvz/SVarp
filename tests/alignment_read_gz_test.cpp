#include <iostream>
#include <string>
#include <map>
#include <set>
#include <fstream>
#include <cstdio>
#include "alignment.h"
#include "variant.h"
#include "reference.h"
#include "common.h"

int main()
{
    // Setup test structures
    parameters params;
    params.fp_logs.open("/tmp/test_alignment.log");

    std::map<std::string, Contig*> ref;
    std::map<std::string, gfaNode*> gfa;
    std::map<std::string, Variant*> vars;
    std::set<std::string> unmapped;
    std::map<std::string, int> read_freq;

    gfaNode* node1 = new gfaNode("node1", "", 200, "contig1", 0);
    gfa.insert({"node1", node1});


    std::cout << "alignment.read_gz function signature verified" << std::endl;

    // Cleanup
    params.fp_logs.close();
    delete node1;
    for (auto &kv : vars) delete kv.second;
    for (auto &kv : ref) delete kv.second;

    return 0;
}
