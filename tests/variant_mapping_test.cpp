#include <iostream>
#include <string>
#include <map>
#include "variant.h"
#include "reference.h"
#include "common.h"

int main()
{
    // Prepare gfa nodes
    std::map<std::string, gfaNode*> gfa;
    gfaNode* n1 = new gfaNode("node1", "", 100, "contig1", 0);
    gfaNode* n2 = new gfaNode("node2", "", 100, "contig1", 100);
    gfa.insert({"node1", n1});
    gfa.insert({"node2", n2});

    // Prepare GAF line spanning two nodes
    Gaf line;
    line.query_name = "readX";
    line.query_length = 1000;
    line.query_start = 300; // ensure not skipped (MIN_READ_START_END_WINDOW == 200)
    line.query_end = 700;
    line.path = ">node1>node2"; // two nodes, both forward strand
    line.path_length = 200;
    line.path_start = 10;  // start within first node
    line.path_end = 160;   // end maps into second node (160 - first_node_len = 60)

    std::map<std::string, Variant*> variations_inter;

    int inserted = mapping_start_end(gfa, line, variations_inter);
    if (inserted != 2) { std::cerr << "Expected 2 inserted variants, got " << inserted << std::endl; return 1; }

    std::string k1 = "node1:10";
    std::string k2 = "node2:60";
    if (variations_inter.find(k1) == variations_inter.end()) { std::cerr << "Missing variant "<<k1<< std::endl; return 1; }
    if (variations_inter.find(k2) == variations_inter.end()) { std::cerr << "Missing variant "<<k2<< std::endl; return 1; }

    if (variations_inter[k1]->reads_untagged.find("readX") == variations_inter[k1]->reads_untagged.end()) { std::cerr << "read not recorded for "<<k1<<std::endl; return 1; }
    if (variations_inter[k2]->reads_untagged.find("readX") == variations_inter[k2]->reads_untagged.end()) { std::cerr << "read not recorded for "<<k2<<std::endl; return 1; }

    // Cleanup
    for (auto &kv : variations_inter) delete kv.second;
    delete n1; delete n2;

    std::cout << "variant.mapping_start_end test passed" << std::endl;
    return 0;
}
