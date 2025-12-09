#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include "variant.h"
#include "reference.h"
#include "common.h"

int main() {
    parameters params;
    params.dist_threshold = 10;
    params.support = 0; // allow single-read clusters for test

    // Setup gfa map
    std::map<std::string, gfaNode*> gfa;
    gfaNode* node1 = new gfaNode("node1", "", 100, "contig1", 0);
    gfa.insert({"node1", node1});

    // Create variants map for merge_svs
    std::map<std::string, Variant*> vars;
    Variant *v1 = new Variant();
    v1->node = "node1";
    v1->pos_in_node = 10;
    v1->reads_untagged.insert("r1");

    Variant *v2 = new Variant();
    v2->node = "node1";
    v2->pos_in_node = 15;
    v2->reads_untagged.insert("r2");

    Variant *v3 = new Variant();
    v3->node = "node1";
    v3->pos_in_node = 50;
    v3->reads_untagged.insert("r3");

    vars.insert({"node1:10", v1});
    vars.insert({"node1:15", v2});
    vars.insert({"node1:50", v3});

    std::map<std::string, std::vector<SVCluster*>> final_svtigs;
    std::map<std::string, std::vector<std::string>> incoming, outgoing;
    int rc = merge_svs(params, gfa, vars, final_svtigs, incoming, outgoing);
    if (rc != RETURN_SUCCESS) { std::cerr << "merge_svs_within_node returned error" << std::endl; return 1; }

    if (final_svtigs.size() != 1) { std::cerr << "Expected 1 node key in final_svtigs, got " << final_svtigs.size() << std::endl; return 1; }
    auto &vclusters = final_svtigs["node1"];
    if (vclusters.size() != 2) { std::cerr << "Expected 2 clusters, got " << vclusters.size() << std::endl; return 1; }

    // Check merged reads in first cluster
    auto &reads0 = vclusters[0]->reads_untagged;
    if (reads0.find("r1") == reads0.end() || reads0.find("r2") == reads0.end()) { std::cerr << "First cluster did not include both r1 and r2" << std::endl; return 1; }

    // Check second cluster contains r3
    auto &reads1 = vclusters[1]->reads_untagged;
    if (reads1.find("r3") == reads1.end()) { std::cerr << "Second cluster did not include r3" << std::endl; return 1; }

    // Cleanup
    delete v1; delete v2; delete v3;
    for (auto &kv : final_svtigs)
        for (auto s : kv.second) delete s;
    delete node1;

    std::cout << "variant.merge_svs_within_node test passed" << std::endl;
    return 0;
}
