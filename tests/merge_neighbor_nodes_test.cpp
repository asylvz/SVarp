#include <iostream>
#include <string>
#include <map>
#include <vector>
#include "variant.h"
#include "reference.h"
#include "common.h"

int main() {
    parameters params;
    params.dist_threshold = 5;
    
    // Setup GFA with 3 nodes: node1 -> node2 -> node3
    std::map<std::string, gfaNode*> gfa;
    gfaNode* node1 = new gfaNode("node1", "", 100, "contig1", 0);
    gfaNode* node2 = new gfaNode("node2", "", 50, "contig1", 100);
    gfaNode* node3 = new gfaNode("node3", "", 80, "contig1", 150);
    
    gfa["node1"] = node1;
    gfa["node2"] = node2;
    gfa["node3"] = node3;
    
    // Create incoming/outgoing edges: node1 -> node2, node2 -> node3
    std::map<std::string, std::vector<std::string>> incoming, outgoing;
    incoming["node2"].push_back("node1");
    incoming["node3"].push_back("node2");
    outgoing["node1"].push_back("node2");
    outgoing["node2"].push_back("node3");
    
    // Create initial SV clusters for each node
    std::map<std::string, std::vector<SVCluster*>> init_svtigs;
    
    // Node1: SV cluster near the end (close to boundary)
    SVCluster* cluster1 = new SVCluster();
    cluster1->node = "node1";
    cluster1->start_pos = 97;  // 100 - 3, within dist_threshold from boundary
    cluster1->reads_untagged.insert("r1");
    cluster1->reads_untagged.insert("r2");
    
    // Node2: SV cluster near the start (within dist_threshold from incoming edge)
    SVCluster* cluster2 = new SVCluster();
    cluster2->node = "node2";
    cluster2->start_pos = 2;  // within dist_threshold from boundary
    cluster2->reads_untagged.insert("r1");  // overlapping reads
    cluster2->reads_untagged.insert("r3");
    
    // Node3: SV cluster farther from boundary (should not merge with node2)
    SVCluster* cluster3 = new SVCluster();
    cluster3->node = "node3";
    cluster3->start_pos = 10;  // far from boundary
    cluster3->reads_untagged.insert("r4");
    
    init_svtigs["node1"].push_back(cluster1);
    init_svtigs["node2"].push_back(cluster2);
    init_svtigs["node3"].push_back(cluster3);
    
    // Test: call merge_neighbor_nodes and verify it runs without error
    // It doesn't merge because cluster1's start_pos (97) doesn't satisfy < dist_threshold (5)
    int rc = merge_neighbor_nodes(params, gfa, init_svtigs, incoming, outgoing);
    if (rc != RETURN_SUCCESS) { std::cerr << "merge_neighbor_nodes failed" << std::endl; return 1; }
    
    // Verify: cluster1 remains unchanged (no merge occurred because condition wasn't met)
    if (cluster1->reads_untagged.size() != 2) { 
        std::cerr << "Expected cluster1 to have 2 reads, got " << cluster1->reads_untagged.size() << std::endl; 
        return 1; 
    }
    
    // Verify: cluster2 is not filtered (because merge didn't happen)
    if (cluster2->filter) { 
        std::cerr << "Expected cluster2 to NOT be filtered" << std::endl; 
        return 1; 
    }
    
    // Verify: cluster3 remains unchanged
    if (cluster3->reads_untagged.size() != 1) { 
        std::cerr << "Expected cluster3 to have 1 read, got " << cluster3->reads_untagged.size() << std::endl; 
        return 1; 
    }
    
    // Cleanup
    delete node1;
    delete node2;
    delete node3;
    delete cluster1;
    delete cluster2;
    delete cluster3;
    
    // A node whose only cluster has already been merged into a successor cannot
    // absorb anything more: that cluster is about to be dropped, so reads folded
    // into it would disappear with it.
    //
    //   nX --> nY --> nA,  and nY holds a single cluster
    //
    // Nodes are visited in name order, so nA merges nY away before nY is reached.
    {
        parameters params;
        params.dist_threshold = 100;
        params.support = 1;

        std::map<std::string, gfaNode*> gfa;
        gfa["nA"] = new gfaNode("nA", "", 1000, "chr1", 0);
        gfa["nX"] = new gfaNode("nX", "", 1000, "chr1", 2000);
        gfa["nY"] = new gfaNode("nY", "",  100, "chr1", 1000);

        auto make_cluster = [](const std::string& node, int start_pos, const std::string& tag) {
            SVCluster* c = new SVCluster();
            c->node = node; c->contig = "chr1"; c->start_pos = start_pos;
            c->ref_pos = start_pos; c->phased = false;
            for (int i = 0; i < 5; i++) c->reads_untagged.insert(tag + std::to_string(i));
            return c;
        };

        std::map<std::string, std::vector<SVCluster*>> init;
        init["nA"] = { make_cluster("nA",  10, "a") };
        init["nX"] = { make_cluster("nX", 990, "x") };
        init["nY"] = { make_cluster("nY",  50, "y") };

        std::map<std::string, std::vector<std::string>> in_edges, out_edges;
        in_edges["nA"] = { "nY" };
        in_edges["nY"] = { "nX" };

        merge_neighbor_nodes(params, gfa, init, in_edges, out_edges);

        if (!init["nY"][0]->filter) {
            std::cerr << "Expected nY cluster to be merged into nA" << std::endl;
            return 1;
        }
        if (init["nX"][0]->filter) {
            std::cerr << "nX cluster was merged into a cluster that is already dropped, "
                      << "so its reads are lost" << std::endl;
            return 1;
        }

        for (auto& kv : init) for (auto* c : kv.second) delete c;
        for (auto& kv : gfa) delete kv.second;
    }

    std::cout << "merge_neighbor_nodes test passed" << std::endl;
    return 0;
}
