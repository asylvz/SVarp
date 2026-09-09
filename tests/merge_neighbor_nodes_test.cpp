// merge_neighbor_nodes: clusters within dist_threshold across a link merge, measured from the
// cluster end that faces the link, for every orientation pair.
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include "variant.h"
#include "reference.h"
#include "common.h"

static SVCluster* make_cluster(const std::string& node, int start_pos, int end_pos, const std::string& tag, int n = 5)
{
    SVCluster* c = new SVCluster();
    c->node = node; c->contig = "chr1"; c->start_pos = start_pos; c->end_pos = end_pos; c->ref_pos = start_pos;
    for (int i = 0; i < n; i++) c->reads_untagged.insert(tag + std::to_string(i));
    return c;
}
static void cleanup(std::map<std::string, std::vector<SVCluster*>>& init, std::map<std::string, gfaNode*>& gfa)
{
    for (auto& kv : init) for (auto* c : kv.second) delete c;
    for (auto& kv : gfa) delete kv.second;
}

int main()
{
    parameters params;
    params.dist_threshold = 100;

    // Test 1: A+ -> B+. The gap is A's length minus the cluster END, not its start (N1):
    // a 95 bp cluster ending 5 bp before the link merges with a cluster 3 bp into B.
    {
        std::map<std::string, gfaNode*> gfa;
        gfa["A"] = new gfaNode("A", "", 1000, "chr1", 0);
        gfa["B"] = new gfaNode("B", "", 500, "chr1", 1000);
        EdgeMap in, out;
        in["B"] = { {"A", '+', '+'} };
        std::map<std::string, std::vector<SVCluster*>> init;
        init["A"] = { make_cluster("A", 100, 110, "far"), make_cluster("A", 900, 995, "a") };
        init["B"] = { make_cluster("B", 3, 3, "b"), make_cluster("B", 400, 400, "bfar") };

        merge_neighbor_nodes(params, gfa, init, in, out);
        if (!init["A"][1]->filter || init["A"][0]->filter) { std::cerr << "Test 1: only A's last cluster should merge" << std::endl; return 1; }
        if (init["B"][0]->reads_untagged.size() != 10 || init["B"][1]->reads_untagged.size() != 5) { std::cerr << "Test 1: reads should land in B's first cluster" << std::endl; return 1; }
        cleanup(init, gfa);
    }

    // Test 2: the same clusters, but the A cluster ends 150 bp before the link: no merge
    {
        std::map<std::string, gfaNode*> gfa;
        gfa["A"] = new gfaNode("A", "", 1000, "chr1", 0);
        gfa["B"] = new gfaNode("B", "", 500, "chr1", 1000);
        EdgeMap in, out;
        in["B"] = { {"A", '+', '+'} };
        std::map<std::string, std::vector<SVCluster*>> init;
        init["A"] = { make_cluster("A", 800, 850, "a") };
        init["B"] = { make_cluster("B", 3, 3, "b") };
        merge_neighbor_nodes(params, gfa, init, in, out);
        if (init["A"][0]->filter) { std::cerr << "Test 2: clusters 153 bp apart merged" << std::endl; return 1; }
        cleanup(init, gfa);
    }

    // Test 3: A- -> B+. A is left through its forward start, so A's FIRST cluster is the neighbour.
    {
        std::map<std::string, gfaNode*> gfa;
        gfa["A"] = new gfaNode("A", "", 1000, "chr1", 0);
        gfa["B"] = new gfaNode("B", "", 500, "chr1", 1000);
        EdgeMap in, out;
        in["B"] = { {"A", '-', '+'} };
        std::map<std::string, std::vector<SVCluster*>> init;
        init["A"] = { make_cluster("A", 4, 10, "a"), make_cluster("A", 990, 999, "aend") };
        init["B"] = { make_cluster("B", 3, 3, "b") };
        merge_neighbor_nodes(params, gfa, init, in, out);
        if (!init["A"][0]->filter || init["A"][1]->filter) { std::cerr << "Test 3: A's first cluster should merge across A- -> B+" << std::endl; return 1; }
        if (init["B"][0]->reads_untagged.size() != 10) { std::cerr << "Test 3: reads not merged" << std::endl; return 1; }
        cleanup(init, gfa);
    }

    // Test 4: A+ -> B-. B is entered through its forward end, so B's LAST cluster is the neighbour.
    {
        std::map<std::string, gfaNode*> gfa;
        gfa["A"] = new gfaNode("A", "", 1000, "chr1", 0);
        gfa["B"] = new gfaNode("B", "", 2000, "chr1", 1000);
        EdgeMap in, out;
        in["B"] = { {"A", '+', '-'} };
        std::map<std::string, std::vector<SVCluster*>> init;
        init["A"] = { make_cluster("A", 990, 998, "a") };
        init["B"] = { make_cluster("B", 5, 5, "bstart"), make_cluster("B", 1990, 1995, "b") };
        merge_neighbor_nodes(params, gfa, init, in, out);
        if (!init["A"][0]->filter) { std::cerr << "Test 4: A's last cluster should merge across A+ -> B-" << std::endl; return 1; }
        if (init["B"][1]->reads_untagged.size() != 10 || init["B"][0]->reads_untagged.size() != 5) { std::cerr << "Test 4: reads should land in B's last cluster" << std::endl; return 1; }
        cleanup(init, gfa);
    }

    // Test 5: a link from a node to itself is ignored
    {
        std::map<std::string, gfaNode*> gfa;
        gfa["A"] = new gfaNode("A", "", 50, "chr1", 0);
        EdgeMap in, out;
        in["A"] = { {"A", '+', '+'} };
        std::map<std::string, std::vector<SVCluster*>> init;
        init["A"] = { make_cluster("A", 2, 48, "a") };
        merge_neighbor_nodes(params, gfa, init, in, out);
        if (init["A"][0]->filter) { std::cerr << "Test 5: self link merged a cluster into itself" << std::endl; return 1; }
        cleanup(init, gfa);
    }

    // Test 6: a cluster already merged into a successor cannot absorb anything more.
    //   nX --> nY --> nA, nY holds a single cluster; nodes are visited in name order,
    //   so nA merges nY away before nY is reached, and nX must stay.
    {
        std::map<std::string, gfaNode*> gfa;
        gfa["nA"] = new gfaNode("nA", "", 1000, "chr1", 0);
        gfa["nX"] = new gfaNode("nX", "", 1000, "chr1", 2000);
        gfa["nY"] = new gfaNode("nY", "",  100, "chr1", 1000);
        EdgeMap in, out;
        in["nA"] = { {"nY", '+', '+'} };
        in["nY"] = { {"nX", '+', '+'} };
        std::map<std::string, std::vector<SVCluster*>> init;
        init["nA"] = { make_cluster("nA",  10,  10, "a") };
        init["nX"] = { make_cluster("nX", 990, 990, "x") };
        init["nY"] = { make_cluster("nY",  50,  50, "y") };
        merge_neighbor_nodes(params, gfa, init, in, out);
        if (!init["nY"][0]->filter) { std::cerr << "Test 6: nY cluster should merge into nA" << std::endl; return 1; }
        if (init["nX"][0]->filter) { std::cerr << "Test 6: nX merged into a dropped cluster" << std::endl; return 1; }
        cleanup(init, gfa);
    }

    std::cout << "merge_neighbor_nodes test passed" << std::endl;
    return 0;
}
