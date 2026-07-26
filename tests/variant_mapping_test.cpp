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

    // --- a boundary node absent from the graph must not crash (null deref) ---
    {
        std::map<std::string, gfaNode*> gfa2;
        gfaNode* g2 = new gfaNode("node2", "", 100, "contig1", 100);
        gfa2.insert({"node2", g2});          // node1 intentionally absent

        Gaf line2;
        line2.query_name = "readY";
        line2.query_length = 1000;
        line2.query_start = 300;             // skip_start = false
        line2.query_end = 700;               // skip_end = false
        line2.path = ">node1>node2";         // first node NOT in the graph
        line2.path_length = 200;
        line2.path_start = 10;
        line2.path_end = 160;

        std::map<std::string, Variant*> vi2;
        mapping_start_end(gfa2, line2, vi2); // must return, not segfault

        for (auto &kv : vi2)
            if (kv.first.empty() || kv.first[0] == ':') {
                std::cerr << "variant created for unresolved node: " << kv.first << std::endl;
                return 1;
            }
        for (auto &kv : vi2) delete kv.second;
        delete g2;
    }

    // --- a reverse-strand breakpoint must not depend on how many nodes follow ---
    {
        std::map<std::string, gfaNode*> gfa3;
        gfaNode* a = new gfaNode("nodeA", "", 20000, "chr1", 1000000);
        gfaNode* b = new gfaNode("nodeB", "", 20000, "chr1", 1020000);
        gfa3.insert({"nodeA", a});
        gfa3.insert({"nodeB", b});

        // 5'-clipped read, so only the start of the alignment is a breakpoint.
        Gaf single;
        single.query_name = "readR";
        single.query_length = 20000;
        single.query_start = 8000;   // skip_start = false
        single.query_end = 20000;    // skip_end = true
        single.path = "<nodeA";
        single.path_length = 20000;
        single.path_start = 2000;
        single.path_end = 14000;

        Gaf multi = single;
        multi.path = "<nodeA<nodeB";
        multi.path_length = 40000;

        std::map<std::string, Variant*> vi_single, vi_multi;
        mapping_start_end(gfa3, single, vi_single);
        mapping_start_end(gfa3, multi, vi_multi);

        // The alignment starts at path offset 2000 of a reverse-complemented
        // node, which is 20000 - 2000 within the node itself.
        const std::string expected = "nodeA:18000";
        if (vi_single.find(expected) == vi_single.end()) {
            std::cerr << "reverse single-node breakpoint wrong; expected " << expected << ", got";
            for (auto &kv : vi_single) std::cerr << " " << kv.first;
            std::cerr << std::endl;
            return 1;
        }
        if (vi_multi.find(expected) == vi_multi.end()) {
            std::cerr << "reverse multi-node breakpoint wrong; expected " << expected << ", got";
            for (auto &kv : vi_multi) std::cerr << " " << kv.first;
            std::cerr << std::endl;
            return 1;
        }

        for (auto &kv : vi_single) delete kv.second;
        for (auto &kv : vi_multi) delete kv.second;
        delete a; delete b;
    }

    // Cleanup
    for (auto &kv : variations_inter) delete kv.second;
    delete n1; delete n2;

    std::cout << "variant.mapping_start_end test passed" << std::endl;
    return 0;
}
