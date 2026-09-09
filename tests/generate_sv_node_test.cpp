#include <iostream>
#include <string>
#include <map>
#include "variant.h"
#include "reference.h"
#include "common.h"

int main() {
    // Test generate_sv_node with single-node path
    std::map<std::string, gfaNode*> gfa;
    
    // Setup GFA nodes
    gfaNode* node1 = new gfaNode("node1", "ACGTACGTACGT", 12, "contig1", 0);
    gfaNode* node2 = new gfaNode("node2", "TGCATGCATGCA", 12, "contig1", 12);
    gfaNode* node3 = new gfaNode("node3", "AAATTTGGGCCC", 12, "contig1", 24);
    
    gfa["node1"] = node1;
    gfa["node2"] = node2;
    gfa["node3"] = node3;
    
    // Positions are 0-based; base_pos arrives 1-based (add_variant passes base_pos + 1)
    // Test 1: Single-node deletion, first deleted base at offset 1, length 3
    Gaf line1;
    line1.query_name = "read1";
    line1.path = ">node1";
    line1.path_length = 12;
    line1.path_start = 0;
    line1.path_end = 12;
    
    Variant* v1 = generate_sv_node(gfa, line1, 2, 3, DELETION);
    if (v1 == nullptr) { std::cerr << "generate_sv_node returned nullptr for test 1" << std::endl; return 1; }
    if (v1->node != "node1") { std::cerr << "Test 1: Expected node1, got " << v1->node << std::endl; return 1; }
    if (v1->pos_in_node != 1) { std::cerr << "Test 1: Expected pos_in_node=1, got " << v1->pos_in_node << std::endl; return 1; }
    if (v1->sv_type != DELETION) { std::cerr << "Test 1: Expected sv_type=DELETION" << std::endl; return 1; }
    if (v1->pos_in_node_end != 4) { std::cerr << "Test 1: Expected pos_in_node_end=4, got " << v1->pos_in_node_end << std::endl; return 1; }
    delete v1;
    
    // Test 2: Single-node insertion before the base at offset 4
    Variant* v2 = generate_sv_node(gfa, line1, 5, 0, INSERTION);
    if (v2 == nullptr) { std::cerr << "generate_sv_node returned nullptr for test 2" << std::endl; return 1; }
    if (v2->node != "node1") { std::cerr << "Test 2: Expected node1, got " << v2->node << std::endl; return 1; }
    if (v2->pos_in_node != 4) { std::cerr << "Test 2: Expected pos_in_node=4, got " << v2->pos_in_node << std::endl; return 1; }
    if (v2->sv_type != INSERTION) { std::cerr << "Test 2: Expected sv_type=INSERTION" << std::endl; return 1; }
    if (v2->pos_in_node_end != 4) { std::cerr << "Test 2: Expected pos_in_node_end=4, got " << v2->pos_in_node_end << std::endl; return 1; }
    delete v2;
    
    // Test 3: Multi-node path, deletion in first node
    Gaf line3;
    line3.query_name = "read3";
    line3.path = ">node1>node2>node3";
    line3.path_length = 36;
    line3.path_start = 0;
    line3.path_end = 36;
    
    Variant* v3 = generate_sv_node(gfa, line3, 1, 2, DELETION);
    if (v3 == nullptr) { std::cerr << "generate_sv_node returned nullptr for test 3" << std::endl; return 1; }
    if (v3->node != "node1") { std::cerr << "Test 3: Expected node1, got " << v3->node << std::endl; return 1; }
    if (v3->pos_in_node != 0) { std::cerr << "Test 3: Expected pos_in_node=0, got " << v3->pos_in_node << std::endl; return 1; }
    if (v3->sv_type != DELETION) { std::cerr << "Test 3: Expected sv_type=DELETION" << std::endl; return 1; }
    delete v3;
    
    // Test 4: Multi-node path, variant spanning into second node
    // base_pos = 15: path offset 14 = 12 (node1) + 2, so offset 2 in node2
    Variant* v4 = generate_sv_node(gfa, line3, 15, 2, INSERTION);
    if (v4 == nullptr) { std::cerr << "generate_sv_node returned nullptr for test 4" << std::endl; return 1; }
    if (v4->node != "node2") { std::cerr << "Test 4: Expected node2, got " << v4->node << std::endl; return 1; }
    if (v4->pos_in_node != 2) { std::cerr << "Test 4: Expected pos_in_node=2, got " << v4->pos_in_node << std::endl; return 1; }
    if (v4->sv_type != INSERTION) { std::cerr << "Test 4: Expected sv_type=INSERTION" << std::endl; return 1; }
    delete v4;
    
    // Test 5: Multi-node, variant in last node
    // base_pos = 28: path offset 27 = 24 (node1 + node2) + 3, so offset 3 in node3
    Variant* v5 = generate_sv_node(gfa, line3, 28, 1, DELETION);
    if (v5 == nullptr) { std::cerr << "generate_sv_node returned nullptr for test 5" << std::endl; return 1; }
    if (v5->node != "node3") { std::cerr << "Test 5: Expected node3, got " << v5->node << std::endl; return 1; }
    if (v5->pos_in_node != 3) { std::cerr << "Test 5: Expected pos_in_node=3, got " << v5->pos_in_node << std::endl; return 1; }
    if (v5->sv_type != DELETION) { std::cerr << "Test 5: Expected sv_type=DELETION" << std::endl; return 1; }
    delete v5;
    
    // Test 6: A deletion starting on the last base of node1 (offset 11) and running into node2 stays on node1.
    Variant* v6 = generate_sv_node(gfa, line3, 11 + 1, 2, DELETION);
    if (v6 == nullptr) { std::cerr << "generate_sv_node returned nullptr for test 6" << std::endl; return 1; }
    if (v6->node != "node1") { std::cerr << "Test 6: Expected node1, got " << v6->node << std::endl; return 1; }
    if (v6->pos_in_node != 11) { std::cerr << "Test 6: Expected pos_in_node=11, got " << v6->pos_in_node << std::endl; return 1; }
    if (v6->pos_in_node_end != 13) { std::cerr << "Test 6: Expected pos_in_node_end=13, got " << v6->pos_in_node_end << std::endl; return 1; }
    delete v6;

    // Test 7: same event from the other strand (<node3<node2<node1, path offset 36 - 13 = 23)
    Gaf line7 = line3;
    line7.path = "<node3<node2<node1";
    Variant* v7 = generate_sv_node(gfa, line7, 23 + 1, 2, DELETION);
    if (v7 == nullptr) { std::cerr << "generate_sv_node returned nullptr for test 7" << std::endl; return 1; }
    if (v7->node != "node1" || v7->pos_in_node != 11 || v7->pos_in_node_end != 13) {
        std::cerr << "Test 7: Expected node1:11-13, got " << v7->node << ":" << v7->pos_in_node << "-" << v7->pos_in_node_end << std::endl; return 1;
    }
    delete v7;

    // Test 8: Reverse-strand insertion before forward base 4 of node1 (path offset 36 - 4 = 32)
    Variant* v8 = generate_sv_node(gfa, line7, 32 + 1, 0, INSERTION);
    if (v8 == nullptr) { std::cerr << "generate_sv_node returned nullptr for test 8" << std::endl; return 1; }
    if (v8->node != "node1" || v8->pos_in_node != 4) {
        std::cerr << "Test 8: Expected node1:4, got " << v8->node << ":" << v8->pos_in_node << std::endl; return 1;
    }
    delete v8;

    // Test 9: A node missing from the graph gives nullptr
    Gaf line9 = line3;
    line9.path = ">node1>nodeX>node3";
    if (generate_sv_node(gfa, line9, 15, 2, DELETION) != nullptr) { std::cerr << "Test 9: Expected nullptr for a missing node" << std::endl; return 1; }

    // Cleanup
    delete node1;
    delete node2;
    delete node3;
    
    std::cout << "generate_sv_node test passed" << std::endl;
    return 0;
}
