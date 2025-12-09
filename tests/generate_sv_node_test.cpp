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
    
    // Test 1: Single-node deletion at position 2, length 3
    Gaf line1;
    line1.query_name = "read1";
    line1.path = ">node1";
    line1.path_length = 12;
    line1.path_start = 0;
    line1.path_end = 12;
    
    Variant* v1 = generate_sv_node(gfa, line1, 2, 3, DELETION);
    if (v1 == nullptr) { std::cerr << "generate_sv_node returned nullptr for test 1" << std::endl; return 1; }
    if (v1->node != "node1") { std::cerr << "Test 1: Expected node1, got " << v1->node << std::endl; return 1; }
    if (v1->pos_in_node != 2) { std::cerr << "Test 1: Expected pos_in_node=2, got " << v1->pos_in_node << std::endl; return 1; }
    if (v1->sv_type != DELETION) { std::cerr << "Test 1: Expected sv_type=DELETION" << std::endl; return 1; }
    if (v1->pos_in_node_end != 5) { std::cerr << "Test 1: Expected pos_in_node_end=5, got " << v1->pos_in_node_end << std::endl; return 1; }
    delete v1;
    
    // Test 2: Single-node insertion at position 5, length 0
    Variant* v2 = generate_sv_node(gfa, line1, 5, 0, INSERTION);
    if (v2 == nullptr) { std::cerr << "generate_sv_node returned nullptr for test 2" << std::endl; return 1; }
    if (v2->node != "node1") { std::cerr << "Test 2: Expected node1, got " << v2->node << std::endl; return 1; }
    if (v2->pos_in_node != 5) { std::cerr << "Test 2: Expected pos_in_node=5, got " << v2->pos_in_node << std::endl; return 1; }
    if (v2->sv_type != INSERTION) { std::cerr << "Test 2: Expected sv_type=INSERTION" << std::endl; return 1; }
    if (v2->pos_in_node_end != 5) { std::cerr << "Test 2: Expected pos_in_node_end=5, got " << v2->pos_in_node_end << std::endl; return 1; }
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
    if (v3->pos_in_node != 1) { std::cerr << "Test 3: Expected pos_in_node=1, got " << v3->pos_in_node << std::endl; return 1; }
    if (v3->sv_type != DELETION) { std::cerr << "Test 3: Expected sv_type=DELETION" << std::endl; return 1; }
    delete v3;
    
    // Test 4: Multi-node path, variant spanning into second node
    // base_pos = 15 means CIGAR position 15, which is 12 (node1 len) + 3, so position 3 in node2
    Variant* v4 = generate_sv_node(gfa, line3, 15, 2, INSERTION);
    if (v4 == nullptr) { std::cerr << "generate_sv_node returned nullptr for test 4" << std::endl; return 1; }
    if (v4->node != "node2") { std::cerr << "Test 4: Expected node2, got " << v4->node << std::endl; return 1; }
    if (v4->pos_in_node != 3) { std::cerr << "Test 4: Expected pos_in_node=3, got " << v4->pos_in_node << std::endl; return 1; }
    if (v4->sv_type != INSERTION) { std::cerr << "Test 4: Expected sv_type=INSERTION" << std::endl; return 1; }
    delete v4;
    
    // Test 5: Multi-node, variant in last node
    // base_pos = 28 means CIGAR position 28, which is 12 (node1) + 12 (node2) + 4, so position 4 in node3
    Variant* v5 = generate_sv_node(gfa, line3, 28, 1, DELETION);
    if (v5 == nullptr) { std::cerr << "generate_sv_node returned nullptr for test 5" << std::endl; return 1; }
    if (v5->node != "node3") { std::cerr << "Test 5: Expected node3, got " << v5->node << std::endl; return 1; }
    if (v5->pos_in_node != 4) { std::cerr << "Test 5: Expected pos_in_node=4, got " << v5->pos_in_node << std::endl; return 1; }
    if (v5->sv_type != DELETION) { std::cerr << "Test 5: Expected sv_type=DELETION" << std::endl; return 1; }
    delete v5;
    
    // Cleanup
    delete node1;
    delete node2;
    delete node3;
    
    std::cout << "generate_sv_node test passed" << std::endl;
    return 0;
}
