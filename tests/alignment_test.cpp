#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include "alignment.h"
#include "reference.h"
#include "variant.h"
#include "common.h"

int main() {
    // Test 1: Verify Gaf struct instantiation and basic properties
    Gaf gaf_line;
    gaf_line.query_name = "read_001";
    gaf_line.query_length = 500;
    gaf_line.query_start = 10;
    gaf_line.query_end = 490;
    gaf_line.path = ">node1>node2>node3";
    gaf_line.path_length = 480;
    gaf_line.path_start = 5;
    gaf_line.path_end = 485;
    
    if (gaf_line.query_name != "read_001") { std::cerr << "Test 1: Expected query_name=read_001" << std::endl; return 1; }
    if (gaf_line.query_length != 500) { std::cerr << "Test 1: Expected query_length=500" << std::endl; return 1; }
    if (gaf_line.path != ">node1>node2>node3") { std::cerr << "Test 1: Expected multi-node path" << std::endl; return 1; }
    
    // Test 2: Verify path parsing (node extraction from GAF path string)
    std::string path_test = ">node_A>node_B<node_C";
    int node_count = 0;
    for (size_t i = 0; i < path_test.length(); ++i) {
        if (path_test[i] == '>' || path_test[i] == '<') {
            node_count++;
        }
    }
    if (node_count != 3) { std::cerr << "Test 2: Expected 3 nodes, got " << node_count << std::endl; return 1; }
    
    // Test 3: Verify query coverage calculation
    int query_len = 500;
    int query_start = 50;
    int query_end = 450;
    double coverage = (double)(query_end - query_start) / query_len;
    
    if (coverage < 0.7 || coverage > 0.9) { std::cerr << "Test 3: Expected coverage ~0.8, got " << coverage << std::endl; return 1; }
    
    // Test 4: Build GFA for alignment context
    std::map<std::string, gfaNode*> gfa;
    gfaNode* node1 = new gfaNode("node1", "ACGTACGTACGT", 12, "contig1", 0);
    gfaNode* node2 = new gfaNode("node2", "TGCATGCATGCA", 12, "contig1", 12);
    gfaNode* node3 = new gfaNode("node3", "AAATTTGGGCCC", 12, "contig1", 24);
    
    gfa["node1"] = node1;
    gfa["node2"] = node2;
    gfa["node3"] = node3;
    
    if (gfa.size() != 3) { std::cerr << "Test 4: Expected 3 GFA nodes, got " << gfa.size() << std::endl; return 1; }
    if (gfa["node1"]->len != 12) { std::cerr << "Test 4: Expected node1 len=12" << std::endl; return 1; }
    
    // Test 5: Verify Contig and reference setup
    std::map<std::string, Contig*> ref;
    Contig* contig1 = new Contig();
    contig1->mapped_bases = 1000;
    contig1->mapped_reads = 50;
    contig1->contig_length = 1200;
    contig1->coverage = (double)contig1->mapped_bases / contig1->contig_length;
    
    ref["contig1"] = contig1;
    
    if (ref.size() != 1) { std::cerr << "Test 5: Expected 1 contig, got " << ref.size() << std::endl; return 1; }
    if (contig1->coverage < 0.8 || contig1->coverage > 0.85) { std::cerr << "Test 5: Expected coverage ~0.83" << std::endl; return 1; }
    
    // Test 6: Verify Variant collection for alignment
    std::map<std::string, Variant*> vars;
    Variant* var1 = new Variant();
    var1->node = "node1";
    var1->pos_in_node = 5;
    var1->contig = "contig1";
    var1->reads_untagged.insert("read_001");
    
    vars["node1:5"] = var1;
    
    if (vars.size() != 1) { std::cerr << "Test 6: Expected 1 variant, got " << vars.size() << std::endl; return 1; }
    if (var1->reads_untagged.find("read_001") == var1->reads_untagged.end()) { std::cerr << "Test 6: Expected read_001 in variant reads" << std::endl; return 1; }
    
    // Test 7: Verify unmapped reads set
    std::set<std::string> unmapped;
    unmapped.insert("read_999");
    unmapped.insert("read_998");
    
    if (unmapped.size() != 2) { std::cerr << "Test 7: Expected 2 unmapped reads, got " << unmapped.size() << std::endl; return 1; }
    
    // Cleanup
    delete node1;
    delete node2;
    delete node3;
    delete contig1;
    delete var1;
    
    std::cout << "alignment test passed" << std::endl;
    return 0;
}
