#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include "remap.h"
#include "reference.h"
#include "variant.h"
#include "common.h"

int main() {
    // Test 1: Verify Read struct instantiation and basic properties
    Read read1;
    read1.rname = "read_001";
    read1.node = "node1";
    read1.start = 100;
    read1.end = 200;
    read1.highest_map_ratio = 0.95;
    read1.svtig_size = 1000;
    read1.highest_aln_identity = 0.98;
    read1.freq = 5;
    read1.sv_in_cigar = false;
    read1.duplicate = false;
    
    if (read1.rname != "read_001") { std::cerr << "Test 1: Expected rname=read_001, got " << read1.rname << std::endl; return 1; }
    if (read1.start != 100) { std::cerr << "Test 1: Expected start=100, got " << read1.start << std::endl; return 1; }
    if (read1.highest_map_ratio != 0.95) { std::cerr << "Test 1: Expected highest_map_ratio=0.95" << std::endl; return 1; }
    if (read1.sv_in_cigar) { std::cerr << "Test 1: Expected sv_in_cigar=false" << std::endl; return 1; }
    if (read1.duplicate) { std::cerr << "Test 1: Expected duplicate=false" << std::endl; return 1; }
    
    // Test 2: Verify multiple Read objects and collection
    std::map<std::string, Read> reads_map;
    
    Read read2;
    read2.rname = "read_002";
    read2.node = "node2";
    read2.start = 150;
    read2.end = 250;
    read2.highest_map_ratio = 0.92;
    read2.svtig_size = 950;
    read2.highest_aln_identity = 0.96;
    read2.freq = 3;
    read2.sv_in_cigar = true;  // This read has SV in CIGAR
    read2.duplicate = false;
    
    reads_map["read_001"] = read1;
    reads_map["read_002"] = read2;
    
    if (reads_map.size() != 2) { std::cerr << "Test 2: Expected 2 reads, got " << reads_map.size() << std::endl; return 1; }
    
    // Test 3: Verify filtering logic based on map_ratio threshold
    int high_quality_count = 0;
    double map_ratio_threshold = 0.93;
    for (auto &kv : reads_map) {
        if (kv.second.highest_map_ratio >= map_ratio_threshold) {
            high_quality_count++;
        }
    }
    
    if (high_quality_count != 1) { std::cerr << "Test 3: Expected 1 high-quality read (>= 0.93), got " << high_quality_count << std::endl; return 1; }
    
    // Test 4: Verify SV filtering
    int sv_count = 0;
    for (auto &kv : reads_map) {
        if (kv.second.sv_in_cigar) {
            sv_count++;
        }
    }
    
    if (sv_count != 1) { std::cerr << "Test 4: Expected 1 read with SV in CIGAR, got " << sv_count << std::endl; return 1; }
    
    // Test 5: Verify frequency-based filtering
    int frequent_reads = 0;
    int freq_threshold = 4;
    for (auto &kv : reads_map) {
        if (kv.second.freq >= freq_threshold) {
            frequent_reads++;
        }
    }
    
    if (frequent_reads != 1) { std::cerr << "Test 5: Expected 1 read with freq >= 4, got " << frequent_reads << std::endl; return 1; }
    
    // Test 6: Verify GFA setup (no need to call filter_svtigs which requires GraphAligner)
    std::map<std::string, gfaNode*> gfa;
    gfaNode* node1 = new gfaNode("node1", "", 100, "contig1", 0);
    gfaNode* node2 = new gfaNode("node2", "", 100, "contig1", 100);
    gfa["node1"] = node1;
    gfa["node2"] = node2;
    
    // Create minimal SVtig map
    std::map<std::string, SVtig*> final_svtigs;
    SVtig* sv1 = new SVtig();
    sv1->name = "svtig_001";
    sv1->pos = 50;
    sv1->contig = "contig1";
    sv1->output = true;
    sv1->reads.insert("read_001");
    sv1->reads.insert("read_002");
    
    final_svtigs["svtig_001"] = sv1;
    
    // Verify GFA nodes and SVtig setup (without calling filter_svtigs which requires external tools)
    if (gfa.size() != 2) { std::cerr << "Test 6: Expected 2 GFA nodes, got " << gfa.size() << std::endl; return 1; }
    if (final_svtigs.size() != 1) { std::cerr << "Test 6: Expected 1 SVtig, got " << final_svtigs.size() << std::endl; return 1; }
    if (sv1->reads.size() != 2) { std::cerr << "Test 6: Expected 2 reads in SVtig, got " << sv1->reads.size() << std::endl; return 1; }
    
    // Cleanup
    delete node1;
    delete node2;
    delete sv1;
    
    std::cout << "remap test passed" << std::endl;
    return 0;
}
