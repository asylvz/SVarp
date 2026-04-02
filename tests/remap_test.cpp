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

    std::cout << "remap basic tests passed" << std::endl;

    // ================================================================
    // remove_duplicates tests
    // ================================================================
    auto make_read = [](const std::string& name, const std::string& node,
                        int start, int end, int svtig_size) -> Read* {
        Read* r = new Read;
        r->rname = name;
        r->node = node;
        r->start = start;
        r->end = end;
        r->svtig_size = svtig_size;
        r->highest_map_ratio = 0.95;
        r->highest_aln_identity = 0.95;
        r->freq = 1;
        r->sv_in_cigar = true;
        r->duplicate = false;
        return r;
    };

    auto make_svtig = [](const std::string& name) -> SVtig* {
        SVtig* s = new SVtig;
        s->name = name;
        s->pos = 0;
        s->contig = "chr1";
        s->output = false;
        return s;
    };

    // Test 7: No duplicates - different regions on same node
    {
        std::vector<Read*> reads;
        std::map<std::string, SVtig*> svtigs;
        int extra = 0;

        reads.push_back(make_read("H1-s100_1", "nodeA", 0, 1000, 500));
        reads.push_back(make_read("H1-s100_2", "nodeA", 5000, 6000, 400));
        reads.push_back(make_read("H1-s100_3", "nodeB", 0, 1000, 300));

        svtigs["H1-s100_1"] = make_svtig("H1-s100_1");
        svtigs["H1-s100_2"] = make_svtig("H1-s100_2");
        svtigs["H1-s100_3"] = make_svtig("H1-s100_3");

        auto result = remove_duplicates(reads, svtigs, extra);
        if (result.first != 0 || result.second != 3) {
            std::cerr << "Test 7 FAILED: expected 0 dup, 3 kept, got " << result.first << " dup, " << result.second << " kept" << std::endl;
            return 1;
        }
        std::cout << "Test 7 passed: no duplicates" << std::endl;
        for (auto* r : reads) delete r;
        for (auto& p : svtigs) delete p.second;
    }

    // Test 8: Two overlapping reads, larger svtig_size wins
    {
        std::vector<Read*> reads;
        std::map<std::string, SVtig*> svtigs;
        int extra = 0;

        reads.push_back(make_read("H1-s200_1", "nodeA", 100, 1100, 800));
        reads.push_back(make_read("H1-s200_2", "nodeA", 100, 1100, 400));

        svtigs["H1-s200_1"] = make_svtig("H1-s200_1");
        svtigs["H1-s200_2"] = make_svtig("H1-s200_2");

        auto result = remove_duplicates(reads, svtigs, extra);
        if (result.first != 1 || result.second != 1) {
            std::cerr << "Test 8 FAILED: expected 1 dup, 1 kept, got " << result.first << " dup, " << result.second << " kept" << std::endl;
            return 1;
        }
        bool larger_kept = false, smaller_dup = false;
        for (auto* r : reads) {
            if (r->svtig_size == 800 && !r->duplicate) larger_kept = true;
            if (r->svtig_size == 400 && r->duplicate) smaller_dup = true;
        }
        if (!larger_kept || !smaller_dup) {
            std::cerr << "Test 8 FAILED: larger svtig should be kept" << std::endl;
            return 1;
        }
        std::cout << "Test 8 passed: larger svtig wins" << std::endl;
        for (auto* r : reads) delete r;
        for (auto& p : svtigs) delete p.second;
    }

    // Test 9: Transitivity - greedy largest-first resolves consistently
    {
        std::vector<Read*> reads;
        std::map<std::string, SVtig*> svtigs;
        int extra = 0;

        reads.push_back(make_read("H1-s300_1", "nodeA", 0, 950, 500));
        reads.push_back(make_read("H1-s300_2", "nodeA", 0, 1000, 900));
        reads.push_back(make_read("H1-s300_3", "nodeA", 50, 1000, 400));

        svtigs["H1-s300_1"] = make_svtig("H1-s300_1");
        svtigs["H1-s300_2"] = make_svtig("H1-s300_2");
        svtigs["H1-s300_3"] = make_svtig("H1-s300_3");

        auto result = remove_duplicates(reads, svtigs, extra);
        if (result.first != 2 || result.second != 1) {
            std::cerr << "Test 9 FAILED: expected 2 dup, 1 kept, got " << result.first << " dup, " << result.second << " kept" << std::endl;
            return 1;
        }
        for (auto* r : reads) {
            if (r->svtig_size == 900 && r->duplicate) {
                std::cerr << "Test 9 FAILED: largest svtig should not be duplicate" << std::endl;
                return 1;
            }
        }
        std::cout << "Test 9 passed: transitivity" << std::endl;
        for (auto* r : reads) delete r;
        for (auto& p : svtigs) delete p.second;
    }

    // Test 10: Different nodes never duplicate each other
    {
        std::vector<Read*> reads;
        std::map<std::string, SVtig*> svtigs;
        int extra = 0;

        reads.push_back(make_read("H1-s400_1", "nodeA", 100, 1100, 500));
        reads.push_back(make_read("H1-s400_2", "nodeB", 100, 1100, 500));

        svtigs["H1-s400_1"] = make_svtig("H1-s400_1");
        svtigs["H1-s400_2"] = make_svtig("H1-s400_2");

        auto result = remove_duplicates(reads, svtigs, extra);
        if (result.first != 0 || result.second != 2) {
            std::cerr << "Test 10 FAILED: different nodes should not be duplicates" << std::endl;
            return 1;
        }
        std::cout << "Test 10 passed: different nodes" << std::endl;
        for (auto* r : reads) delete r;
        for (auto& p : svtigs) delete p.second;
    }

    // Test 11: Fragmented assembly (extra underscore suffix)
    {
        std::vector<Read*> reads;
        std::map<std::string, SVtig*> svtigs;
        int extra = 0;

        reads.push_back(make_read("H1-s500_1", "nodeA", 100, 1100, 600));
        reads.push_back(make_read("H1-s500_1_2", "nodeA", 5000, 6000, 300));

        svtigs["H1-s500_1"] = make_svtig("H1-s500_1");

        auto result = remove_duplicates(reads, svtigs, extra);
        if (result.first != 0 || result.second != 2 || extra != 1) {
            std::cerr << "Test 11 FAILED: expected 0 dup, 2 kept, 1 extra, got "
                      << result.first << " dup, " << result.second << " kept, " << extra << " extra" << std::endl;
            return 1;
        }
        if (svtigs.count("H1-s500_1_2") != 1 || !svtigs["H1-s500_1_2"]->output) {
            std::cerr << "Test 11 FAILED: fragmented svtig not added correctly" << std::endl;
            return 1;
        }
        std::cout << "Test 11 passed: fragmented assembly" << std::endl;
        for (auto* r : reads) delete r;
        for (auto& p : svtigs) delete p.second;
    }

    // Test 12: Pre-marked duplicates are skipped
    {
        std::vector<Read*> reads;
        std::map<std::string, SVtig*> svtigs;
        int extra = 0;

        Read* r1 = make_read("H1-s600_1", "nodeA", 100, 1100, 600);
        Read* r2 = make_read("H1-s600_2", "nodeA", 100, 1100, 500);
        r2->duplicate = true;

        reads.push_back(r1);
        reads.push_back(r2);

        svtigs["H1-s600_1"] = make_svtig("H1-s600_1");
        svtigs["H1-s600_2"] = make_svtig("H1-s600_2");

        auto result = remove_duplicates(reads, svtigs, extra);
        if (result.first != 0 || result.second != 1) {
            std::cerr << "Test 12 FAILED: expected 0 new dup, 1 kept" << std::endl;
            return 1;
        }
        std::cout << "Test 12 passed: pre-marked duplicates" << std::endl;
        for (auto* r : reads) delete r;
        for (auto& p : svtigs) delete p.second;
    }

    // Test 13: Partial overlap below threshold - both kept
    {
        std::vector<Read*> reads;
        std::map<std::string, SVtig*> svtigs;
        int extra = 0;

        reads.push_back(make_read("H1-s700_1", "nodeA", 0, 1000, 500));
        reads.push_back(make_read("H1-s700_2", "nodeA", 900, 1900, 400));

        svtigs["H1-s700_1"] = make_svtig("H1-s700_1");
        svtigs["H1-s700_2"] = make_svtig("H1-s700_2");

        auto result = remove_duplicates(reads, svtigs, extra);
        if (result.first != 0 || result.second != 2) {
            std::cerr << "Test 13 FAILED: partial overlap should keep both" << std::endl;
            return 1;
        }
        std::cout << "Test 13 passed: partial overlap" << std::endl;
        for (auto* r : reads) delete r;
        for (auto& p : svtigs) delete p.second;
    }

    std::cout << "All remap tests passed" << std::endl;
    return 0;
}
