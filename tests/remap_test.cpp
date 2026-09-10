#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <fstream>
#include <cstdio>
#include <htslib/faidx.h>
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

    // Test 14: A kept svtig carries the path and map ratio of its remapping
    {
        std::vector<Read*> reads;
        std::map<std::string, SVtig*> svtigs;
        int extra = 0;

        Read* r = make_read("H1-s800_1", ">s800<s801", 100, 1100, 500);
        r->cov = 0.42;
        reads.push_back(r);

        svtigs["H1-s800_1"] = make_svtig("H1-s800_1");

        remove_duplicates(reads, svtigs, extra);
        if (svtigs["H1-s800_1"]->remap_path != ">s800<s801") {
            std::cerr << "Test 14 FAILED: expected remap_path=>s800<s801, got "
                      << svtigs["H1-s800_1"]->remap_path << std::endl;
            return 1;
        }
        if (svtigs["H1-s800_1"]->map_ratio != 0.42) {
            std::cerr << "Test 14 FAILED: expected map_ratio=0.42, got "
                      << svtigs["H1-s800_1"]->map_ratio << std::endl;
            return 1;
        }
        std::cout << "Test 14 passed: remapping recorded on svtig" << std::endl;
        for (auto* p : reads) delete p;
        for (auto& p : svtigs) delete p.second;
    }

    // Test 15: Svtigs from a fragmented assembly carry it as well
    {
        std::vector<Read*> reads;
        std::map<std::string, SVtig*> svtigs;
        int extra = 0;

        reads.push_back(make_read("H1-s900_1", ">s900", 100, 1100, 600));
        Read* frag = make_read("H1-s900_1_2", ">s901", 5000, 6000, 300);
        frag->cov = 0.33;
        reads.push_back(frag);

        svtigs["H1-s900_1"] = make_svtig("H1-s900_1");

        remove_duplicates(reads, svtigs, extra);
        if (svtigs.count("H1-s900_1_2") != 1) {
            std::cerr << "Test 15 FAILED: fragmented svtig missing" << std::endl;
            return 1;
        }
        if (svtigs["H1-s900_1_2"]->remap_path != ">s901" || svtigs["H1-s900_1_2"]->map_ratio != 0.33) {
            std::cerr << "Test 15 FAILED: fragmented svtig lost its remapping" << std::endl;
            return 1;
        }
        std::cout << "Test 15 passed: fragmented svtig keeps remapping" << std::endl;
        for (auto* p : reads) delete p;
        for (auto& p : svtigs) delete p.second;
    }

    // Test 16: Header reports the remapping only when the svtig was remapped
    {
        SVtig* s = make_svtig("H1-s1000_400");
        s->pos = 1234567;
        s->reads.insert("read_a");
        s->reads.insert("read_b");

        std::string bare = svtig_header(s);
        if (bare != ">H1-s1000_400 contig=chr1 pos=1234567 support=2") {
            std::cerr << "Test 16 FAILED: unexpected header without remapping: " << bare << std::endl;
            return 1;
        }

        s->remap_path = ">s1000<s1001";
        s->map_ratio = 0.75;
        std::string mapped = svtig_header(s);
        if (mapped != ">H1-s1000_400 contig=chr1 pos=1234567 support=2 path=>s1000<s1001 graph_cov=0.75 max_gap=0 max_indel=0 graph_explained=no") {
            std::cerr << "Test 16 FAILED: unexpected header with remapping: " << mapped << std::endl;
            return 1;
        }
        std::cout << "Test 16 passed: header reports remapping" << std::endl;
        delete s;
    }

    // Test 17: The two haplotypes of one locus are alleles, not duplicates
    {
        std::vector<Read*> reads;
        std::map<std::string, SVtig*> svtigs;
        int extra = 0;

        // Same graph path, overlapping footprints: a homozygous SV assembled twice.
        reads.push_back(make_read("H1-s1234_500", ">s1234>s1235", 400, 3400, 3000));
        reads.push_back(make_read("H2-s1234_500", ">s1234>s1235", 418, 3392, 2974));

        svtigs["H1-s1234_500"] = make_svtig("H1-s1234_500");
        svtigs["H2-s1234_500"] = make_svtig("H2-s1234_500");

        auto result = remove_duplicates(reads, svtigs, extra);
        if (result.first != 0 || result.second != 2) {
            std::cerr << "Test 17 FAILED: expected 0 dup, 2 kept, got "
                      << result.first << " dup, " << result.second << " kept" << std::endl;
            return 1;
        }
        if (!svtigs["H1-s1234_500"]->output || !svtigs["H2-s1234_500"]->output) {
            std::cerr << "Test 17 FAILED: a haplotype lost its svtig" << std::endl;
            return 1;
        }
        std::cout << "Test 17 passed: haplotypes are not duplicates" << std::endl;
        for (auto* r : reads) delete r;
        for (auto& p : svtigs) delete p.second;
    }

    // Test 18: The untagged copy survives alongside both haplotypes
    {
        std::vector<Read*> reads;
        std::map<std::string, SVtig*> svtigs;
        int extra = 0;

        reads.push_back(make_read("H1-s5000_7", ">s5000", 400, 3400, 3000));
        reads.push_back(make_read("H2-s5000_7", ">s5000", 410, 3390, 2900));
        reads.push_back(make_read("None-s5000_7", ">s5000", 405, 3395, 2950));

        svtigs["H1-s5000_7"] = make_svtig("H1-s5000_7");
        svtigs["H2-s5000_7"] = make_svtig("H2-s5000_7");
        svtigs["None-s5000_7"] = make_svtig("None-s5000_7");

        auto result = remove_duplicates(reads, svtigs, extra);
        if (result.first != 0 || result.second != 3) {
            std::cerr << "Test 18 FAILED: expected 0 dup, 3 kept, got "
                      << result.first << " dup, " << result.second << " kept" << std::endl;
            return 1;
        }
        std::cout << "Test 18 passed: untagged svtig kept with both haplotypes" << std::endl;
        for (auto* r : reads) delete r;
        for (auto& p : svtigs) delete p.second;
    }

    // Test 19: MINSVSIZE is the smallest callable SV, so a realignment that
    // diverges by exactly that much still counts as carrying one
    {
        if (!cigar_has_sv("100M50D100M")) {
            std::cerr << "Test 19 FAILED: 50 bp deletion not counted as an SV" << std::endl;
            return 1;
        }
        if (!cigar_has_sv("100M50I100M")) {
            std::cerr << "Test 19 FAILED: 50 bp insertion not counted as an SV" << std::endl;
            return 1;
        }
        if (cigar_has_sv("100M49D100M")) {
            std::cerr << "Test 19 FAILED: 49 bp deletion counted as an SV" << std::endl;
            return 1;
        }
        if (!cigar_has_sv("100M120I100M")) {
            std::cerr << "Test 19 FAILED: 120 bp insertion not counted as an SV" << std::endl;
            return 1;
        }
        if (cigar_has_sv("300M")) {
            std::cerr << "Test 19 FAILED: a gapless alignment counted as an SV" << std::endl;
            return 1;
        }
        std::cout << "Test 19 passed: SV bound matches MINSVSIZE" << std::endl;
    }

    // Test 20: alternative alleles read off a remap path
    {
        std::map<std::string, gfaNode*> gfa;
        auto add = [&gfa](const std::string& name, int rank, const std::string& contig) {
            gfaNode* n = new gfaNode();
            n->name = name;
            n->rank = rank;
            n->contig = contig;
            n->len = 100;
            gfa[name] = n;
        };
        add("r1", 0, "CHM13#0#chr21");
        add("r2", 0, "CHM13#0#chr21");
        add("r3", 0, "CHM13#0#chr21");
        add("a1", 1, "HG01252#2#JBHIHZ010000014.1");
        add("a2", 1, "HG01252#2#JBHIHZ010000014.1");
        add("a3", 1, "NA21102#1#JBIREH010000004.1");
        add("u1", -1, "CHM13#0#chr21");   // GFA carried no SR tag

        struct { const char* path; const char* want; const char* why; } cases[] = {
            {">r1>r2>r3", "", "a path that stays on the reference reports nothing"},
            {">r1>a1>r2", "a1:HG01252#2#JBHIHZ010000014.1", "a lone allele needs no range"},
            {">r1>a1>a2>r2", "a1-a2:HG01252#2#JBHIHZ010000014.1",
             "consecutive nodes from one donor are one allele"},
            {">r1>a1>a2>a3>r2",
             "a1-a2:HG01252#2#JBHIHZ010000014.1;a3:NA21102#1#JBIREH010000004.1",
             "a change of donor splits the block"},
            {">r1>a1>r2>a3>r3",
             "a1:HG01252#2#JBHIHZ010000014.1;a3:NA21102#1#JBIREH010000004.1",
             "a reference node between two alleles splits them"},
            {"<r3<a1<r1", "a1:HG01252#2#JBHIHZ010000014.1",
             "orientation does not change which nodes are alleles"},
            {">r1>a1>a1>r2", "a1:HG01252#2#JBHIHZ010000014.1",
             "a node repeated back to back is still one allele"},
            {">r1>missing>r2", "", "a node absent from the graph withholds the whole list"},
            {"chr21:100-900", "", "a stable-sequence path carries no nodes"},
            {"", "", "an empty path reports nothing"},
            {">r1>u1>r2", "", "without an SR tag a node cannot be called an allele"},
        };

        for (const auto& c : cases) {
            std::string got = collect_alt_nodes(c.path, gfa);
            if (got != c.want) {
                std::cerr << "Test 20 FAILED (" << c.why << "): path=" << c.path
                          << " expected \"" << c.want << "\" got \"" << got << "\"" << std::endl;
                return 1;
            }
        }
        std::cout << "Test 20 passed: alt_nodes read off the remap path" << std::endl;

        // Test 21: filling alt_nodes and printing it in the header
        {
            std::map<std::string, SVtig*> svtigs;
            auto mk = [&svtigs](const std::string& name, bool output,
                                const std::string& path, double ratio) {
                SVtig* s = new SVtig();
                s->name = name;
                s->contig = "CHM13#0#chr21";
                s->pos = 100;
                s->output = output;
                s->remap_path = path;
                s->map_ratio = ratio;
                svtigs[name] = s;
                return s;
            };
            SVtig* kept = mk("kept", true, ">r1>a1>r2", 0.9);
            SVtig* dropped = mk("dropped", false, ">r1>a1>r2", 0.9);
            SVtig* noremap = mk("noremap", true, "", -1);

            fill_alt_nodes(svtigs, gfa);

            if (kept->alt_nodes != "a1:HG01252#2#JBHIHZ010000014.1") {
                std::cerr << "Test 21 FAILED: kept svtig got \"" << kept->alt_nodes << "\"" << std::endl;
                return 1;
            }
            if (!dropped->alt_nodes.empty()) {
                std::cerr << "Test 21 FAILED: filtered svtig was filled" << std::endl;
                return 1;
            }
            if (!noremap->alt_nodes.empty()) {
                std::cerr << "Test 21 FAILED: svtig without a path was filled" << std::endl;
                return 1;
            }

            std::string h = svtig_header(kept);
            if (h.find(" alt_nodes=a1:HG01252#2#JBHIHZ010000014.1") == std::string::npos) {
                std::cerr << "Test 21 FAILED: header missing alt_nodes: " << h << std::endl;
                return 1;
            }
            if (svtig_header(noremap).find("alt_nodes=") != std::string::npos) {
                std::cerr << "Test 21 FAILED: unremapped svtig printed alt_nodes" << std::endl;
                return 1;
            }
            kept->alt_nodes.clear();
            if (svtig_header(kept).find("alt_nodes=") != std::string::npos) {
                std::cerr << "Test 21 FAILED: empty alt_nodes still printed" << std::endl;
                return 1;
            }
            for (auto& s : svtigs)
                delete s.second;
            std::cout << "Test 21 passed: alt_nodes filled and printed" << std::endl;
        }

        for (auto& n : gfa)
            delete n.second;
    }

    // Tests 22-25: update_read picks the representative record
    {
        auto gaf = [](int qs, int qe, const std::string& path, int ps, int pe) {
            Gaf g; g.query_start = qs; g.query_end = qe; g.query_length = 10000;
            g.path = path; g.path_start = ps; g.path_end = pe; return g;
        };
        // 22: a later primary record with an SV sets sv_in_cigar; freq counts primaries
        Read r; r.freq = 0; r.highest_map_ratio = 0;
        update_read(&r, gaf(0, 4000, ">a", 0, 4000), true, false, 0.4);
        update_read(&r, gaf(4000, 9000, ">b", 0, 5000), true, true, 0.5);
        if (!r.sv_in_cigar || r.freq != 2 || r.node != ">b" || r.highest_map_ratio != 0.5) {
            std::cerr << "Test 22 FAILED: later primary record not folded in" << std::endl; return 1;
        }
        std::cout << "Test 22 passed: later primary record updates sv_in_cigar" << std::endl;
        // 23: a low-MAPQ record never overrides a primary one
        update_read(&r, gaf(0, 9900, ">c", 0, 9900), false, true, 0.99);
        if (r.node != ">b" || r.highest_map_ratio != 0.5 || r.freq != 2) {
            std::cerr << "Test 23 FAILED: low-MAPQ record overrode the primary" << std::endl; return 1;
        }
        std::cout << "Test 23 passed: low-MAPQ record ignored" << std::endl;
        // 24: low-MAPQ records stand in while nothing better exists, then yield
        Read q; q.freq = 0; q.highest_map_ratio = 0;
        update_read(&q, gaf(0, 9900, ">c", 0, 9900), false, true, 0.99);
        if (q.node != ">c" || q.freq != 0 || q.sv_in_cigar) {
            std::cerr << "Test 24 FAILED: fallback record not recorded as expected" << std::endl; return 1;
        }
        update_read(&q, gaf(0, 3000, ">d", 0, 3000), true, false, 0.3);
        if (q.node != ">d" || q.highest_map_ratio != 0.3 || q.freq != 1) {
            std::cerr << "Test 24 FAILED: primary record did not replace the fallback" << std::endl; return 1;
        }
        std::cout << "Test 24 passed: fallback yields to primary" << std::endl;
        // 25: the representative is the longest query span, not the longest path span
        Read w; w.freq = 0; w.highest_map_ratio = 0;
        update_read(&w, gaf(0, 3000, ">e", 0, 8000), true, false, 0.3);
        update_read(&w, gaf(3000, 9000, ">f", 0, 6000), true, false, 0.6);
        if (w.node != ">f" || w.start != 0 || w.end != 6000) {
            std::cerr << "Test 25 FAILED: representative chosen by path span" << std::endl; return 1;
        }
        std::cout << "Test 25 passed: representative by query span" << std::endl;
    }

    // Tests 26-29: merged indels and the explained-by-graph decision
    {
        if (merged_indel("100=30I5=25I100=") != 55) { std::cerr << "Test 26 FAILED: split insertion not merged" << std::endl; return 1; }
        if (merged_indel("100=30I30=25I100=") != 30) { std::cerr << "Test 26 FAILED: distant insertions merged" << std::endl; return 1; }
        if (merged_indel("100=30I5=25D100=") != 30) { std::cerr << "Test 26 FAILED: insertion and deletion merged" << std::endl; return 1; }
        std::cout << "Test 26 passed: merged indel" << std::endl;

        Read full; full.svtig_size = 10000; full.ivals = {{0, 6000}, {5500, 10000}};
        graph_fit(&full);
        if (full.cov != 1.0 || full.max_gap != 0 || !explained_by_graph(&full)) { std::cerr << "Test 27 FAILED: fully covered svtig not explained" << std::endl; return 1; }
        std::cout << "Test 27 passed: covered svtig explained" << std::endl;

        Read gap; gap.svtig_size = 10000; gap.ivals = {{0, 4000}, {4600, 10000}};
        graph_fit(&gap);
        if (gap.max_gap != 600 || explained_by_graph(&gap)) { std::cerr << "Test 28 FAILED: 600 bp hole means the graph lacks the allele" << std::endl; return 1; }
        Read ind; ind.svtig_size = 10000; ind.ivals = {{0, 10000}}; ind.max_indel = 120;
        graph_fit(&ind);
        if (explained_by_graph(&ind)) { std::cerr << "Test 28 FAILED: 120 bp indel means the graph lacks the allele" << std::endl; return 1; }
        // coverage is the output gate, not part of the explained decision
        Read ends; ends.svtig_size = 10000; ends.ivals = {{2000, 8000}};
        graph_fit(&ends);
        if (ends.cov != 0.6 || !explained_by_graph(&ends)) { std::cerr << "Test 28 FAILED: 60% coverage is a gate matter, the aligned part is explained" << std::endl; return 1; }
        parameters gate; if (ends.cov >= gate.min_graph_cov) { std::cerr << "Test 28 FAILED: 60% coverage should fail the default gate" << std::endl; return 1; }
        Read none; none.svtig_size = 10000;
        graph_fit(&none);
        if (none.cov != 0 || none.cov >= gate.min_graph_cov) { std::cerr << "Test 28 FAILED: uncovered svtig passed the gate" << std::endl; return 1; }
        if (gate.min_svtig_len != 5000 || !gate.skip_untagged) { std::cerr << "Test 28 FAILED: default gate values" << std::endl; return 1; }
        std::cout << "Test 28 passed: holes and indels mark novel alleles, coverage gates the output" << std::endl;

        Read small; small.svtig_size = 10000; small.ivals = {{0, 5000}, {5040, 10000}}; small.max_indel = 49;
        graph_fit(&small);
        if (!explained_by_graph(&small)) { std::cerr << "Test 29 FAILED: sub-SV hole and indel should count as explained" << std::endl; return 1; }
        std::cout << "Test 29 passed: sub-SV differences do not count" << std::endl;
    }

    // Test 34: haplotype comes from the name prefix, not from H1/H2 inside the node name
    {
        if (haplotype_of("H1-sH2x_5") != "H1" || haplotype_of("None-sH1_3") != "None" || haplotype_of("H2-s1_1_2") != "H2" || haplotype_of("s1_1") != "") {
            std::cerr << "Test 34 FAILED: haplotype_of" << std::endl; return 1;
        }
        const char* fa = "/tmp/test_svarp_final_svtigs.fa";
        { std::ofstream f(fa); f << ">H1-sH2x_5\nACGTACGT\n>H2-s1_1\nAAAACCCC\n>None-sH1_3\nGGGGTTTT\n"; }
        faidx_t* fai = fai_load(fa);
        if (!fai) { std::cerr << "Test 34 FAILED: fai_load" << std::endl; return 1; }
        std::map<std::string, SVtig*> svtigs;
        for (const char* n : {"H1-sH2x_5", "H2-s1_1", "None-sH1_3"}) { svtigs[n] = make_svtig(n); svtigs[n]->output = true; }
        int ok = 1;
        for (const char* hap : {"H1", "H2", "None"}) {
            std::string out = std::string("/tmp/test_svarp_final_") + hap + ".fa";
            int n = write_final_svtigs(fai, svtigs, out, hap);
            std::ifstream in(out); std::string line; std::getline(in, line);
            if (n != 1 || line.rfind(std::string(">") + hap + "-", 0) != 0) { std::cerr << "Test 34 FAILED: " << hap << " file got " << n << " svtigs, first " << line << std::endl; ok = 0; }
            std::remove(out.c_str());
        }
        fai_destroy(fai);
        std::remove(fa); std::remove("/tmp/test_svarp_final_svtigs.fa.fai");
        for (auto& p : svtigs) delete p.second;
        if (!ok) return 1;
        std::cout << "Test 34 passed: haplotype files follow the name prefix" << std::endl;
    }

    // Test 35: reference_colinear - only an unbroken walk over reference nodes counts as the reference itself
    {
        std::map<std::string, gfaNode*> gfa;
        auto add = [&](const std::string& n, int len, const std::string& contig, int off, int rank) {
            gfa[n] = new gfaNode(n, std::string(len, 'A'), len, contig, off); gfa[n]->rank = rank; };
        add("r1", 100, "chr1", 0, 0); add("r2", 50, "chr1", 100, 0); add("r3", 200, "chr1", 150, 0);
        add("a1", 60, "HG002#1#chr1", 5000, 1); add("c1", 100, "chr2", 0, 0); add("u1", 100, "chr1", 350, -1);
        struct Case { const char* path; bool expect; const char* why; };
        Case cases[] = {
            {">r1>r2>r3", true, "forward walk"}, {"<r3<r2<r1", true, "reverse walk"}, {">r2", true, "single node"},
            {">r1>r3", false, "skipped node = deletion edge"}, {">r1>a1>r3", false, "alt node"},
            {">r1<r2>r3", false, "inversion"}, {">r3>r1", false, "out of order"}, {">r1>r2>r3>u1", false, "node without SR tag"},
            {">r1>c1", false, "contig change"}, {">r1>zz", false, "unknown node"}, {"chr1:1-500", false, "stable name path"},
        };
        for (auto& c : cases)
            if (reference_colinear(c.path, gfa) != c.expect) { std::cerr << "Test 35 FAILED: " << c.path << " (" << c.why << ")" << std::endl; return 1; }
        for (auto& kv : gfa) delete kv.second;
        std::cout << "Test 35 passed: reference-colinear paths" << std::endl;
    }

    std::cout << "All remap tests passed" << std::endl;
    return 0;
}
