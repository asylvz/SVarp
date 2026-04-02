#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <cstdio>
#include <zlib.h>
#include "reference.h"
#include "common.h"

int main()
{
    // Test 1: read_gfa with a small compressed GFA file
    {
        const char* gfa_path = "/tmp/test_svarp.gfa.gz";
        gzFile gz = gzopen(gfa_path, "wb");
        if (!gz) { std::cerr << "Test 1: Failed to create gz file" << std::endl; return 1; }

        // Write GFA lines: S (segment) and L (link)
        std::string s1 = "S\ts1\tACGTACGT\tLN:i:8\tSN:Z:chr1\tSO:i:0\tSR:i:0\n";
        std::string s2 = "S\ts2\tTGCATGCA\tLN:i:8\tSN:Z:chr1\tSO:i:8\tSR:i:0\n";
        std::string s3 = "S\ts3\tAAAATTTT\tLN:i:8\tSN:Z:chr2\tSO:i:0\tSR:i:0\n";
        std::string l1 = "L\ts1\t+\ts2\t+\t0M\n";
        std::string l2 = "L\ts2\t+\ts3\t+\t0M\n";

        gzwrite(gz, s1.c_str(), s1.size());
        gzwrite(gz, s2.c_str(), s2.size());
        gzwrite(gz, l1.c_str(), l1.size());
        gzwrite(gz, s3.c_str(), s3.size());
        gzwrite(gz, l2.c_str(), l2.size());
        gzclose(gz);

        parameters params;
        params.ref_graph = gfa_path;

        std::map<std::string, Contig*> ref;
        std::map<std::string, gfaNode*> gfa;
        std::map<std::string, std::vector<std::string>> incoming, outgoing;

        int rc = read_gfa(params, ref, gfa, incoming, outgoing);
        if (rc != RETURN_SUCCESS) { std::cerr << "Test 1: read_gfa failed" << std::endl; return 1; }

        // Should have 3 nodes
        if (gfa.size() != 3) { std::cerr << "Test 1: Expected 3 nodes, got " << gfa.size() << std::endl; return 1; }

        // Verify node properties
        if (gfa["s1"]->len != 8) { std::cerr << "Test 1: s1 len should be 8" << std::endl; return 1; }
        if (gfa["s1"]->contig != "chr1") { std::cerr << "Test 1: s1 contig should be chr1" << std::endl; return 1; }
        if (gfa["s1"]->offset != 0) { std::cerr << "Test 1: s1 offset should be 0" << std::endl; return 1; }
        if (gfa["s2"]->offset != 8) { std::cerr << "Test 1: s2 offset should be 8" << std::endl; return 1; }
        if (gfa["s3"]->contig != "chr2") { std::cerr << "Test 1: s3 contig should be chr2" << std::endl; return 1; }

        // Verify contigs: chr1 (8+8=16), chr2 (8), overall (24)
        if (ref.find("chr1") == ref.end()) { std::cerr << "Test 1: chr1 contig missing" << std::endl; return 1; }
        if (ref["chr1"]->contig_length != 16) { std::cerr << "Test 1: chr1 len should be 16, got " << ref["chr1"]->contig_length << std::endl; return 1; }
        if (ref["chr2"]->contig_length != 8) { std::cerr << "Test 1: chr2 len should be 8" << std::endl; return 1; }
        if (ref["overall"]->contig_length != 24) { std::cerr << "Test 1: overall len should be 24, got " << ref["overall"]->contig_length << std::endl; return 1; }

        // Verify links
        if (outgoing.find("s1") == outgoing.end() || outgoing["s1"][0] != "s2") {
            std::cerr << "Test 1: s1->s2 link missing" << std::endl; return 1;
        }
        if (incoming.find("s2") == incoming.end() || incoming["s2"][0] != "s1") {
            std::cerr << "Test 1: s2<-s1 incoming missing" << std::endl; return 1;
        }

        std::remove(gfa_path);
        for (auto &kv : gfa) delete kv.second;
        for (auto &kv : ref) delete kv.second;
    }

    // Test 2: read_gfa with plain (uncompressed) GFA via gzopen
    {
        // gzopen can transparently read uncompressed files
        const char* gfa_path = "/tmp/test_svarp_plain.gfa";
        FILE* f = fopen(gfa_path, "w");
        fprintf(f, "S\tn1\tACGT\tLN:i:4\tSN:Z:chrX\tSO:i:0\tSR:i:0\n");
        fprintf(f, "S\tn2\tTGCA\tLN:i:4\tSN:Z:chrX\tSO:i:4\tSR:i:0\n");
        fclose(f);

        parameters params;
        params.ref_graph = gfa_path;

        std::map<std::string, Contig*> ref;
        std::map<std::string, gfaNode*> gfa;
        std::map<std::string, std::vector<std::string>> incoming, outgoing;

        int rc = read_gfa(params, ref, gfa, incoming, outgoing);
        if (rc != RETURN_SUCCESS) { std::cerr << "Test 2: read_gfa failed on plain file" << std::endl; return 1; }
        if (gfa.size() != 2) { std::cerr << "Test 2: Expected 2 nodes, got " << gfa.size() << std::endl; return 1; }
        if (ref["chrX"]->contig_length != 8) { std::cerr << "Test 2: chrX len should be 8" << std::endl; return 1; }

        std::remove(gfa_path);
        for (auto &kv : gfa) delete kv.second;
        for (auto &kv : ref) delete kv.second;
    }

    // Test 3: GFA with long lines (exceeding buffer)
    {
        const char* gfa_path = "/tmp/test_svarp_longline.gfa.gz";
        gzFile gz = gzopen(gfa_path, "wb");

        // Create a sequence longer than CHUNK (1MB)
        std::string long_seq(2000000, 'A'); // 2MB sequence
        std::string line = "S\tlong_node\t" + long_seq + "\tLN:i:2000000\tSN:Z:chr1\tSO:i:0\tSR:i:0\n";
        gzwrite(gz, line.c_str(), line.size());
        gzclose(gz);

        parameters params;
        params.ref_graph = gfa_path;

        std::map<std::string, Contig*> ref;
        std::map<std::string, gfaNode*> gfa;
        std::map<std::string, std::vector<std::string>> incoming, outgoing;

        int rc = read_gfa(params, ref, gfa, incoming, outgoing);
        if (rc != RETURN_SUCCESS) { std::cerr << "Test 3: read_gfa failed on long line" << std::endl; return 1; }
        if (gfa.size() != 1) { std::cerr << "Test 3: Expected 1 node, got " << gfa.size() << std::endl; return 1; }
        if (gfa["long_node"]->len != 2000000) { std::cerr << "Test 3: len should be 2000000" << std::endl; return 1; }
        if (gfa["long_node"]->sequence.size() != 2000000) { std::cerr << "Test 3: sequence size mismatch" << std::endl; return 1; }

        std::remove(gfa_path);
        for (auto &kv : gfa) delete kv.second;
        for (auto &kv : ref) delete kv.second;
    }

    // Test 4: Empty GFA file
    {
        const char* gfa_path = "/tmp/test_svarp_empty.gfa.gz";
        gzFile gz = gzopen(gfa_path, "wb");
        gzclose(gz);

        parameters params;
        params.ref_graph = gfa_path;

        std::map<std::string, Contig*> ref;
        std::map<std::string, gfaNode*> gfa;
        std::map<std::string, std::vector<std::string>> incoming, outgoing;

        int rc = read_gfa(params, ref, gfa, incoming, outgoing);
        if (rc != RETURN_SUCCESS) { std::cerr << "Test 4: read_gfa failed on empty" << std::endl; return 1; }
        if (gfa.size() != 0) { std::cerr << "Test 4: Expected 0 nodes" << std::endl; return 1; }

        std::remove(gfa_path);
        for (auto &kv : ref) delete kv.second;
    }

    // Test 5: GFA with malformed S lines (< 6 fields) should be skipped
    {
        const char* gfa_path = "/tmp/test_svarp_malformed.gfa.gz";
        gzFile gz = gzopen(gfa_path, "wb");
        std::string bad = "S\ts_bad\tACGT\n"; // only 3 fields
        std::string good = "S\ts_good\tACGT\tLN:i:4\tSN:Z:chr1\tSO:i:0\tSR:i:0\n";
        gzwrite(gz, bad.c_str(), bad.size());
        gzwrite(gz, good.c_str(), good.size());
        gzclose(gz);

        parameters params;
        params.ref_graph = gfa_path;

        std::map<std::string, Contig*> ref;
        std::map<std::string, gfaNode*> gfa;
        std::map<std::string, std::vector<std::string>> incoming, outgoing;

        int rc = read_gfa(params, ref, gfa, incoming, outgoing);
        if (rc != RETURN_SUCCESS) { std::cerr << "Test 5: read_gfa failed" << std::endl; return 1; }
        if (gfa.size() != 1) { std::cerr << "Test 5: Expected 1 node (bad skipped), got " << gfa.size() << std::endl; return 1; }
        if (gfa.find("s_good") == gfa.end()) { std::cerr << "Test 5: s_good should exist" << std::endl; return 1; }

        std::remove(gfa_path);
        for (auto &kv : gfa) delete kv.second;
        for (auto &kv : ref) delete kv.second;
    }

    // Test 6: contig_coverage with single-node path
    {
        std::map<std::string, gfaNode*> gfa;
        gfaNode* n1 = new gfaNode("s1", "", 1000, "chr1", 0);
        gfa["s1"] = n1;

        std::map<std::string, Contig*> ref;
        ref["chr1"] = new Contig();
        ref["overall"] = new Contig();

        Gaf g;
        g.path = ">s1";
        g.path_start = 100;
        g.path_end = 500;

        int mapped = contig_coverage(ref, gfa, g);
        if (mapped != 400) { std::cerr << "Test 6: Expected 400 mapped bases, got " << mapped << std::endl; return 1; }
        if (ref["chr1"]->mapped_bases != 400) { std::cerr << "Test 6: chr1 mapped_bases mismatch" << std::endl; return 1; }
        if (ref["chr1"]->mapped_reads != 1) { std::cerr << "Test 6: chr1 mapped_reads should be 1" << std::endl; return 1; }

        for (auto &kv : gfa) delete kv.second;
        for (auto &kv : ref) delete kv.second;
    }

    // Test 7: contig_coverage with multi-node path
    {
        std::map<std::string, gfaNode*> gfa;
        gfaNode* n1 = new gfaNode("s1", "", 100, "chr1", 0);
        gfaNode* n2 = new gfaNode("s2", "", 200, "chr1", 100);
        gfaNode* n3 = new gfaNode("s3", "", 150, "chr1", 300);
        gfa["s1"] = n1;
        gfa["s2"] = n2;
        gfa["s3"] = n3;

        std::map<std::string, Contig*> ref;
        ref["chr1"] = new Contig();
        ref["overall"] = new Contig();

        Gaf g;
        g.path = ">s1>s2>s3";
        g.path_start = 50;
        g.path_end = 400;
        // total_path_length = 400-50 = 350
        // first node: 100-50=50 bases
        // middle node: 200 bases
        // last node: 350-50-200=100 bases

        int mapped = contig_coverage(ref, gfa, g);
        if (mapped != 350) { std::cerr << "Test 7: Expected 350, got " << mapped << std::endl; return 1; }
        if (ref["chr1"]->mapped_reads != 3) { std::cerr << "Test 7: mapped_reads should be 3 (one per node)" << std::endl; return 1; }

        for (auto &kv : gfa) delete kv.second;
        for (auto &kv : ref) delete kv.second;
    }

    std::cout << "reference test passed" << std::endl;
    return 0;
}
