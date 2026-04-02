#include <iostream>
#include <string>
#include <map>
#include <set>
#include <fstream>
#include <cstdio>
#include <zlib.h>
#include "alignment.h"
#include "variant.h"
#include "reference.h"
#include "common.h"

// Reset global counters defined in alignment.cpp
extern int primary_cnt, secondary_cnt, inter_cnt, intra_cnt, insertion_cnt, deletion_cnt, mismatch_cnt;

static void reset_counters() {
    primary_cnt = 0; secondary_cnt = 0; inter_cnt = 0; intra_cnt = 0;
    insertion_cnt = 0; deletion_cnt = 0; mismatch_cnt = 0;
}

int main()
{
    // Test 1: read_gz with a small compressed GAF file
    {
        reset_counters();

        // Create a small GFA-like graph
        std::map<std::string, gfaNode*> gfa;
        gfaNode* n1 = new gfaNode("s1", "ACGT", 1000, "chr1", 0);
        gfaNode* n2 = new gfaNode("s2", "TGCA", 500, "chr1", 1000);
        gfa["s1"] = n1;
        gfa["s2"] = n2;

        std::map<std::string, Contig*> ref;
        Contig* c1 = new Contig();
        c1->contig_length = 1500;
        ref["chr1"] = c1;
        Contig* overall = new Contig();
        overall->contig_length = 1500;
        ref["overall"] = overall;

        // Write a small compressed GAF file
        const char* gz_path = "/tmp/test_svarp_read_gz.gaf.gz";
        gzFile gz = gzopen(gz_path, "wb");
        if (!gz) { std::cerr << "Test 1: Failed to create gz file" << std::endl; return 1; }

        // Two valid primary reads mapping to s1, one with SV in cigar
        // read1: simple match
        std::string line1 = "read1\t1000\t0\t1000\t+\t>s1\t1000\t0\t1000\t990\t1000\t60\ttp:A:P\tcg:Z:1000M\n";
        // read2: has a 100bp insertion (>= MINSVSIZE=50)
        std::string line2 = "read2\t1000\t0\t1000\t+\t>s1\t1000\t0\t900\t890\t900\t60\ttp:A:P\tcg:Z:400M100I500M\n";
        // read3: same as read1 (makes read1 multi-mapped -> inter-alignment SV)
        std::string line3 = "read1\t1000\t0\t1000\t+\t>s2\t500\t0\t500\t490\t500\t60\ttp:A:P\tcg:Z:500M\n";
        // read4: low mapq, should be filtered
        std::string line4 = "read4\t500\t0\t500\t+\t>s1\t1000\t0\t500\t490\t500\t2\ttp:A:P\tcg:Z:500M\n";
        // read5: secondary alignment, should be filtered
        std::string line5 = "read5\t800\t0\t800\t+\t>s1\t1000\t0\t800\t790\t800\t60\ttp:A:S\tcg:Z:800M\n";

        gzwrite(gz, line1.c_str(), line1.size());
        gzwrite(gz, line2.c_str(), line2.size());
        gzwrite(gz, line3.c_str(), line3.size());
        gzwrite(gz, line4.c_str(), line4.size());
        gzwrite(gz, line5.c_str(), line5.size());
        gzclose(gz);

        parameters params;
        params.gaf = gz_path;
        params.fp_logs.open("/tmp/test_read_gz.log");

        std::map<std::string, Variant*> vars;
        std::set<std::string> unmapped;
        std::map<std::string, int> read_freq;

        int rc = read_gz(params, ref, gfa, vars, unmapped, read_freq);
        if (rc != RETURN_SUCCESS) { std::cerr << "Test 1: read_gz failed" << std::endl; return 1; }

        // read1 appears twice -> should be in read_freq
        if (read_freq.find("read1") == read_freq.end()) {
            std::cerr << "Test 1: read1 should be in read_freq (multi-mapped)" << std::endl; return 1;
        }
        if (read_freq["read1"] != 2) {
            std::cerr << "Test 1: read1 freq should be 2, got " << read_freq["read1"] << std::endl; return 1;
        }

        // read2 has 100bp insertion >= MINSVSIZE, should create a variant
        if (vars.empty()) {
            std::cerr << "Test 1: Expected at least one variant from read2's 100bp insertion" << std::endl; return 1;
        }

        // primary_cnt should be 3 (read1, read2, read1 second mapping)
        // read4 (low mapq) and read5 (secondary) should be filtered
        if (primary_cnt != 3) {
            std::cerr << "Test 1: Expected 3 primary, got " << primary_cnt << std::endl; return 1;
        }

        params.fp_logs.close();
        std::remove(gz_path);
        std::remove("/tmp/test_read_gz.log");

        for (auto &kv : vars) delete kv.second;
        for (auto &kv : ref) delete kv.second;
        for (auto &kv : gfa) delete kv.second;
    }

    // Test 2: read_gz with empty file
    {
        reset_counters();

        const char* gz_path = "/tmp/test_svarp_empty.gaf.gz";
        gzFile gz = gzopen(gz_path, "wb");
        gzclose(gz);

        std::map<std::string, gfaNode*> gfa;
        std::map<std::string, Contig*> ref;
        Contig* overall = new Contig();
        ref["overall"] = overall;

        parameters params;
        params.gaf = gz_path;
        params.fp_logs.open("/tmp/test_empty_gz.log");

        std::map<std::string, Variant*> vars;
        std::set<std::string> unmapped;
        std::map<std::string, int> read_freq;

        int rc = read_gz(params, ref, gfa, vars, unmapped, read_freq);
        if (rc != RETURN_SUCCESS) { std::cerr << "Test 2: read_gz failed on empty file" << std::endl; return 1; }
        if (!vars.empty()) { std::cerr << "Test 2: Expected no variants" << std::endl; return 1; }
        if (primary_cnt != 0) { std::cerr << "Test 2: Expected 0 primary" << std::endl; return 1; }

        params.fp_logs.close();
        std::remove(gz_path);
        std::remove("/tmp/test_empty_gz.log");
        for (auto &kv : ref) delete kv.second;
    }

    std::cout << "alignment.read_gz test passed" << std::endl;
    return 0;
}
