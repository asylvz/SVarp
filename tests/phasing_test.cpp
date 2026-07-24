#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <cstdio>
#include "phasing.h"
#include "variant.h"
#include "common.h"

int main() {
    // Test 1: read_phase_file with a TSV file
    {
        const char* tsv_path = "/tmp/test_svarp_phase.tsv";
        std::ofstream f(tsv_path);
        f << "# header comment\n";
        f << "read_001\tH1\tps_100\tchr1\n";
        f << "read_002\tH2\tps_100\tchr1\n";
        f << "read_003\tH1\tps_200\tchr1\n";
        f << "read_004\tnone\tnone\tchr1\n";
        f.close();

        parameters params;
        params.phase_tags = tsv_path;
        std::map<std::string, phase*> phased_reads;

        int rc = read_phase_file(params, phased_reads);
        if (rc != RETURN_SUCCESS) { std::cerr << "Test 1: read_phase_file failed" << std::endl; return 1; }
        if (phased_reads.size() != 4) { std::cerr << "Test 1: Expected 4 reads, got " << phased_reads.size() << std::endl; return 1; }
        if (phased_reads["read_001"]->haplotype != "H1") { std::cerr << "Test 1: read_001 should be H1" << std::endl; return 1; }
        if (phased_reads["read_002"]->haplotype != "H2") { std::cerr << "Test 1: read_002 should be H2" << std::endl; return 1; }
        if (phased_reads["read_004"]->haplotype != "none") { std::cerr << "Test 1: read_004 should be none" << std::endl; return 1; }

        // Verify comment line was skipped
        if (phased_reads.find("#") != phased_reads.end()) { std::cerr << "Test 1: comment should be skipped" << std::endl; return 1; }

        for (auto &kv : phased_reads) delete kv.second;
        std::remove(tsv_path);
    }

    // Test 2: phase_svs assigns reads to haplotypes
    {
        // Create phase data
        std::map<std::string, phase*> phased_reads;

        phase* p1 = new phase(); p1->read_name = "read_A"; p1->haplotype = "H1"; p1->phase_set = "ps1"; p1->contig = "chr1";
        phase* p2 = new phase(); p2->read_name = "read_B"; p2->haplotype = "H2"; p2->phase_set = "ps1"; p2->contig = "chr1";
        phase* p3 = new phase(); p3->read_name = "read_C"; p3->haplotype = "H1"; p3->phase_set = "ps1"; p3->contig = "chr1";
        phase* p4 = new phase(); p4->read_name = "read_D"; p4->haplotype = "none"; p4->phase_set = "none"; p4->contig = "chr1";

        phased_reads["read_A"] = p1;
        phased_reads["read_B"] = p2;
        phased_reads["read_C"] = p3;
        phased_reads["read_D"] = p4;

        // Create SV clusters with untagged reads
        std::map<std::string, std::vector<SVCluster*>> vars;
        SVCluster* sv1 = new SVCluster();
        sv1->node = "node1";
        sv1->start_pos = 100;
        sv1->contig = "chr1";
        sv1->reads_untagged.insert("read_A");
        sv1->reads_untagged.insert("read_B");
        sv1->reads_untagged.insert("read_C");
        sv1->reads_untagged.insert("read_D");
        sv1->reads_untagged.insert("read_E"); // not in phase file

        vars["node1"].push_back(sv1);

        // Run phase_svs
        phase_svs(phased_reads, vars);

        // Check: read_A and read_C should move to reads_h1
        if (sv1->reads_h1.size() != 2) { std::cerr << "Test 2: Expected 2 H1 reads, got " << sv1->reads_h1.size() << std::endl; return 1; }
        if (sv1->reads_h1.find("read_A") == sv1->reads_h1.end()) { std::cerr << "Test 2: read_A should be in H1" << std::endl; return 1; }
        if (sv1->reads_h1.find("read_C") == sv1->reads_h1.end()) { std::cerr << "Test 2: read_C should be in H1" << std::endl; return 1; }

        // Check: read_B should move to reads_h2
        if (sv1->reads_h2.size() != 1) { std::cerr << "Test 2: Expected 1 H2 read, got " << sv1->reads_h2.size() << std::endl; return 1; }

        // Check: read_D (none) and read_E (not in file) should stay untagged
        if (sv1->reads_untagged.size() != 2) { std::cerr << "Test 2: Expected 2 untagged, got " << sv1->reads_untagged.size() << std::endl; return 1; }
        if (sv1->reads_untagged.find("read_E") == sv1->reads_untagged.end()) { std::cerr << "Test 2: read_E should stay untagged" << std::endl; return 1; }

        for (auto &kv : phased_reads) delete kv.second;
        delete sv1;
    }

    // Test 3: Empty phase file
    {
        parameters params;
        params.phase_tags = "";
        std::map<std::string, phase*> phased_reads;

        int rc = read_phase_file(params, phased_reads);
        if (rc != RETURN_ERROR) { std::cerr << "Test 3: Expected error for empty path" << std::endl; return 1; }
    }

    // Test 4: a non-existent (but non-empty) phase path must not silently succeed
    {
        parameters params;
        params.phase_tags = "/tmp/svarp_no_such_phase_file_12345.tsv";
        std::map<std::string, phase*> phased_reads;

        int rc = read_phase_file(params, phased_reads);
        if (rc != RETURN_ERROR) { std::cerr << "Test 4: Expected error for missing phase file, got SUCCESS (silent unphased run)" << std::endl; return 1; }
    }

    // Test 5: malformed lines (empty, or fewer than 4 fields) are skipped
    {
        const char* tsv_path = "/tmp/test_svarp_phase_malformed.tsv";
        std::ofstream f(tsv_path);
        f << "read_001\tH1\tps_100\tchr1\n";  // valid
        f << "\n";                            // empty line
        f << "read_002\tH2\n";                // only 2 fields
        f.close();

        parameters params;
        params.phase_tags = tsv_path;
        std::map<std::string, phase*> phased_reads;

        int rc = read_phase_file(params, phased_reads);
        if (rc != RETURN_SUCCESS) { std::cerr << "Test 5: read_phase_file failed" << std::endl; return 1; }
        if (phased_reads.size() != 1) { std::cerr << "Test 5: Expected only the valid line, got " << phased_reads.size() << std::endl; return 1; }

        for (auto &kv : phased_reads) delete kv.second;
        std::remove(tsv_path);
    }

    std::cout << "phasing test passed" << std::endl;
    return 0;
}
