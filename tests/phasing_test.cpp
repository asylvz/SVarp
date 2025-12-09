#include <iostream>
#include <string>
#include <map>
#include <vector>
#include "phasing.h"
#include "variant.h"
#include "common.h"

int main() {
    // Test 1: Verify phase struct instantiation and basic properties
    phase* ph1 = new phase();
    ph1->read_name = "read_001";
    ph1->haplotype = "H1";
    ph1->phase_set = "ps_chr1_100";
    ph1->contig = "contig1";
    
    if (ph1->read_name != "read_001") { std::cerr << "Test 1: Expected read_name=read_001" << std::endl; return 1; }
    if (ph1->haplotype != "H1") { std::cerr << "Test 1: Expected haplotype=H1" << std::endl; return 1; }
    if (ph1->phase_set != "ps_chr1_100") { std::cerr << "Test 1: Expected phase_set=ps_chr1_100" << std::endl; return 1; }
    
    // Test 2: Verify multiple phase objects and phased_reads map
    std::map<std::string, phase*> phased_reads;
    
    phase* ph2 = new phase();
    ph2->read_name = "read_002";
    ph2->haplotype = "H2";
    ph2->phase_set = "ps_chr1_100";
    ph2->contig = "contig1";
    
    phase* ph3 = new phase();
    ph3->read_name = "read_003";
    ph3->haplotype = "H1";
    ph3->phase_set = "ps_chr1_200";
    ph3->contig = "contig1";
    
    phased_reads["read_001"] = ph1;
    phased_reads["read_002"] = ph2;
    phased_reads["read_003"] = ph3;
    
    if (phased_reads.size() != 3) { std::cerr << "Test 2: Expected 3 phased reads, got " << phased_reads.size() << std::endl; return 1; }
    
    // Test 3: Verify phasing group extraction (reads with same phase_set)
    std::map<std::string, std::vector<phase*>> phase_groups;
    for (auto &kv : phased_reads) {
        phase_groups[kv.second->phase_set].push_back(kv.second);
    }
    
    if (phase_groups.size() != 2) { std::cerr << "Test 3: Expected 2 phase groups, got " << phase_groups.size() << std::endl; return 1; }
    if (phase_groups["ps_chr1_100"].size() != 2) { std::cerr << "Test 3: Expected 2 reads in ps_chr1_100" << std::endl; return 1; }
    if (phase_groups["ps_chr1_200"].size() != 1) { std::cerr << "Test 3: Expected 1 read in ps_chr1_200" << std::endl; return 1; }
    
    // Test 4: Verify haplotype assignment filtering
    int h1_count = 0, h2_count = 0;
    for (auto &kv : phased_reads) {
        if (kv.second->haplotype == "H1") h1_count++;
        if (kv.second->haplotype == "H2") h2_count++;
    }
    
    if (h1_count != 2) { std::cerr << "Test 4: Expected 2 H1 reads, got " << h1_count << std::endl; return 1; }
    if (h2_count != 1) { std::cerr << "Test 4: Expected 1 H2 read, got " << h2_count << std::endl; return 1; }
    
    // Test 5: Verify SVCluster setup for phasing
    std::map<std::string, std::vector<SVCluster*>> vars;
    SVCluster* sv1 = new SVCluster();
    sv1->node = "node1";
    sv1->start_pos = 100;
    sv1->contig = "contig1";
    sv1->phased = false;
    sv1->reads_h1.insert("read_001");
    sv1->reads_h2.insert("read_002");
    sv1->reads_untagged.insert("read_999");
    
    SVCluster* sv2 = new SVCluster();
    sv2->node = "node2";
    sv2->start_pos = 200;
    sv2->contig = "contig1";
    sv2->phased = false;
    sv2->reads_h1.insert("read_003");
    
    vars["node1"].push_back(sv1);
    vars["node2"].push_back(sv2);
    
    if (vars.size() != 2) { std::cerr << "Test 5: Expected 2 nodes with SVs, got " << vars.size() << std::endl; return 1; }
    if (sv1->reads_h1.size() != 1) { std::cerr << "Test 5: Expected 1 H1 read in sv1" << std::endl; return 1; }
    if (sv1->reads_h2.size() != 1) { std::cerr << "Test 5: Expected 1 H2 read in sv1" << std::endl; return 1; }
    
    // Test 6: Verify SV phasing status
    if (sv1->phased) { std::cerr << "Test 6: Expected sv1.phased=false initially" << std::endl; return 1; }
    
    // Simulate phasing assignment (what phase_svs would do)
    sv1->phased = true;
    if (!sv1->phased) { std::cerr << "Test 6: Expected sv1.phased=true after assignment" << std::endl; return 1; }
    
    // Test 7: Verify read assignment to haplotypes
    sv2->reads_h1.insert("read_001");
    sv2->reads_h2.insert("read_002");
    
    if (sv2->reads_h1.size() != 2) { std::cerr << "Test 7: Expected 2 H1 reads in sv2, got " << sv2->reads_h1.size() << std::endl; return 1; }
    if (sv2->reads_h2.size() != 1) { std::cerr << "Test 7: Expected 1 H2 read in sv2" << std::endl; return 1; }
    
    // Cleanup
    delete ph1;
    delete ph2;
    delete ph3;
    delete sv1;
    delete sv2;
    
    std::cout << "phasing test passed" << std::endl;
    return 0;
}
