#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <fstream>
#include <cstdlib>
#include "assembly.h"
#include "reference.h"
#include "variant.h"
#include "common.h"

int main() {
    // Test 1: Verify Assembly class instantiation and basic properties
    Assembly assembly_obj;
    if (assembly_obj.filter_hicov != 0) { std::cerr << "Test 1: Expected filter_hicov=0, got " << assembly_obj.filter_hicov << std::endl; return 1; }
    if (assembly_obj.filter_lowcov != 0) { std::cerr << "Test 1: Expected filter_lowcov=0, got " << assembly_obj.filter_lowcov << std::endl; return 1; }
    if (assembly_obj.filter_support != 0) { std::cerr << "Test 1: Expected filter_support=0, got " << assembly_obj.filter_support << std::endl; return 1; }
    if (assembly_obj.unassembled_cnt != 0) { std::cerr << "Test 1: Expected unassembled_cnt=0, got " << assembly_obj.unassembled_cnt << std::endl; return 1; }
    
    // Test 2: Verify raw_svtigs set operations (core Assembly functionality)
    assembly_obj.raw_svtigs.insert("svtig_001");
    assembly_obj.raw_svtigs.insert("svtig_002");
    assembly_obj.raw_svtigs.insert("svtig_003");
    
    if (assembly_obj.raw_svtigs.size() != 3) { std::cerr << "Test 2: Expected 3 raw_svtigs, got " << assembly_obj.raw_svtigs.size() << std::endl; return 1; }
    if (assembly_obj.raw_svtigs.find("svtig_001") == assembly_obj.raw_svtigs.end()) { std::cerr << "Test 2: Expected svtig_001 in raw_svtigs" << std::endl; return 1; }
    if (assembly_obj.raw_svtigs.find("svtig_002") == assembly_obj.raw_svtigs.end()) { std::cerr << "Test 2: Expected svtig_002 in raw_svtigs" << std::endl; return 1; }
    if (assembly_obj.raw_svtigs.find("svtig_003") == assembly_obj.raw_svtigs.end()) { std::cerr << "Test 2: Expected svtig_003 in raw_svtigs" << std::endl; return 1; }
    
    // Test 3: Verify incrementing counters
    assembly_obj.filter_hicov++;
    assembly_obj.filter_lowcov += 2;
    assembly_obj.filter_support += 3;
    assembly_obj.unassembled_cnt++;
    
    if (assembly_obj.filter_hicov != 1) { std::cerr << "Test 3: Expected filter_hicov=1, got " << assembly_obj.filter_hicov << std::endl; return 1; }
    if (assembly_obj.filter_lowcov != 2) { std::cerr << "Test 3: Expected filter_lowcov=2, got " << assembly_obj.filter_lowcov << std::endl; return 1; }
    if (assembly_obj.filter_support != 3) { std::cerr << "Test 3: Expected filter_support=3, got " << assembly_obj.filter_support << std::endl; return 1; }
    if (assembly_obj.unassembled_cnt != 1) { std::cerr << "Test 3: Expected unassembled_cnt=1, got " << assembly_obj.unassembled_cnt << std::endl; return 1; }
    
    // Test 4: Splitting a cluster's reads by haplotype must not change whether
    // the cluster is plausible for the region's coverage.
    // contig_depth is above MAX_CONTIG_DEPTH so both outcomes return before any
    // assembler runs; only which filter fires differs.
    {
        parameters params;
        params.support = 5;              // support threshold becomes 3
        params.log_path = "/tmp/";

        SVCluster cluster;
        cluster.node = "s1";
        cluster.contig = "chr1";
        cluster.ref_pos = 1000;
        cluster.start_pos = 100;
        for (int i = 0; i < 6; i++)  cluster.reads_h1.insert("h1_read_" + std::to_string(i));
        for (int i = 0; i < 24; i++) cluster.reads_h2.insert("h2_read_" + std::to_string(i));

        Assembly asm_phased;
        faidx_t* idx = nullptr;
        std::map<std::string, SVtig*> svtigs;
        std::string name = "H1-s1_100";
        double contig_depth = 150.0;     // 30 cluster reads, so 150 <= 5 * 30
        SVCluster* svp = &cluster;

        asm_phased.final_assembly(params, idx, cluster.reads_h1, name, contig_depth, svp, svtigs);

        if (asm_phased.filter_lowcov != 0) {
            std::cerr << "Test 4: H1 half judged low coverage against the whole contig depth" << std::endl;
            return 1;
        }
        if (asm_phased.filter_hicov != 1) {
            std::cerr << "Test 4: expected the cluster to reach the contig depth check, got hicov="
                      << asm_phased.filter_hicov << std::endl;
            return 1;
        }
        for (auto& kv : svtigs) delete kv.second;
    }

    // Test 5: Without a phase file the whole cluster is the read set, so the
    // same cluster is still judged low coverage.
    {
        parameters params;
        params.support = 5;
        params.log_path = "/tmp/";

        SVCluster cluster;
        cluster.node = "s1";
        cluster.contig = "chr1";
        cluster.ref_pos = 1000;
        cluster.start_pos = 100;
        for (int i = 0; i < 6; i++) cluster.reads_untagged.insert("read_" + std::to_string(i));

        Assembly asm_unphased;
        faidx_t* idx = nullptr;
        std::map<std::string, SVtig*> svtigs;
        std::string name = "s1_100";
        double contig_depth = 150.0;     // only 6 reads, so 150 > 5 * 6
        SVCluster* svp = &cluster;

        asm_unphased.final_assembly(params, idx, cluster.reads_untagged, name, contig_depth, svp, svtigs);

        if (asm_unphased.filter_lowcov != 1) {
            std::cerr << "Test 5: expected low coverage, got lowcov="
                      << asm_unphased.filter_lowcov << std::endl;
            return 1;
        }
        for (auto& kv : svtigs) delete kv.second;
    }

    std::cout << "assembly test passed" << std::endl;
    return 0;
}
