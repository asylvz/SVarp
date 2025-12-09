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
    
    std::cout << "assembly test passed" << std::endl;
    return 0;
}
