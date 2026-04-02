#include <iostream>
#include <vector>
#include <string>
#include "common.h"

int main() {
    // Test 1: Standard mixed CIGAR
    {
        std::vector<int> cigarLen;
        std::vector<char> cigarOp;
        int cnt = decompose_cigars(std::string("10M5I3D"), cigarLen, cigarOp);
        if (cnt != 3) { std::cerr << "Test 1: Expected 3 ops, got " << cnt << std::endl; return 1; }
        if (cigarLen[0] != 10 || cigarLen[1] != 5 || cigarLen[2] != 3) { std::cerr << "Test 1: Lens mismatch" << std::endl; return 1; }
        if (cigarOp[0] != 'M' || cigarOp[1] != 'I' || cigarOp[2] != 'D') { std::cerr << "Test 1: Ops mismatch" << std::endl; return 1; }
    }

    // Test 2: Single operation
    {
        std::vector<int> cigarLen;
        std::vector<char> cigarOp;
        int cnt = decompose_cigars(std::string("100M"), cigarLen, cigarOp);
        if (cnt != 1) { std::cerr << "Test 2: Expected 1 op, got " << cnt << std::endl; return 1; }
        if (cigarLen[0] != 100 || cigarOp[0] != 'M') { std::cerr << "Test 2: Mismatch" << std::endl; return 1; }
    }

    // Test 3: Large numbers
    {
        std::vector<int> cigarLen;
        std::vector<char> cigarOp;
        int cnt = decompose_cigars(std::string("50000M200I30000D"), cigarLen, cigarOp);
        if (cnt != 3) { std::cerr << "Test 3: Expected 3 ops, got " << cnt << std::endl; return 1; }
        if (cigarLen[0] != 50000) { std::cerr << "Test 3: Expected 50000, got " << cigarLen[0] << std::endl; return 1; }
        if (cigarLen[1] != 200) { std::cerr << "Test 3: Expected 200, got " << cigarLen[1] << std::endl; return 1; }
        if (cigarLen[2] != 30000) { std::cerr << "Test 3: Expected 30000, got " << cigarLen[2] << std::endl; return 1; }
    }

    // Test 4: All CIGAR operation types
    {
        std::vector<int> cigarLen;
        std::vector<char> cigarOp;
        int cnt = decompose_cigars(std::string("5M3I2D1N4S6H7X"), cigarLen, cigarOp);
        if (cnt != 7) { std::cerr << "Test 4: Expected 7 ops, got " << cnt << std::endl; return 1; }
        if (cigarOp[0] != 'M' || cigarOp[1] != 'I' || cigarOp[2] != 'D' ||
            cigarOp[3] != 'N' || cigarOp[4] != 'S' || cigarOp[5] != 'H' || cigarOp[6] != 'X') {
            std::cerr << "Test 4: Ops mismatch" << std::endl; return 1;
        }
    }

    // Test 5: Empty CIGAR string
    {
        std::vector<int> cigarLen;
        std::vector<char> cigarOp;
        int cnt = decompose_cigars(std::string(""), cigarLen, cigarOp);
        if (cnt != 0) { std::cerr << "Test 5: Expected 0 ops for empty, got " << cnt << std::endl; return 1; }
    }

    // Test 6: Long complex CIGAR (realistic)
    {
        std::vector<int> cigarLen;
        std::vector<char> cigarOp;
        int cnt = decompose_cigars(std::string("150M50I200M75D100M"), cigarLen, cigarOp);
        if (cnt != 5) { std::cerr << "Test 6: Expected 5 ops, got " << cnt << std::endl; return 1; }
        if (cigarLen[0] != 150 || cigarLen[1] != 50 || cigarLen[2] != 200 ||
            cigarLen[3] != 75 || cigarLen[4] != 100) {
            std::cerr << "Test 6: Lens mismatch" << std::endl; return 1;
        }
    }

    std::cout << "decompose_cigars test passed" << std::endl;
    return 0;
}
