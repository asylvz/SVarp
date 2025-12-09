#include <iostream>
#include <vector>
#include <string>
#include "common.h"

int main() {
    std::vector<int> cigarLen;
    std::vector<char> cigarOp;
    int cnt = decompose_cigars(std::string("10M5I3D"), cigarLen, cigarOp);
    if (cnt != 3) { std::cerr << "Expected 3 ops, got " << cnt << std::endl; return 1; }
    if (cigarLen.size() != 3 || cigarOp.size() != 3) { std::cerr << "Sizes mismatch" << std::endl; return 1; }
    if (cigarLen[0] != 10 || cigarLen[1] != 5 || cigarLen[2] != 3) { std::cerr << "Lens mismatch" << std::endl; return 1; }
    if (cigarOp[0] != 'M' || cigarOp[1] != 'I' || cigarOp[2] != 'D') { std::cerr << "Ops mismatch" << std::endl; return 1; }
    std::cout << "decompose_cigars test passed" << std::endl;
    return 0;
}
