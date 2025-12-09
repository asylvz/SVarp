#include <iostream>
#include <string>
#include "common.h"

int main() {
    // overlap_ratio
    double r = overlap_ratio(0, 10, 5, 15);
    if (abs(r - 0.5) > 1e-6) { std::cerr << "overlap_ratio mismatch: " << r << std::endl; return 1; }

    // reverse_complement
    std::string s = "ACTG";
    std::string rc = reverse_complement(s);
    if (rc != "CAGT") { std::cerr << "reverse_complement mismatch: " << rc << std::endl; return 1; }

    std::cout << "common_utils test passed" << std::endl;
    return 0;
}
