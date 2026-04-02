#include <iostream>
#include <string>
#include <cmath>
#include "common.h"

int main() {
    // Test 1: Basic overlap (50%)
    {
        double r = overlap_ratio(0, 10, 5, 15);
        if (fabs(r - 0.5) > 1e-6) { std::cerr << "Test 1: overlap_ratio mismatch: " << r << std::endl; return 1; }
    }

    // Test 2: No overlap
    {
        double r = overlap_ratio(0, 10, 20, 30);
        if (fabs(r) > 1e-6) { std::cerr << "Test 2: Expected 0 overlap, got " << r << std::endl; return 1; }
    }

    // Test 3: Identical ranges (100% overlap)
    {
        double r = overlap_ratio(5, 15, 5, 15);
        if (fabs(r - 1.0) > 1e-6) { std::cerr << "Test 3: Expected 1.0 overlap, got " << r << std::endl; return 1; }
    }

    // Test 4: One range completely contains the other
    {
        double r = overlap_ratio(0, 100, 20, 30);
        // overlap=10, x_len=100, y_len=10, ratio should be max(2*10/110, 10/100, 10/10) = max(0.18, 0.1, 1.0) = 1.0
        if (fabs(r - 1.0) > 1e-6) { std::cerr << "Test 4: Expected 1.0, got " << r << std::endl; return 1; }
    }

    // Test 5: Adjacent ranges (touching but not overlapping)
    {
        double r = overlap_ratio(0, 10, 10, 20);
        if (fabs(r) > 1e-6) { std::cerr << "Test 5: Expected 0 for adjacent, got " << r << std::endl; return 1; }
    }

    // Test 6: Basic reverse complement
    {
        std::string s = "ACTG";
        std::string rc = reverse_complement(s);
        if (rc != "CAGT") { std::cerr << "Test 6: reverse_complement mismatch: " << rc << std::endl; return 1; }
    }

    // Test 7: Reverse complement of single base
    {
        std::string s = "A";
        std::string rc = reverse_complement(s);
        if (rc != "T") { std::cerr << "Test 7: Expected T, got " << rc << std::endl; return 1; }
    }

    // Test 8: Reverse complement of palindrome
    {
        std::string s = "AATT";
        std::string rc = reverse_complement(s);
        if (rc != "AATT") { std::cerr << "Test 8: Expected AATT, got " << rc << std::endl; return 1; }
    }

    // Test 9: Reverse complement double application returns original
    {
        std::string s = "ACGTACGT";
        std::string original = s;
        reverse_complement(s);
        reverse_complement(s);
        if (s != original) { std::cerr << "Test 9: Double reverse_complement should return original, got " << s << std::endl; return 1; }
    }

    // Test 10: Longer sequence
    {
        std::string s = "AAACCCGGGTTT";
        std::string rc = reverse_complement(s);
        if (rc != "AAACCCGGGTTT") { std::cerr << "Test 10: Expected AAACCCGGGTTT, got " << rc << std::endl; return 1; }
    }

    std::cout << "common_utils test passed" << std::endl;
    return 0;
}
