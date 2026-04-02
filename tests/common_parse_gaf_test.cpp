#include <iostream>
#include <string>
#include <cmath>
#include "common.h"

int main() {
    // Test 1: Standard GAF line with all fields
    {
        Gaf g;
        std::string line = "read1\t100\t0\t100\t+\tpathA\t1000\t10\t90\t80\t100\t60\ttp:A:P\tcg:Z:10M5I10M\tAS:f:123.4";
        int res = parse_gaf_line(line, g);
        if (res != RETURN_SUCCESS) { std::cerr << "Test 1: parse failed" << std::endl; return 1; }
        if (g.query_name != "read1") { std::cerr << "Test 1: query_name mismatch" << std::endl; return 1; }
        if (g.query_length != 100) { std::cerr << "Test 1: query_length mismatch" << std::endl; return 1; }
        if (g.query_start != 0) { std::cerr << "Test 1: query_start mismatch" << std::endl; return 1; }
        if (g.query_end != 100) { std::cerr << "Test 1: query_end mismatch" << std::endl; return 1; }
        if (g.strand != "+") { std::cerr << "Test 1: strand mismatch" << std::endl; return 1; }
        if (g.path != "pathA") { std::cerr << "Test 1: path mismatch" << std::endl; return 1; }
        if (g.path_length != 1000) { std::cerr << "Test 1: path_length mismatch" << std::endl; return 1; }
        if (g.path_start != 10 || g.path_end != 90) { std::cerr << "Test 1: path coords mismatch" << std::endl; return 1; }
        if (g.mapping_quality != 60) { std::cerr << "Test 1: mapq mismatch" << std::endl; return 1; }
        if (g.cigar != "10M5I10M") { std::cerr << "Test 1: cigar mismatch" << std::endl; return 1; }
        if (fabs(g.aln_score - 123.4f) > 0.1) { std::cerr << "Test 1: aln_score mismatch" << std::endl; return 1; }
        if (!g.is_primary) { std::cerr << "Test 1: should be primary" << std::endl; return 1; }
    }

    // Test 2: Too few fields should fail
    {
        Gaf g;
        std::string line = "read1\t100\t0\t100\t+";
        int res = parse_gaf_line(line, g);
        if (res != RETURN_ERROR) { std::cerr << "Test 2: Expected error for short line" << std::endl; return 1; }
    }

    // Test 3: Secondary alignment (tp:A:S)
    {
        Gaf g;
        std::string line = "read2\t200\t10\t190\t-\t>s1>s2\t500\t0\t480\t170\t180\t30\ttp:A:S\tcg:Z:180M";
        int res = parse_gaf_line(line, g);
        if (res != RETURN_SUCCESS) { std::cerr << "Test 3: parse failed" << std::endl; return 1; }
        if (g.is_primary) { std::cerr << "Test 3: should be secondary" << std::endl; return 1; }
        if (g.strand != "-") { std::cerr << "Test 3: strand should be -" << std::endl; return 1; }
    }

    // Test 4: Line without optional tags (no cigar, no aln_score)
    {
        Gaf g;
        std::string line = "read3\t300\t5\t295\t+\t>node1\t300\t0\t290\t280\t290\t50";
        int res = parse_gaf_line(line, g);
        if (res != RETURN_SUCCESS) { std::cerr << "Test 4: parse failed" << std::endl; return 1; }
        if (g.cigar != "") { std::cerr << "Test 4: cigar should be empty" << std::endl; return 1; }
        if (g.aln_score != 0.0f) { std::cerr << "Test 4: aln_score should be 0" << std::endl; return 1; }
        if (!g.is_primary) { std::cerr << "Test 4: should be primary by default" << std::endl; return 1; }
    }

    // Test 5: Query name with space (should truncate at space)
    {
        Gaf g;
        std::string line = "read4 extra_info\t150\t0\t150\t+\t>s1\t200\t0\t150\t140\t150\t40";
        int res = parse_gaf_line(line, g);
        if (res != RETURN_SUCCESS) { std::cerr << "Test 5: parse failed" << std::endl; return 1; }
        if (g.query_name != "read4") { std::cerr << "Test 5: Expected 'read4', got '" << g.query_name << "'" << std::endl; return 1; }
    }

    // Test 6: Multi-node path
    {
        Gaf g;
        std::string line = "read5\t5000\t100\t4900\t+\t>s100>s101>s102>s103\t6000\t50\t5950\t4700\t4800\t60\tcg:Z:1000M100I3000M500D400M";
        int res = parse_gaf_line(line, g);
        if (res != RETURN_SUCCESS) { std::cerr << "Test 6: parse failed" << std::endl; return 1; }
        if (g.path != ">s100>s101>s102>s103") { std::cerr << "Test 6: path mismatch" << std::endl; return 1; }
        if (g.cigar != "1000M100I3000M500D400M") { std::cerr << "Test 6: cigar mismatch" << std::endl; return 1; }
    }

    // Test 7: Empty line should fail
    {
        Gaf g;
        std::string line = "";
        int res = parse_gaf_line(line, g);
        if (res != RETURN_ERROR) { std::cerr << "Test 7: Expected error for empty line" << std::endl; return 1; }
    }

    std::cout << "parse_gaf_line test passed" << std::endl;
    return 0;
}
