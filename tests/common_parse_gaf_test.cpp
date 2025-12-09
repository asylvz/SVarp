#include <iostream>
#include <string>
#include "common.h"

int main() {
    Gaf g;
    std::string line = "read1\t100\t0\t100\t+\tpathA\t1000\t10\t90\t80\t100\t60\ttp:A:P\tcg:Z:10M5I10M\tAS:f:123.4";

    int res = parse_gaf_line(line, g);
    if (res != RETURN_SUCCESS) {
        std::cerr << "parse_gaf_line returned error" << std::endl;
        return 1;
    }

    if (g.query_name != "read1") { std::cerr << "query_name mismatch: " << g.query_name << std::endl; return 1; }
    if (g.query_length != 100) { std::cerr << "query_length mismatch" << std::endl; return 1; }
    if (g.query_start != 0) { std::cerr << "query_start mismatch" << std::endl; return 1; }
    if (g.query_end != 100) { std::cerr << "query_end mismatch" << std::endl; return 1; }
    if (g.strand != "+") { std::cerr << "strand mismatch" << std::endl; return 1; }
    if (g.path != "pathA") { std::cerr << "path mismatch" << std::endl; return 1; }
    if (g.path_length != 1000) { std::cerr << "path_length mismatch" << std::endl; return 1; }
    if (g.path_start != 10 || g.path_end != 90) { std::cerr << "path coords mismatch" << std::endl; return 1; }
    if (g.mapping_quality != 60) { std::cerr << "mapq mismatch" << std::endl; return 1; }
    if (g.cigar != "10M5I10M") { std::cerr << "cigar mismatch: " << g.cigar << std::endl; return 1; }
    if (abs(g.aln_score - 123.4f) > 0.001) { std::cerr << "aln_score mismatch: " << g.aln_score << std::endl; return 1; }

    std::cout << "parse_gaf_line test passed" << std::endl;
    return 0;
}
