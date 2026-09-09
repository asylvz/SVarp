// parse_command_line: defaults and range checks on numeric options
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <fstream>
#include <getopt.h>
#include "cmdline.h"
#include "common.h"

static const std::string DIR = "/tmp/svarp_cmdline_test";

static int parse(std::vector<std::string> extra, parameters& params)
{
    std::filesystem::create_directories(DIR);
    for (const char* f : {"g.gfa", "a.gaf", "r.fa"}) std::ofstream(DIR + "/" + f).put('\n');   // the parser checks that inputs exist
    std::vector<std::string> args = {"svarp", "--graph", DIR + "/g.gfa", "--gaf", DIR + "/a.gaf", "--fasta", DIR + "/r.fa", "--out", DIR};
    args.insert(args.end(), extra.begin(), extra.end());
    std::vector<char*> argv;
    for (auto& a : args) argv.push_back(const_cast<char*>(a.c_str()));
    optind = 1;
#if defined(__APPLE__) || defined(__FreeBSD__)
    optreset = 1;
#endif
    params = parameters();
    return parse_command_line((int)argv.size(), argv.data(), params);
}

int main()
{
    parameters params;
    if (parse({}, params) != RETURN_SUCCESS) { std::cerr << "Test 1: minimal command line rejected" << std::endl; return 1; }
    if (params.support != 5 || params.threads != 16 || params.asm_jobs != 8 || params.dist_threshold != 100 || params.min_clip != 500 || params.min_graph_cov != 0.90 || params.min_svtig_len != 5000 || !params.skip_untagged) {
        std::cerr << "Test 1: unexpected defaults" << std::endl; return 1;
    }
    if (parse({"--support", "3", "--threads", "4", "--min-identity", "0.8", "--min-graph-cov", "0.95", "--pc", "0.9", "--min-svtig-len", "1000", "--keep-untagged"}, params) != RETURN_SUCCESS || params.support != 3 || params.threads != 4 || params.asm_jobs != 4 || params.min_graph_cov != 0.95 || params.min_svtig_len != 1000 || params.skip_untagged) {
        std::cerr << "Test 2: valid values rejected" << std::endl; return 1;
    }
    if (parse({"--keep-remap", "--write-unmapped"}, params) != RETURN_SUCCESS || !params.keep_remap || !params.write_unmapped) {
        std::cerr << "Test 2: flags not set" << std::endl; return 1;
    }
    std::vector<std::vector<std::string>> bad = {
        {"--support", "-1"}, {"--threads", "0"}, {"--dist-threshold", "-5"}, {"--as", "-10"},
        {"--min-identity", "1.5"}, {"--min-graph-cov", "-0.1"}, {"--min-graph-cov", "1.5"}, {"--min-svtig-len", "-5"}, {"--pc", "0"}, {"--pc", "1.2"}, {"--asm-jobs", "0"}, {"--support", "five"},
    };
    for (auto& b : bad) {
        if (parse(b, params) != RETURN_ERROR) { std::cerr << "Test 3: accepted " << b[0] << " " << b[1] << std::endl; return 1; }
    }
    std::filesystem::remove_all(DIR);
    std::cout << "cmdline test passed" << std::endl;
    return 0;
}
