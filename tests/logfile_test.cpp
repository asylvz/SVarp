#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include "common.h"

int main() {
    parameters params;
    // Ensure build directory exists
    std::filesystem::create_directories("build/test_tmp/");
    std::string logpath = "build/test_tmp/";
    params.log_path = logpath;
    params.sample_name = "logfile_test";

    // Open the log file
    std::string logfilename = params.log_path + params.sample_name + ".log";
    if (!params.fp_logs.open(logfilename)) {
        std::cerr << "ERROR: Could not open log file: " << logfilename << std::endl;
        return 1;
    }

    // Write a line and flush using std::endl
    params.fp_logs << "HelloTestLine" << std::endl;

    // Use fp_svtigs as well
    std::string svtig_file = params.log_path + "test_svtigs.fa";
    params.fp_svtigs.open(svtig_file);
    params.fp_svtigs << ">test_svtig" << std::endl;
    params.fp_svtigs << "ACGT" << std::endl;
    params.fp_svtigs.close();

    // Close logs to ensure flush/destruction
    params.fp_logs.close();

    // Verify logfile exists and contains the line
    std::ifstream f(logfilename);
    if (!f) {
        std::cerr << "ERROR: Cannot open logfile " << logfilename << std::endl;
        return 1;
    }

    std::string content; bool found = false;
    while (std::getline(f, content)) {
        if (content.find("HelloTestLine") != std::string::npos) {
            found = true;
            break;
        }
    }
    f.close();

    if (!found) {
        std::cerr << "ERROR: Expected text 'HelloTestLine' not found in log file." << std::endl;
        return 1;
    }

    // Verify svtig file
    std::ifstream s(svtig_file);
    if (!s) {
        std::cerr << "ERROR: Cannot open svtig file " << svtig_file << std::endl;
        return 1;
    }
    bool foundFA = false;
    while (std::getline(s, content)) {
        if (content == ">test_svtig") {
            foundFA = true; break;
        }
    }
    s.close();

    if (!foundFA) {
        std::cerr << "ERROR: Expected FASTA header not found in svtig file." << std::endl;
        return 1;
    }

    // Cleanup
    std::filesystem::remove(logfilename);
    std::filesystem::remove(svtig_file);
    std::filesystem::remove_all(params.log_path);

    std::cout << "LogFile test passed" << std::endl;
    return 0;
}
