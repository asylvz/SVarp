#include <iostream>
#include <string>
#include <filesystem>
#include <fstream>
#include "common.h"

int main() {
    parameters params;
    std::filesystem::create_directories("build/test_tmp2/");
    params.log_path = "build/test_tmp2/";
    params.sample_name = "run_and_log_test";
    params.fp_logs.open(params.log_path + params.sample_name + ".log");

    // Run a command that fails on most systems
    std::string cmd = "sh -c 'exit 3'";
    int rc = run_and_log(cmd, params, "test-fail", 0, 1, false);
    if (rc == 0) { std::cerr << "Expected non-zero rc from failing command" << std::endl; return 1; }

    params.fp_logs.close();

    // Verify the log contains the failure line
    std::ifstream f(params.log_path + params.sample_name + ".log");
    std::string line; bool found = false;
    while (std::getline(f, line)) {
        if (line.find("test-fail") != std::string::npos) { found = true; break; }
    }
    f.close();
    if (!found) { std::cerr << "Log entry not found for failing command" << std::endl; return 1; }

    std::filesystem::remove_all(params.log_path);

    std::cout << "run_and_log test passed" << std::endl;
    return 0;
}
