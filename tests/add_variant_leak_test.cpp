// add_variant must not leak the Variant returned by generate_sv_node when a
// repeat SV signal hits a locus already in the map. LeakSanitizer is not
// available on macOS, so global operator new/delete are replaced with a
// counter; the test asserts net live allocations return to zero after cleanup.
//
// The scenario calls add_variant twice for the same (node, position): the first
// stores the Variant, the second adds its read to the existing entry and must
// free the Variant it allocated.

#include <cstdlib>
#include <new>
#include <map>
#include <string>
#include <iostream>
#include "variant.h"
#include "reference.h"
#include "common.h"

static long g_live = 0;

void* operator new(std::size_t n) { g_live++; return std::malloc(n ? n : 1); }
void* operator new[](std::size_t n) { g_live++; return std::malloc(n ? n : 1); }
void  operator delete(void* p) noexcept { if (p) { g_live--; std::free(p); } }
void  operator delete[](void* p) noexcept { if (p) { g_live--; std::free(p); } }
void  operator delete(void* p, std::size_t) noexcept { if (p) { g_live--; std::free(p); } }
void  operator delete[](void* p, std::size_t) noexcept { if (p) { g_live--; std::free(p); } }

// defined in alignment.cpp (not declared in a header)
int add_variant(std::map<std::string, gfaNode*>& gfa,
                std::map<std::string, Variant*>& vars,
                Gaf& line, int pos_in_cigar, char var_type, int var_len);

int main()
{
    std::map<std::string, gfaNode*> gfa;
    gfa.insert({"n1", new gfaNode("n1", "ACGT", 100, "contig1", 0)});

    Gaf line;
    line.query_name = "readX";
    line.query_length = 1000;
    line.query_start = 0;
    line.query_end = 500;
    line.path = ">n1";
    line.path_length = 100;
    line.path_start = 0;
    line.path_end = 80;

    std::map<std::string, Variant*> vars;

    long base = g_live;

    // Same locus twice: the second call must free the Variant it allocates
    add_variant(gfa, vars, line, 10, INSERTION, 60);
    add_variant(gfa, vars, line, 10, INSERTION, 60);

    // Free everything reachable
    for (auto& kv : vars) delete kv.second;
    vars.clear();

    long delta = g_live - base;

    for (auto& kv : gfa) delete kv.second;
    gfa.clear();

    if (delta != 0) {
        std::cerr << "add_variant leaked; net live allocations = " << delta << std::endl;
        return 1;
    }
    std::cout << "add_variant leak test passed" << std::endl;
    return 0;
}
