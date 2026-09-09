// generate_sv_node vs a reference path walk on both strands; each event must match from the reverse-complement path too.
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <cstdlib>
#include "variant.h"
#include "reference.h"
#include "common.h"

struct Step { char strand; std::string node; int len; };

// path position p (0-based) -> node index and forward coordinate (0-based)
static std::pair<int,int> ref_map(const std::vector<Step>& path, int p)
{
    int acc = 0;
    for (size_t i = 0; i < path.size(); i++) {
        const Step& s = path[i];
        if (p < acc + s.len) {
            int off = p - acc;
            return {(int)i, (s.strand == '>') ? off : (s.len - 1 - off)};
        }
        acc += s.len;
    }
    return {-1, -1};
}

static std::string path_string(const std::vector<Step>& path)
{
    std::string s;
    for (auto& st : path) { s += st.strand; s += st.node; }
    return s;
}

static std::vector<Step> reverse_complement(const std::vector<Step>& path)
{
    std::vector<Step> rc(path.rbegin(), path.rend());
    for (auto& s : rc) s.strand = (s.strand == '>') ? '<' : '>';
    return rc;
}

int main()
{
    std::map<std::string, gfaNode*> gfa;
    std::vector<int> lens = {40, 25, 60};
    std::vector<std::string> names = {"n1", "n2", "n3"};
    int off = 0;
    for (int i = 0; i < 3; i++) { gfa[names[i]] = new gfaNode(names[i], std::string(lens[i], 'A'), lens[i], "c", off); off += lens[i]; }

    std::vector<std::vector<char>> strand_sets = {{'>','>','>'}, {'<','<','<'}, {'>','<','>'}, {'<','>','<'}, {'>'}, {'<'}};
    int checked = 0, bad = 0, spanning = 0;
    for (auto& ss : strand_sets) {
        std::vector<Step> path;
        for (size_t i = 0; i < ss.size(); i++) path.push_back({ss[i], names[i], lens[i]});
        std::vector<Step> rc = reverse_complement(path);
        std::string pstr = path_string(path), rcstr = path_string(rc);
        int plen = 0; for (auto& s : path) plen += s.len;
        for (int path_start : {0, 3}) {
            Gaf g; g.query_name = "r"; g.path = pstr; g.path_length = plen; g.path_start = path_start; g.path_end = plen;
            Gaf gr = g; gr.path = rcstr; gr.path_start = 0;
            for (int var_len : {5, 30, 7}) {
                char type = (var_len == 7) ? INSERTION : DELETION;   // 30 bp deletions cover n2 entirely
                // base_pos: path bases consumed before the op; add_variant passes base_pos + 1
                for (int base_pos = 0; base_pos + (type == DELETION ? var_len : 0) <= plen - path_start; base_pos += 4) {
                    int p0 = path_start + base_pos;                              // first deleted base / base after the insertion
                    int p1 = (type == DELETION) ? p0 + var_len - 1 : p0;         // last deleted base
                    auto m0 = ref_map(path, p0), m1 = ref_map(path, p1);
                    if (type == INSERTION && p0 == plen) continue;              // insertion after the path end: no anchor base
                    Variant* v = generate_sv_node(gfa, g, base_pos + 1, var_len, type);
                    checked++;
                    if (!v) { bad++; std::cout << "NULL path=" << pstr << " ps=" << path_start << " type=" << type << " base_pos=" << base_pos << "\n"; continue; }

                    std::string why;
                    if (type == INSERTION) {
                        const Step& s = path[m0.first];
                        int exp = (s.strand == '>') ? m0.second : m0.second + 1;   // boundary before the forward base after the insertion
                        if (v->node != s.node || v->pos_in_node != exp || v->pos_in_node_end != exp) why = "insertion: expected " + s.node + ":" + std::to_string(exp);
                    } else if (m0.first == m1.first) {
                        const Step& s = path[m0.first];
                        int exp = std::min(m0.second, m1.second);
                        if (v->node != s.node || v->pos_in_node != exp || v->pos_in_node_end != exp + var_len) why = "deletion: expected " + s.node + ":" + std::to_string(exp);
                    } else {
                        // across a boundary: one of the two end nodes, position inside the deleted part of that node
                        spanning++;
                        const Step& s0 = path[m0.first]; const Step& s1 = path[m1.first];
                        bool ok = false;
                        if (v->node == s0.node) ok = (v->pos_in_node == ((s0.strand == '>') ? m0.second : 0));
                        else if (v->node == s1.node) ok = (v->pos_in_node == ((s1.strand == '<') ? m1.second : 0));
                        if (!ok || v->pos_in_node_end != v->pos_in_node + var_len) why = "spanning deletion: expected " + s0.node + " or " + s1.node;
                    }

                    // the same event on the reverse-complement path
                    int rp0 = (type == DELETION) ? plen - 1 - p1 : plen - p0;
                    Variant* w = (rp0 < plen) ? generate_sv_node(gfa, gr, rp0 + 1, var_len, type) : nullptr;
                    if (!w) { if (rp0 < plen) why += " | rc: NULL"; }
                    else if (w->node != v->node || w->pos_in_node != v->pos_in_node || w->pos_in_node_end != v->pos_in_node_end) {
                        // an insertion exactly on a node boundary is anchored to the next base in path order on either strand
                        bool boundary_ins = (type == INSERTION && m0.second == ((path[m0.first].strand == '>') ? 0 : lens[m0.first] - 1));
                        if (!boundary_ins) why += " | rc: got " + w->node + ":" + std::to_string(w->pos_in_node);
                    }
                    if (!why.empty()) {
                        bad++;
                        std::cout << "MISMATCH path=" << pstr << " ps=" << path_start << " type=" << type << " base_pos=" << base_pos
                                  << " | got " << v->node << ":" << v->pos_in_node << "-" << v->pos_in_node_end << " | " << why << "\n";
                    }
                    delete v; delete w;
                }
            }
        }
    }
    std::cout << "checked=" << checked << " spanning=" << spanning << " mismatches=" << bad << "\n";
    for (auto& kv : gfa) delete kv.second;
    if (bad == 0) std::cout << "generate_sv_node strand test passed" << std::endl;
    return bad ? 1 : 0;
}
