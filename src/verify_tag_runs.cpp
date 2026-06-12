// verify_tag_runs: for each tag run in a tag array, verify that every BWT
// position inside the run resolves (via r-index locate + gafpack-style
// cum_bp walk on the GFA path) to the same (node, offset) that the tag array
// stores for that run.
//
// Purpose: diagnose tag-array vs GFA inconsistencies that could cause the
// lightweight pipeline (one record per tag run) to produce different
// (read, read_st, node, offset) tuples than the all-positions pipeline
// (one record per BWT position) after gafpack's dedup.
//
// Algorithm (per tag run):
//   1. Read stored (node_T, offset_T, strand_T) and BWT interval [bwt_lo, bwt_hi].
//   2. SA = locate_sa_value(bwt_lo).
//   3. For each BWT position p in [bwt_lo, bwt_hi]:
//        seq_id  = r_index.seqId(SA), path_bp = r_index.seqOffset(SA)
//        path_idx = seq_id / 2, rev_bucket = seq_id & 1
//        Walk the (reversed if rev_bucket) step list to find the step
//        containing path_bp -> (walked_node, walked_offset).
//        Compare against (node_T, offset_T); count mismatches.
//        SA = locateNext(SA) (skipped after last iter).
//
// Strand intentionally NOT compared: per find_mems.cpp:35-39 the tag's
// is_rev(graph_pos) and the GAF strand (bucket-parity XOR step-orientation)
// are independent concepts.
//
// CLI:
//   verify_tag_runs <r_index> <tag_array> <gfa> <paths_file>
//                   [--sample N | --all]
//                   [--max-mismatch-dump K]  (default 50)
//                   [--seed S]               (default 42)

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <random>
#include <atomic>
#include <mutex>
#include <chrono>
#include <algorithm>
#include <cstdlib>

#include "pangenome_index/r-index.hpp"
#include "pangenome_index/tag_arrays.hpp"
#include <gbwtgraph/utils.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace panindexer;
using namespace gbwtgraph;

// ---------- locate helper (mirrors find_mems.cpp:locate_sa_value) ----------
static inline size_type locate_sa_value(const FastLocate& r_index, size_t bwt_pos) {
    size_t run_id = 0, offset_of_first = 0;
    r_index.run_id_and_offset_at(bwt_pos, run_id, offset_of_first);
    size_type sa = r_index.getSample(run_id);
    while (offset_of_first < bwt_pos) {
        sa = r_index.locateNext(sa);
        offset_of_first++;
    }
    return sa;
}

// ---------- per-path layout (one per GFA path, indexed by seq_id/2) ----------
struct PathLayout {
    string name;
    vector<uint32_t> step_node_ids;  // forward step order
    vector<uint64_t> cum_bp;         // size = n_steps + 1, cum_bp[0] = 0
};

// ---------- mismatch record (kept small to bound memory) --------------------
struct Mismatch {
    uint64_t run_id;
    uint64_t bwt_pos;
    uint64_t seq_id;
    uint64_t path_bp;
    uint64_t tag_node;
    uint32_t tag_offset;
    uint64_t walked_node;
    uint64_t walked_offset;
    bool path_bp_out_of_range;  // true if path_bp >= path_total_bp
};

// =============================================================================
// GFA parsing: segment lengths + per-path step lists
// =============================================================================

static void load_segment_lengths(const string& gfa_path,
                                 unordered_map<uint64_t, uint32_t>& node_len) {
    ifstream in(gfa_path);
    if (!in) { cerr << "ERROR: cannot open GFA: " << gfa_path << endl; exit(1); }
    string line;
    size_t loaded = 0;
    while (getline(in, line)) {
        if (line.empty() || line[0] != 'S') continue;
        size_t tab1 = line.find('\t');
        if (tab1 == string::npos) continue;
        size_t tab2 = line.find('\t', tab1 + 1);
        if (tab2 == string::npos) continue;
        size_t tab3 = line.find('\t', tab2 + 1);
        if (tab3 == string::npos) tab3 = line.size();
        uint64_t nid = strtoull(line.c_str() + tab1 + 1, nullptr, 10);
        uint32_t slen = static_cast<uint32_t>(tab3 - (tab2 + 1));
        node_len[nid] = slen;
        if (++loaded % 1000000 == 0) {
            cerr << "  loaded " << loaded << " segment lengths" << endl;
        }
    }
    cerr << "  total segments: " << node_len.size() << endl;
}

// Read paths_file into an in-order vector of names (line N -> path_idx N)
static vector<string> read_path_names(const string& paths_file) {
    ifstream pf(paths_file);
    if (!pf) { cerr << "ERROR: cannot open paths file: " << paths_file << endl; exit(1); }
    vector<string> names;
    string line;
    while (getline(pf, line)) {
        if (!line.empty()) names.push_back(line);
    }
    cerr << "  paths file has " << names.size() << " path names" << endl;
    return names;
}

// Parse GFA P-lines, populate per-path step+cum_bp for paths present in
// path_layouts (indexed by the path_names_file order).
static void load_path_layouts(const string& gfa_path,
                              const vector<string>& path_names,
                              const unordered_map<uint64_t, uint32_t>& node_len,
                              vector<PathLayout>& layouts) {
    // Build name -> idx map for fast lookup; use it once.
    unordered_map<string, size_t> name_to_idx;
    name_to_idx.reserve(path_names.size() * 2);
    for (size_t i = 0; i < path_names.size(); i++) {
        name_to_idx.emplace(path_names[i], i);
        layouts[i].name = path_names[i];
    }

    ifstream in(gfa_path);
    if (!in) { cerr << "ERROR: cannot open GFA: " << gfa_path << endl; exit(1); }
    string line;
    size_t matched = 0;
    while (getline(in, line)) {
        if (line.empty() || line[0] != 'P') continue;
        // P<TAB>name<TAB>steps[<TAB>overlaps]
        size_t tab1 = line.find('\t');
        if (tab1 == string::npos) continue;
        size_t tab2 = line.find('\t', tab1 + 1);
        if (tab2 == string::npos) continue;
        size_t tab3 = line.find('\t', tab2 + 1);
        if (tab3 == string::npos) tab3 = line.size();

        string name = line.substr(tab1 + 1, tab2 - tab1 - 1);
        auto it = name_to_idx.find(name);
        if (it == name_to_idx.end()) continue;
        size_t idx = it->second;

        // Parse steps: comma-separated "<id>+" / "<id>-"
        PathLayout& L = layouts[idx];
        L.cum_bp.push_back(0);
        size_t pos = tab2 + 1;
        size_t end = tab3;
        size_t step_start = pos;
        while (step_start < end) {
            size_t step_end = step_start;
            while (step_end < end && line[step_end] != ',') step_end++;
            // Step format: digits then a trailing + or -
            if (step_end > step_start + 1) {
                size_t digit_end = step_end - 1;  // exclude the +/- char
                uint64_t nid = strtoull(line.c_str() + step_start, nullptr, 10);
                (void)digit_end;
                auto nl_it = node_len.find(nid);
                uint32_t nlen = (nl_it != node_len.end()) ? nl_it->second : 0;
                L.step_node_ids.push_back(static_cast<uint32_t>(nid));
                L.cum_bp.push_back(L.cum_bp.back() + nlen);
            }
            step_start = step_end + 1;
        }
        matched++;
        if (matched % 100 == 0) {
            cerr << "  parsed " << matched << " path P-lines" << endl;
        }
    }
    cerr << "  matched " << matched << " of " << path_names.size() << " expected paths" << endl;
}

// =============================================================================
// gafpack-style cum_bp walk: find step index containing path_bp
// =============================================================================

// Returns step_idx, and writes the local offset to *out_off.
// path_bp must be < cum_bp.back() (caller's responsibility).
static inline size_t step_for_path_bp(const vector<uint64_t>& cum_bp,
                                      uint64_t path_bp,
                                      uint64_t& out_off) {
    // Binary search: find largest i such that cum_bp[i] <= path_bp.
    // cum_bp is strictly non-decreasing.
    size_t lo = 0, hi = cum_bp.size() - 1;  // last entry is path_total_bp
    while (lo + 1 < hi) {
        size_t mid = (lo + hi) >> 1;
        if (cum_bp[mid] <= path_bp) lo = mid;
        else hi = mid;
    }
    out_off = path_bp - cum_bp[lo];
    return lo;
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char** argv) {
    if (argc < 5) {
        cerr << "Usage: " << argv[0]
             << " <r_index> <tag_array> <gfa> <paths_file>"
             << " [--sample N | --all] [--max-mismatch-dump K] [--seed S]" << endl;
        return 1;
    }

    string ri_path    = argv[1];
    string tags_path  = argv[2];
    string gfa_path   = argv[3];
    string paths_file = argv[4];

    bool   do_all          = false;
    size_t sample_n        = 100000;
    size_t max_mismatch    = 50;
    uint64_t seed          = 42;

    for (int i = 5; i < argc; i++) {
        string a = argv[i];
        if (a == "--all") do_all = true;
        else if (a == "--sample" && i + 1 < argc) { sample_n = strtoull(argv[++i], nullptr, 10); }
        else if (a == "--max-mismatch-dump" && i + 1 < argc) { max_mismatch = strtoull(argv[++i], nullptr, 10); }
        else if (a == "--seed" && i + 1 < argc) { seed = strtoull(argv[++i], nullptr, 10); }
        else { cerr << "ERROR: unknown arg: " << a << endl; return 1; }
    }

    cerr << "# Loading r-index: " << ri_path << endl;
    FastLocate r_index;
    {
        ifstream in(ri_path, ios::binary);
        if (!in) { cerr << "ERROR: cannot open r-index" << endl; return 1; }
        r_index.load_encoded(in);
    }

    cerr << "# Loading tag array: " << tags_path << endl;
    TagArray tag_array;
    {
        ifstream in(tags_path, ios::binary);
        if (!in) { cerr << "ERROR: cannot open tag array" << endl; return 1; }
        tag_array.load_compressed_tags_compact(in);
    }

    cerr << "# Loading segment lengths from GFA: " << gfa_path << endl;
    unordered_map<uint64_t, uint32_t> node_len;
    load_segment_lengths(gfa_path, node_len);

    cerr << "# Reading paths file: " << paths_file << endl;
    vector<string> path_names = read_path_names(paths_file);
    vector<PathLayout> path_layouts(path_names.size());

    cerr << "# Parsing GFA P-lines for path step lists" << endl;
    load_path_layouts(gfa_path, path_names, node_len, path_layouts);

    // Collect tag runs into a vector for indexing + sampling.
    cerr << "# Iterating tag runs from tag array" << endl;
    struct TagRun {
        pos_t pos;
        uint64_t bwt_start;
        uint64_t bwt_end;   // inclusive
        uint64_t run_length;
    };
    vector<TagRun> tag_runs;
    tag_runs.reserve(1u << 20);
    tag_array.for_each_run_compact_with_bwt(
        [&](pos_t p, uint64_t run_length, size_t bwt_start, size_t bwt_end) {
            tag_runs.push_back({p, bwt_start, bwt_end - 1, run_length});
        });
    cerr << "  total tag runs: " << tag_runs.size() << endl;

    // Build sampling indices.
    vector<uint64_t> sample_idx;
    if (do_all || sample_n >= tag_runs.size()) {
        sample_idx.resize(tag_runs.size());
        for (uint64_t i = 0; i < tag_runs.size(); i++) sample_idx[i] = i;
        cerr << "# Mode: --all (" << sample_idx.size() << " runs)" << endl;
    } else {
        sample_idx.reserve(sample_n);
        mt19937_64 rng(seed);
        uniform_int_distribution<uint64_t> dist(0, tag_runs.size() - 1);
        // Reservoir-ish: just pick N uniformly random indices with replacement
        // (collisions rare at N << total); fine for sampling diagnostic.
        for (uint64_t i = 0; i < sample_n; i++) sample_idx.push_back(dist(rng));
        cerr << "# Mode: --sample " << sample_n << " (seed=" << seed << ")" << endl;
    }

    // Per-thread accumulators.
    atomic<uint64_t> runs_checked{0};
    atomic<uint64_t> runs_fully_consistent{0};
    atomic<uint64_t> runs_with_mismatch{0};
    atomic<uint64_t> total_positions_checked{0};
    atomic<uint64_t> total_mismatches{0};
    atomic<uint64_t> path_bp_oob_count{0};

    vector<Mismatch> mismatch_dump;
    mismatch_dump.reserve(max_mismatch);
    mutex dump_mutex;

    auto t0 = chrono::steady_clock::now();

    #pragma omp parallel for schedule(dynamic, 1024)
    for (size_t s = 0; s < sample_idx.size(); s++) {
        const uint64_t run_id = sample_idx[s];
        const TagRun& tr = tag_runs[run_id];
        const uint64_t tag_node   = id(tr.pos);
        const uint32_t tag_offset = static_cast<uint32_t>(offset(tr.pos));

        size_type sa = locate_sa_value(r_index, tr.bwt_start);
        bool any_mismatch = false;
        uint64_t local_oob = 0;

        for (uint64_t p = tr.bwt_start; p <= tr.bwt_end; p++) {
            size_type seq_id  = r_index.seqId(sa);
            size_type path_bp = r_index.seqOffset(sa);

            size_t path_idx = static_cast<size_t>(seq_id) >> 1;
            bool   rev      = (seq_id & 1) != 0;

            uint64_t walked_node = 0;
            uint64_t walked_off  = 0;
            bool oob = false;

            if (path_idx >= path_layouts.size() ||
                path_layouts[path_idx].cum_bp.empty()) {
                // path not loaded -- treat as OOB
                oob = true;
            } else {
                const PathLayout& L = path_layouts[path_idx];
                uint64_t total_bp = L.cum_bp.back();
                if (path_bp >= total_bp) {
                    oob = true;
                    local_oob++;
                } else {
                    uint64_t off_in_step = 0;
                    if (rev) {
                        // For reverse bucket: walk on the reversed step list.
                        // Equivalent: find the step from the END.
                        // The simplest correct equivalence: the reversed-step
                        // cum_bp at position k corresponds to the forward
                        // cum_bp from the end. We just transform path_bp:
                        // path_bp_fwd = total_bp - 1 - path_bp (point-symmetric).
                        // Then the forward step containing that point gives
                        // the same node; the local offset is mirrored within
                        // the node.
                        uint64_t fwd_bp = total_bp - 1 - path_bp;
                        size_t step = step_for_path_bp(L.cum_bp, fwd_bp, off_in_step);
                        walked_node = L.step_node_ids[step];
                        uint64_t node_length = L.cum_bp[step + 1] - L.cum_bp[step];
                        // Mirror offset within node so it's bucket-relative.
                        walked_off = (node_length == 0) ? 0
                                   : (node_length - 1 - off_in_step);
                    } else {
                        size_t step = step_for_path_bp(L.cum_bp, path_bp, off_in_step);
                        walked_node = L.step_node_ids[step];
                        walked_off  = off_in_step;
                    }
                }
            }

            bool match = !oob &&
                         (walked_node == tag_node) &&
                         (walked_off  == static_cast<uint64_t>(tag_offset));
            if (!match) {
                any_mismatch = true;
                total_mismatches.fetch_add(1, memory_order_relaxed);
                // Sample a few full-detail mismatches
                if (mismatch_dump.size() < max_mismatch) {
                    lock_guard<mutex> lk(dump_mutex);
                    if (mismatch_dump.size() < max_mismatch) {
                        mismatch_dump.push_back({
                            run_id, p,
                            static_cast<uint64_t>(seq_id),
                            static_cast<uint64_t>(path_bp),
                            tag_node, tag_offset,
                            walked_node, walked_off,
                            oob
                        });
                    }
                }
            }

            total_positions_checked.fetch_add(1, memory_order_relaxed);
            if (p < tr.bwt_end) sa = r_index.locateNext(sa);
        }

        if (local_oob > 0) path_bp_oob_count.fetch_add(local_oob, memory_order_relaxed);

        runs_checked.fetch_add(1, memory_order_relaxed);
        if (any_mismatch) runs_with_mismatch.fetch_add(1, memory_order_relaxed);
        else              runs_fully_consistent.fetch_add(1, memory_order_relaxed);

        if ((runs_checked.load(memory_order_relaxed) & 0xFFFF) == 0) {
            #pragma omp critical
            {
                cerr << "  progress: " << runs_checked.load()
                     << " / " << sample_idx.size()
                     << "  mismatches=" << total_mismatches.load() << endl;
            }
        }
    }

    auto t1 = chrono::steady_clock::now();
    double secs = chrono::duration<double>(t1 - t0).count();

    // -------- Output --------
    cout << "# === SUMMARY ===" << endl;
    cout << "# r_index:                " << ri_path << endl;
    cout << "# tag_array:              " << tags_path << endl;
    cout << "# gfa:                    " << gfa_path << endl;
    cout << "# paths_file:             " << paths_file << endl;
    cout << "# total_tag_runs:         " << tag_runs.size() << endl;
    cout << "# mode:                   " << (do_all ? "ALL" : "SAMPLE") << endl;
    cout << "# runs_checked:           " << runs_checked.load() << endl;
    cout << "# runs_fully_consistent:  " << runs_fully_consistent.load()
         << "  (" << (100.0 * runs_fully_consistent.load() / max((uint64_t)1, runs_checked.load()))
         << "%)" << endl;
    cout << "# runs_with_mismatch:     " << runs_with_mismatch.load()
         << "  (" << (100.0 * runs_with_mismatch.load() / max((uint64_t)1, runs_checked.load()))
         << "%)" << endl;
    cout << "# total_positions_checked:" << total_positions_checked.load() << endl;
    cout << "# total_mismatches:       " << total_mismatches.load()
         << "  (" << (100.0 * total_mismatches.load() / max((uint64_t)1, total_positions_checked.load()))
         << "% of positions)" << endl;
    cout << "# path_bp_out_of_range:   " << path_bp_oob_count.load() << endl;
    cout << "# wall_seconds:           " << secs << endl;
    cout << endl;

    cout << "# === FIRST " << mismatch_dump.size() << " MISMATCHES ===" << endl;
    cout << "# run_id\tbwt_pos\tseq_id\tpath_bp\ttag_node\ttag_offset\twalked_node\twalked_offset\tpath_bp_oob" << endl;
    for (const auto& m : mismatch_dump) {
        cout << m.run_id << '\t' << m.bwt_pos << '\t' << m.seq_id << '\t'
             << m.path_bp << '\t' << m.tag_node << '\t' << m.tag_offset << '\t'
             << m.walked_node << '\t' << m.walked_offset << '\t'
             << (m.path_bp_out_of_range ? "1" : "0") << '\n';
    }

    return 0;
}