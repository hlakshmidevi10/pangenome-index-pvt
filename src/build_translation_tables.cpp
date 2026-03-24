/**
 * build_translation_tables: build and store Translation Table 1 and Table 2 from a GBZ file.
 *
 * Table 1 maps (named_path, global_interval) → (path_id, local_interval).
 * Table 2 maps (source_path_id, target_haplotype) → sorted IntervalMapping segments.
 *
 * Both tables are built together in two phases so that data computed for Table 1
 * (extracted paths, path lengths, base names) is reused when building Table 2.
 *
 * Usage:
 *   build_translation_tables <graph.gbz> <output.t1> <output.t2> [options]
 *
 * Options:
 *   --dump          Print a human-readable listing of both tables to stderr.
 *   --debug         Print verbose progress to stderr.
 *   --only-table1   Build and store only Table 1 (skip Table 2).
 *   --help          Show this help.
 *
 * See TRANSLATION_TABLES_ALGORITHM.md for full algorithm description.
 */

#include "pangenome_index/translation_tables.hpp"
#include <sdsl/simple_sds.hpp>
#include <gbwt/gbwt.h>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/gbz.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <omp.h>
#include <tuple>

using namespace std;
using namespace std::chrono;

// ============================================================
// Utilities
// ============================================================

static void usage(const char* prog) {
    cerr << "Usage: " << prog << " <graph.gbz> <output.t1> <output.t2> [options]" << endl;
    cerr << endl;
    cerr << "  graph.gbz     GBZ file (GBWT index + GBWTGraph)" << endl;
    cerr << "  output.t1     Binary Translation Table 1 output" << endl;
    cerr << "  output.t2     Binary Translation Table 2 output" << endl;
    cerr << endl;
    cerr << "Options:" << endl;
    cerr << "  --dump          Print table contents to stderr after building" << endl;
    cerr << "  --debug         Enable verbose progress output" << endl;
    cerr << "  --only-table1   Build only Table 1 (skip Table 2)" << endl;
    cerr << "  --threads N     Number of threads (default: OMP_NUM_THREADS or max)" << endl;
    cerr << "  --help          Show this help" << endl;
}

/// Build the "sample#phase#contig" base name from GBWT metadata and a PathName.
static string build_base_name(const gbwt::Metadata& meta, const gbwt::PathName& pn) {
    string sample = (meta.samples() > 0 && pn.sample < meta.samples())
                    ? meta.sample(pn.sample) : to_string(pn.sample);
    string contig = (meta.contigs() > 0 && pn.contig < meta.contigs())
                    ? meta.contig(pn.contig) : to_string(pn.contig);
    return sample + "#" + to_string(pn.phase) + "#" + contig;
}

/// Build the "sample#phase" haplotype name (no contig) from a PathName.
static string build_haplotype_name(const gbwt::Metadata& meta, const gbwt::PathName& pn) {
    string sample = (meta.samples() > 0 && pn.sample < meta.samples())
                    ? meta.sample(pn.sample) : to_string(pn.sample);
    return sample + "#" + to_string(pn.phase);
}

// ============================================================
// Path metadata (lightweight; no node sequence)
// ============================================================

struct PathMetadata {
    size_t path_id        = 0;
    string base_name;
    string haplotype_name;
    string contig_name;
    size_t subpath_start  = 0;
    size_t length        = 0;  ///< length in bases (0 = skip)
};

// Full path data including node sequence; used only when loading a path for Phase 2.
struct PathData {
    size_t              path_id        = 0;
    string              base_name;
    string              haplotype_name;
    string              contig_name;
    size_t              subpath_start  = 0;
    size_t              length         = 0;
    gbwt::vector_type   nodes;
    vector<size_t>      node_offsets;
};

/// Load one path's node sequence and offsets from GBWT/graph (used in Phase 2 per-contig).
static PathData load_path_data(size_t path_id,
                               const gbwt::GBWT& gbwt_index,
                               const gbwtgraph::GBWTGraph& graph,
                               const gbwt::Metadata& meta,
                               const PathMetadata& pm)
{
    PathData pd;
    pd.path_id        = path_id;
    pd.base_name      = pm.base_name;
    pd.haplotype_name = pm.haplotype_name;
    pd.contig_name    = pm.contig_name;
    pd.subpath_start  = pm.subpath_start;
    pd.length         = pm.length;

    pd.nodes = gbwt_index.extract(gbwt::Path::encode(path_id, false));
    pd.node_offsets.reserve(pd.nodes.size() + 1);
    pd.node_offsets.push_back(0);
    size_t total = 0;
    for (gbwt::node_type node : pd.nodes) {
        if (node == gbwt::ENDMARKER) break;
        total += graph.get_length(
            graph.get_handle(gbwt::Node::id(node), gbwt::Node::is_reverse(node)));
        pd.node_offsets.push_back(total);
    }
    return pd;
}

// ============================================================
// Phase 1: Build Table 1 and collect path metadata only (no node sequences)
// ============================================================

static vector<PathMetadata> build_table1_and_collect_metadata(
        const gbwt::GBWT&        gbwt_index,
        const gbwtgraph::GBWTGraph& graph,
        const gbwt::Metadata&    meta,
        panindexer::TranslationTable1& table1,
        bool debug)
{
    size_t num_paths = meta.paths();
    vector<PathMetadata> all_meta(num_paths);

    // Parallel: extract each path only to compute length; do not store nodes.
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t path_id = 0; path_id < num_paths; ++path_id) {
        gbwt::PathName pn = meta.path(path_id);

        PathMetadata pm;
        pm.path_id        = path_id;
        pm.base_name      = build_base_name(meta, pn);
        pm.haplotype_name = build_haplotype_name(meta, pn);
        pm.contig_name    = (meta.contigs() > 0 && pn.contig < meta.contigs())
                            ? meta.contig(pn.contig) : to_string(pn.contig);
        pm.subpath_start  = static_cast<size_t>(pn.count);

        gbwt::vector_type nodes = gbwt_index.extract(gbwt::Path::encode(path_id, false));
        size_t total = 0;
        for (gbwt::node_type node : nodes) {
            if (node == gbwt::ENDMARKER) break;
            total += graph.get_length(
                graph.get_handle(gbwt::Node::id(node), gbwt::Node::is_reverse(node)));
        }
        pm.length = total;
        all_meta[path_id] = std::move(pm);
    }

    // Sequential: add to Table 1 in path_id order.
    size_t num_added = 0, num_skipped = 0;
    for (size_t path_id = 0; path_id < num_paths; ++path_id) {
        PathMetadata& pm = all_meta[path_id];
        if (pm.length == 0) {
            if (debug) {
                cerr << "  [skip] path_id=" << path_id
                     << " (" << pm.base_name << "[" << pm.subpath_start << "]) length=0" << endl;
            }
            ++num_skipped;
            continue;
        }
        table1.add_subpath(pm.base_name, path_id, pm.subpath_start, pm.length);
        ++num_added;
        if (debug) {
            cerr << "  [T1] path_id=" << path_id << " \"" << pm.base_name << "\""
                 << " offset=" << pm.subpath_start << " len=" << pm.length << endl;
        } else if (path_id % 500 == 0) {
            cerr << "  Phase1: " << path_id << "/" << num_paths << "\r" << flush;
        }
    }
    cerr << endl;
    cerr << "Table 1: " << num_added << " subpaths added, "
         << num_skipped << " skipped, "
         << table1.num_names() << " named paths." << endl;

    return all_meta;
}

// ============================================================
// Phase 2: Build Table 2 using PathData collected in Phase 1
// ============================================================

/**
 * For a given (src_path, tgt_path) pair, find all contiguous matching segments.
 *
 * Strategy:
 *   Walk the src_path node by node. For each node, check if tgt_path also
 *   visits the same graph node. If it does, and it appears at the same
 *   relative position (consecutive common-node run), extend the current
 *   segment; otherwise flush the current segment and start a new one.
 *
 * This is an O(src_len * tgt_len) approach in the worst case, but in practice
 * most paths share only a small fraction of nodes, and we exit early per node.
 *
 * We use the pre-extracted tgt node list: build a map from node_id → list of
 * (tgt_node_idx, local_base_offset) for the target, then scan the source.
 */
static void find_common_segments(
        const PathData& src,
        const PathData& tgt,
        const gbwtgraph::GBWTGraph& graph,
        panindexer::TranslationTable2& table2,
        const string& tgt_haplotype,
        bool debug)
{
    if (src.nodes.empty() || tgt.nodes.empty()) return;

    // Build node_id → list of (tgt_node_idx) for the target path.
    // We use node_id (orientation-independent) to handle potential strand flips,
    // but we require orientation to match for a valid colinear alignment.
    unordered_map<gbwt::node_type, vector<size_t>> tgt_node_positions;
    tgt_node_positions.reserve(tgt.nodes.size());
    for (size_t ti = 0; ti < tgt.nodes.size(); ++ti) {
        if (tgt.nodes[ti] == gbwt::ENDMARKER) break;
        tgt_node_positions[tgt.nodes[ti]].push_back(ti);
    }

    // Walk the source path and find matching target positions.
    // We track "active segments": positions in the target where we are
    // currently building a colinear segment with the source.
    // Key: tgt_idx of the *next expected* tgt node; Value: (src_start, tgt_start).
    // We only ever match each target node once, so active segments are disjoint.

    struct ActiveSeg {
        size_t src_start   = 0;
        size_t tgt_start   = 0;
        size_t tgt_idx_next = 0;  ///< index into tgt.nodes we expect next
    };
    // At most one active segment: guarantees no overlapping segments (same src can't map to two target runs).
    ActiveSeg active_seg;
    bool has_active = false;
    size_t last_tgt_end = 0;  // only start a new run at target position >= this (avoid reusing target range)

    auto flush_segment = [&](size_t src_end, size_t tgt_end, const ActiveSeg& seg) {
        if (src_end <= seg.src_start || tgt_end <= seg.tgt_start) return;
        panindexer::IntervalMapping im;
        im.src_start   = seg.src_start;
        im.src_end     = src_end;
        im.tgt_path_id = tgt.path_id;
        table2.add_mapping(src.path_id, tgt_haplotype, im);
        if (debug) {
            cerr << "    seg: src=[" << im.src_start << "," << im.src_end << ")"
                 << " tgt_pid=" << im.tgt_path_id << endl;
        }
    };

    for (size_t si = 0; si < src.nodes.size(); ++si) {
        gbwt::node_type src_node = src.nodes[si];
        if (src_node == gbwt::ENDMARKER) break;

        auto tgt_it = tgt_node_positions.find(src_node);
        size_t src_start_here = src.node_offsets[si];

        // 1) If we have an active segment and the next expected target node is this source node, extend it.
        if (has_active && active_seg.tgt_idx_next < tgt.nodes.size() &&
            tgt.nodes[active_seg.tgt_idx_next] == src_node) {
            active_seg.tgt_idx_next++;
            if (active_seg.tgt_idx_next >= tgt.nodes.size()) {
                size_t tgt_end = tgt.node_offsets[active_seg.tgt_idx_next];
                if (tgt_end > last_tgt_end) last_tgt_end = tgt_end;
                flush_segment(src.node_offsets[si + 1], tgt_end, active_seg);
                has_active = false;
            }
            continue;
        }

        // 2) Can't extend: flush current segment (if any) with correct length, then maybe start a new run.
        if (has_active) {
            size_t tgt_end = tgt.node_offsets[active_seg.tgt_idx_next];
            if (tgt_end > last_tgt_end) last_tgt_end = tgt_end;
            size_t matched_len = tgt_end - active_seg.tgt_start;
            size_t src_end = active_seg.src_start + matched_len;
            flush_segment(src_end, tgt_end, active_seg);
            has_active = false;
        }

        // 3) Start a new segment at this position if this node appears in the target at or after last_tgt_end.
        if (tgt_it != tgt_node_positions.end()) {
            size_t ti_min = static_cast<size_t>(-1);  // pick smallest ti with tgt.node_offsets[ti] >= last_tgt_end
            for (size_t ti : tgt_it->second) {
                if (tgt.node_offsets[ti] >= last_tgt_end && (ti_min == static_cast<size_t>(-1) || ti < ti_min)) {
                    ti_min = ti;
                }
            }
            if (ti_min == static_cast<size_t>(-1)) continue;  // no valid target position
            active_seg.src_start    = src_start_here;
            active_seg.tgt_start    = tgt.node_offsets[ti_min];
            active_seg.tgt_idx_next = ti_min + 1;
            has_active = true;
            if (active_seg.tgt_idx_next >= tgt.nodes.size()) {
                size_t tgt_end = tgt.node_offsets[active_seg.tgt_idx_next];
                if (tgt_end > last_tgt_end) last_tgt_end = tgt_end;
                flush_segment(src.node_offsets[si + 1], tgt_end, active_seg);
                has_active = false;
            }
        }
    }

    if (has_active) {
        size_t tgt_end = tgt.node_offsets[active_seg.tgt_idx_next];
        size_t matched_len = tgt_end - active_seg.tgt_start;
        size_t src_end = active_seg.src_start + matched_len;
        flush_segment(src_end, tgt_end, active_seg);
    }
}

/// Same as find_common_segments but appends (src_path_id, tgt_haplotype, IntervalMapping) to out (for parallel build).
using Table2Record = std::tuple<size_t, std::string, panindexer::IntervalMapping>;
static void find_common_segments_to_vector(
        const PathData& src,
        const PathData& tgt,
        const gbwtgraph::GBWTGraph& graph,
        const string& tgt_haplotype,
        vector<Table2Record>* out)
{
    if (src.nodes.empty() || tgt.nodes.empty() || !out) return;

    unordered_map<gbwt::node_type, vector<size_t>> tgt_node_positions;
    tgt_node_positions.reserve(tgt.nodes.size());
    for (size_t ti = 0; ti < tgt.nodes.size(); ++ti) {
        if (tgt.nodes[ti] == gbwt::ENDMARKER) break;
        tgt_node_positions[tgt.nodes[ti]].push_back(ti);
    }

    struct ActiveSeg {
        size_t src_start   = 0;
        size_t tgt_start   = 0;
        size_t tgt_idx_next = 0;
    };
    ActiveSeg active_seg;
    bool has_active = false;
    size_t last_tgt_end = 0;

    auto flush_segment = [&](size_t src_end, size_t tgt_end, const ActiveSeg& seg) {
        if (src_end <= seg.src_start || tgt_end <= seg.tgt_start) return;
        panindexer::IntervalMapping im;
        im.src_start   = seg.src_start;
        im.src_end     = src_end;
        im.tgt_path_id = tgt.path_id;
        out->emplace_back(src.path_id, tgt_haplotype, im);
    };

    for (size_t si = 0; si < src.nodes.size(); ++si) {
        gbwt::node_type src_node = src.nodes[si];
        if (src_node == gbwt::ENDMARKER) break;

        auto tgt_it = tgt_node_positions.find(src_node);
        size_t src_start_here = src.node_offsets[si];

        if (has_active && active_seg.tgt_idx_next < tgt.nodes.size() &&
            tgt.nodes[active_seg.tgt_idx_next] == src_node) {
            active_seg.tgt_idx_next++;
            if (active_seg.tgt_idx_next >= tgt.nodes.size()) {
                size_t tgt_end = tgt.node_offsets[active_seg.tgt_idx_next];
                if (tgt_end > last_tgt_end) last_tgt_end = tgt_end;
                flush_segment(src.node_offsets[si + 1], tgt_end, active_seg);
                has_active = false;
            }
            continue;
        }

        if (has_active) {
            size_t tgt_end = tgt.node_offsets[active_seg.tgt_idx_next];
            if (tgt_end > last_tgt_end) last_tgt_end = tgt_end;
            size_t matched_len = tgt_end - active_seg.tgt_start;
            size_t src_end = active_seg.src_start + matched_len;
            flush_segment(src_end, tgt_end, active_seg);
            has_active = false;
        }

        if (tgt_it != tgt_node_positions.end()) {
            size_t ti_min = static_cast<size_t>(-1);
            for (size_t ti : tgt_it->second) {
                if (tgt.node_offsets[ti] >= last_tgt_end && (ti_min == static_cast<size_t>(-1) || ti < ti_min)) {
                    ti_min = ti;
                }
            }
            if (ti_min == static_cast<size_t>(-1)) continue;
            active_seg.src_start    = src_start_here;
            active_seg.tgt_start    = tgt.node_offsets[ti_min];
            active_seg.tgt_idx_next = ti_min + 1;
            has_active = true;
            if (active_seg.tgt_idx_next >= tgt.nodes.size()) {
                size_t tgt_end = tgt.node_offsets[active_seg.tgt_idx_next];
                if (tgt_end > last_tgt_end) last_tgt_end = tgt_end;
                flush_segment(src.node_offsets[si + 1], tgt_end, active_seg);
                has_active = false;
            }
        }
    }

    if (has_active) {
        size_t tgt_end = tgt.node_offsets[active_seg.tgt_idx_next];
        size_t matched_len = tgt_end - active_seg.tgt_start;
        size_t src_end = active_seg.src_start + matched_len;
        flush_segment(src_end, tgt_end, active_seg);
    }
}

static void build_table2(
        const gbwt::GBWT&          gbwt_index,
        const gbwtgraph::GBWTGraph& graph,
        const gbwt::Metadata&      meta,
        const vector<PathMetadata>& path_meta,
        panindexer::TranslationTable2& table2,
        bool debug)
{
    // Group path_ids by contig (metadata only; no node data yet).
    unordered_map<string, vector<size_t>> contig_to_path_ids;
    for (size_t path_id = 0; path_id < path_meta.size(); ++path_id) {
        const PathMetadata& pm = path_meta[path_id];
        if (pm.length == 0) continue;
        contig_to_path_ids[pm.contig_name].push_back(path_id);
    }

    size_t total_pairs = 0;
    for (const auto& kv : contig_to_path_ids) {
        size_t n = kv.second.size();
        if (n < 2) continue;
        total_pairs += n * (n - 1);
    }
    cerr << "Table 2: " << contig_to_path_ids.size() << " contigs, "
         << total_pairs << " (src, tgt) pairs to process." << endl;

    size_t contig_idx = 0;
    size_t num_contigs = contig_to_path_ids.size();

    for (const auto& kv : contig_to_path_ids) {
        const vector<size_t>& path_ids = kv.second;
        if (path_ids.size() < 2) continue;

        // Load only this contig's paths (peak memory = one contig's paths, not all).
        vector<PathData> contig_paths;
        contig_paths.reserve(path_ids.size());
        for (size_t path_id : path_ids) {
            contig_paths.push_back(
                load_path_data(path_id, gbwt_index, graph, meta, path_meta[path_id]));
        }

        // Build pairs (si, ti) = indices into contig_paths; skip same haplotype.
        vector<std::pair<size_t, size_t>> pairs;
        for (size_t si = 0; si < contig_paths.size(); ++si) {
            for (size_t ti = 0; ti < contig_paths.size(); ++ti) {
                if (si == ti) continue;
                if (contig_paths[si].haplotype_name == contig_paths[ti].haplotype_name) continue;
                pairs.emplace_back(si, ti);
            }
        }

        #pragma omp parallel
        {
            vector<Table2Record> local;
            local.reserve(256);
            #pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < pairs.size(); ++i) {
                size_t si = pairs[i].first, ti = pairs[i].second;
                find_common_segments_to_vector(
                    contig_paths[si], contig_paths[ti], graph,
                    contig_paths[ti].haplotype_name, &local);
            }
            #pragma omp critical
            {
                for (const auto& rec : local) {
                    table2.add_mapping(std::get<0>(rec), std::get<1>(rec), std::get<2>(rec));
                }
            }
        }

        // Discard contig_paths before next contig (free memory).
        contig_paths.clear();
        contig_paths.shrink_to_fit();

        ++contig_idx;
        if (contig_idx % 100 == 0 || contig_idx == num_contigs) {
            cerr << "  Phase2: contig " << contig_idx << "/" << num_contigs << "\r" << flush;
        }
    }
    cerr << endl;

    table2.finalize();
    cerr << "Table 2: " << table2.num_entries() << " (src_path, tgt_haplotype) entries, "
         << table2.total_segments() << " total segments." << endl;
}

// ============================================================
// main
// ============================================================

int main(int argc, char** argv) {
    if (argc < 4) {
        usage(argv[0]);
        return 1;
    }

    string gbz_file    = argv[1];
    string output_t1   = argv[2];
    string output_t2   = argv[3];
    bool dump          = false;
    bool debug         = false;
    bool only_table1   = false;
    int num_threads    = 0;  // 0 = use default (OMP_NUM_THREADS / omp_get_max_threads())

    for (int i = 4; i < argc; ++i) {
        string arg = argv[i];
        if      (arg == "--dump")        dump = true;
        else if (arg == "--debug")       debug = true;
        else if (arg == "--only-table1") only_table1 = true;
        else if (arg == "--threads") {
            if (i + 1 >= argc) { cerr << "Error: --threads requires N" << endl; return 1; }
            num_threads = atoi(argv[++i]);
            if (num_threads < 1) num_threads = 1;
        }
        else if (arg == "--help" || arg == "-h") { usage(argv[0]); return 0; }
        else { cerr << "Unknown argument: " << arg << endl; usage(argv[0]); return 1; }
    }

    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
    }
    cerr << "Using " << (num_threads > 0 ? num_threads : omp_get_max_threads()) << " thread(s)." << endl;

    auto t_start = high_resolution_clock::now();

    // ---- Load GBZ ----
    cerr << "Loading GBZ: " << gbz_file << " ..." << endl;
    gbwtgraph::GBZ gbz;
    try {
        sdsl::simple_sds::load_from(gbz, gbz_file);
    } catch (const exception& e) {
        cerr << "Error loading GBZ: " << e.what() << endl;
        return 1;
    }
    const gbwt::GBWT&          gbwt_index = gbz.index;
    const gbwtgraph::GBWTGraph& graph      = gbz.graph;

    auto t_loaded = high_resolution_clock::now();
    cerr << "GBZ loaded in "
         << duration_cast<milliseconds>(t_loaded - t_start).count()
         << " ms. Paths: " << gbwt_index.sequences() / 2
         << ", nodes: " << graph.get_node_count() << endl;

    if (!gbwt_index.hasMetadata()) {
        cerr << "Error: GBWT index has no metadata." << endl;
        return 1;
    }
    const gbwt::Metadata& meta = gbwt_index.metadata;
    cerr << "Metadata: " << meta.paths() << " paths, "
         << meta.samples() << " samples, "
         << meta.contigs() << " contigs" << endl;

    // ---- Phase 1: Table 1 (metadata only; no path node sequences stored) ----
    cerr << "\n=== Phase 1: Building Translation Table 1 ===" << endl;
    panindexer::TranslationTable1 table1;
    auto t_p1_start = high_resolution_clock::now();

    vector<PathMetadata> path_meta = build_table1_and_collect_metadata(
        gbwt_index, graph, meta, table1, debug);

    auto t_p1_end = high_resolution_clock::now();
    cerr << "Phase 1 done in "
         << duration_cast<milliseconds>(t_p1_end - t_p1_start).count() << " ms." << endl;

    // ---- Phase 2: Table 2 (stream by contig; load path nodes only for current contig) ----
    panindexer::TranslationTable2 table2;
    if (!only_table1) {
        cerr << "\n=== Phase 2: Building Translation Table 2 ===" << endl;
        auto t_p2_start = high_resolution_clock::now();

        build_table2(gbwt_index, graph, meta, path_meta, table2, debug);

        auto t_p2_end = high_resolution_clock::now();
        cerr << "Phase 2 done in "
             << duration_cast<milliseconds>(t_p2_end - t_p2_start).count() << " ms." << endl;
    }

    // ---- Dump (optional) ----
    if (dump) {
        cerr << "\n=== Table 1 contents ===" << endl;
        for (const string& name : table1.names()) {
            auto subpaths = table1.subpaths(name);
            cerr << name << "  (" << subpaths.size() << " subpath(s))" << endl;
            for (const auto& sp : subpaths) {
                cerr << "  path_id=" << sp.path_id
                     << "  offset=" << sp.subpath_start
                     << "  len=" << sp.length
                     << "  global=[" << sp.subpath_start << "," << sp.end() << ")" << endl;
            }
        }

        if (!only_table1) {
            cerr << "\n=== Table 2 contents ===" << endl;
            for (const auto& [src_pid, tgt_hap] : table2.keys()) {
                auto segs = table2.segments(src_pid, tgt_hap);
                cerr << "src_path=" << src_pid << " tgt_hap=\"" << tgt_hap
                     << "\"  (" << segs.size() << " seg(s))" << endl;
                for (const auto& seg : segs) {
                    cerr << "  src=[" << seg.src_start << "," << seg.src_end << ")"
                         << " tgt_pid=" << seg.tgt_path_id << endl;
                }
            }
        }
    }

    // ---- Serialize Table 1 ----
    cerr << "\nWriting Table 1 to: " << output_t1 << " ..." << endl;
    {
        ofstream out(output_t1, ios::binary);
        if (!out) { cerr << "Error: cannot open " << output_t1 << endl; return 1; }
        table1.serialize(out);
    }

    // ---- Serialize Table 2 ----
    if (!only_table1) {
        cerr << "Writing Table 2 to: " << output_t2 << " ..." << endl;
        {
            ofstream out(output_t2, ios::binary);
            if (!out) { cerr << "Error: cannot open " << output_t2 << endl; return 1; }
            table2.serialize(out);
        }
    }

    auto t_end = high_resolution_clock::now();
    cerr << "\nDone. Total time: "
         << duration_cast<milliseconds>(t_end - t_start).count() << " ms." << endl;

    return 0;
}
