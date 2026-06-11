#include "pangenome_index/algorithm.hpp"
#include "pangenome_index/tag_arrays.hpp"
#include "pangenome_index/light_tag_index.hpp"
#include <chrono>
#include <unordered_set>
#include <sys/resource.h>
#include <unistd.h>

using namespace std;
using namespace gbwtgraph;
using namespace panindexer;

#ifndef TIME
#define TIME 1
#endif

#ifndef DEBUG
#define DEBUG 0
#endif

const std::string SEPARATOR = "\t";

// ============================================================================
// v2 record format (find_mems → gafpack contract). See:
//   mem-projection/pangenome-pipeline/PLAN_find_mems_binary_io_v2.md
//
// node_id and offset are NO LONGER stored on disk. gafpack derives them by
// walking the path's cum_bp prefix-sum and locating the step that contains
// path_bp. Records are sorted by path_bp within each seq_id bucket so the
// path-walker can advance a single monotonic cursor (linear merge, O(1)
// amortized per record).
//
// NOTE: is_rev(graph_pos) is NOT stored. It is the strand of the BWT hit
// against the underlying *node* sequence, which is independent of the path's
// bucket parity (seq_id & 1 = forward/reverse path strand). gafpack derives
// the GAF strand from the bucket parity XOR the GFA step's +/- orientation
// (see traverse_nodes, main.rs). Earlier v2 drafts stored is_rev as a sanity
// bit and asserted equality with bucket parity; this was incorrect and
// silently dropped ~1.4M valid MEM hits on yeast-235.
// ============================================================================

// In-memory MEM hit; accumulated during read processing, bucket-sorted at end.
// seq_id is kept here for bucketing; stripped before write.
struct PackedEntry {
    uint32_t seq_id;
    uint32_t path_bp;   // bp offset of hit within seq_id's text (= r_index.seqOffset(sa_value))
    uint32_t match_len;
    uint32_t read_st;
    uint32_t read_id;
};

// On-disk record for _path_pos_v2.bin (consumed by gafpack). Little-endian,
// 16 bytes, naturally 8-byte aligned (no padding).
struct BinRecordV2 {
    uint32_t path_bp;
    uint32_t match_len;
    uint32_t read_st;
    uint32_t read_id;
};
static_assert(sizeof(BinRecordV2) == 16, "BinRecordV2 must be 16 bytes");

// Helper function to get current memory usage in MB
double get_memory_usage_mb() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
#ifdef __linux__
    return usage.ru_maxrss / 1024.0; // Linux reports in KB
#else
    return usage.ru_maxrss / (1024.0 * 1024.0); // macOS reports in bytes
#endif
}

std::string format_number(size_t n) {
    if (n >= 1000000000) {
        return std::to_string(n / 1000000000) + "." + std::to_string((n % 1000000000) / 100000000) + "B";
    } else if (n >= 1000000) {
        return std::to_string(n / 1000000) + "." + std::to_string((n % 1000000) / 100000) + "M";
    } else if (n >= 1000) {
        return std::to_string(n / 1000) + "." + std::to_string((n % 1000) / 100) + "K";
    }
    return std::to_string(n);
}

// Profiling structure for per-MEM statistics
struct MEMProfilingStats {
    double tag_query_time = 0.0;
    double locate_time = 0.0;
    double first_locate_time = 0.0;    // Time for initial locate_sa_value call
    double locate_next_time = 0.0;     // Aggregated time for all locateNext calls
    double file_write_time = 0.0;
    size_t tag_runs = 0;
    size_t locate_operations = 0;
    size_t first_locate_calls = 0;     // Number of locate_sa_value calls (1 per MEM)
    size_t locate_next_calls = 0;      // Number of locateNext calls
    size_t entries_written = 0;
};

// Global profiling structure
struct ProfilingData {
    // R-index
    double rindex_load_time = 0.0;
    double rindex_load_memory_mb = 0.0;
    
    // Tag index
    double tag_index_load_time = 0.0;
    double tag_index_load_memory_mb = 0.0;
    
    // Reads
    size_t total_reads = 0;
    double total_reads_processing_time = 0.0;
    size_t total_mems_outputted = 0;
    double total_mem_finding_time = 0.0;
    double total_mem_processing_time = 0.0;  // Time spent processing all MEMs (tag queries + locate + write)
    
    // MEM statistics
    size_t total_mem_length = 0;      // Sum of all MEM lengths
    size_t total_mem_occurrences = 0;  // Sum of all MEM occurrences (mem.size)
    double total_n_over_r_ratio = 0.0; // Sum of n/r ratios (n=mem.size, r=decoded_runs.size)
    
    // Per-MEM operations (aggregated)
    double total_tag_query_time = 0.0;
    double total_locate_time = 0.0;
    double total_first_locate_time = 0.0;   // Aggregated time for initial locate_sa_value calls
    double total_locate_next_time = 0.0;    // Aggregated time for all locateNext calls
    double total_file_write_time = 0.0;
    size_t total_tag_runs = 0;
    size_t total_locate_operations = 0;
    size_t total_first_locate_calls = 0;    // Total number of locate_sa_value calls
    size_t total_locate_next_calls = 0;     // Total number of locateNext calls
    size_t total_entries_written = 0;
    
    // Sorting
    double total_sort_time = 0.0;
    
    // Final memory
    double peak_memory_mb = 0.0;
};


// TODO Optimizations:
//  1. remove SeqID counter - instead compute running length when processing final output
//  2. Optimize for large no. of reads / streaming reads scenario - sorting n chunks and merging n sorted chunks
//  3. total entries at the start of the file to efficiently reserve space
//  4. Encode all entries in bytes instead of csv
//  5. Filtering: Order runs by length, this will far spurious matches are not picked

// // Helper function to locate sequence ID for a single BWT position
// // This is optimized for single-position queries to avoid vector allocation overhead
// // Returns the sequence ID for the given BWT position
// inline size_type locate_single_position(const FastLocate& r_index, size_type bwt_pos) {
//     size_type first = r_index.NO_POSITION;
//     size_type offset_of_first = bwt_pos;
//
//     // Find the nearest run start and get the sample
//     auto iter = r_index.blocks_start_pos.predecessor(bwt_pos);
//     size_type run_id = 0;
//     size_t cur_pos = 0;
//
//     auto run_num = r_index.blocks[iter->first].run_id_at(bwt_pos - iter->second, cur_pos);
//     run_id = iter->first * r_index.block_size + run_num;
//     offset_of_first = iter->second + cur_pos;
//     first = r_index.getSample(run_id);
//
//     // Iterate until the exact position
//     while (offset_of_first < bwt_pos) {
//         first = r_index.locateNext(first);
//         offset_of_first++;
//     }
//
//     // Return just the sequence ID
//     return r_index.seqId(first);
// }

// Helper function to get the SA value (not seq_id) for a single BWT position
// Returns the packed SA value that can be used with locateNext()
inline size_type locate_sa_value(const FastLocate& r_index, size_type bwt_pos) {
    size_t run_id = 0;
    size_t offset_of_first = 0;
    r_index.run_id_and_offset_at(bwt_pos, run_id, offset_of_first);
    size_type first = r_index.getSample(run_id);

    // Iterate until the exact position
    while (offset_of_first < bwt_pos) {
        first = r_index.locateNext(first);
        offset_of_first++;
    }

    return first;
}

// NOTE (v2): this "Option 1" variant is NOT updated for the v2 record format.
// It is currently unreferenced (see main() — the call site at line ~744 uses
// dump_mem_info_unique_runs instead). If reactivated, it must be ported to
// push PackedEntry (v2 layout) instead of writing the legacy text format.
//
// Helper function to dump MEM information (Original Algorithm)
// Locates SA entries for ALL BWT positions in the MEM interval
void dump_mem_info(const MEM& mem, const int read_id, TagArray& tag_array, FastLocate& r_index,
    vector<size_t> &seq_id_counter, MEMProfilingStats& mem_stats,
    std::ofstream* output_file = nullptr, bool debug_stats = false) {

    // Get BWT range for this MEM
    size_t bwt_start = mem.bwt_start;
    size_t bwt_end = mem.bwt_start + mem.size - 1;
    size_t mem_length = mem.end - mem.start;

#if TIME
    auto tag_query_start = chrono::high_resolution_clock::now();
#endif

    // Query tag array to get the runs for this BWT interval
    std::vector<std::pair<pos_t, uint16_t>> decoded_runs;
    tag_array.query_compressed_decoded_runs(mem.bwt_start, mem.bwt_start + mem.size - 1, decoded_runs);

#if TIME
    auto tag_query_end = chrono::high_resolution_clock::now();
    std::chrono::duration<double> tag_query_duration = tag_query_end - tag_query_start;
    mem_stats.tag_query_time = tag_query_duration.count();
    mem_stats.tag_runs = decoded_runs.size();
#endif

#if TIME
    auto locate_start = chrono::high_resolution_clock::now();
#endif

    // Get sequence IDs from r_index.locate()
    gbwt::range_type range(bwt_start, bwt_end);
    std::vector<size_type> sa_values = r_index.locate(range);
    
    // Track number of locate operations (one batch locate for the entire range)
    mem_stats.locate_operations = sa_values.size();

#if TIME
    auto locate_end = chrono::high_resolution_clock::now();
    std::chrono::duration<double> locate_duration = locate_end - locate_start;
    mem_stats.locate_time = locate_duration.count();
#endif

#if TIME
    auto write_start = chrono::high_resolution_clock::now();
#endif

    // Statistics tracking when debug is enabled
    std::unordered_map<size_type, size_t> seq_id_occurrence_count;
    size_t total_positions = 0;
    size_t unique_seq_ids = 0;
    size_t duplicate_occurrences = 0;

    // Output information for each BWT position in the interval
    size_t sa_index = 0;
    std::unordered_set<size_type> seen_seq_ids;
    for (size_t run_idx = 0; run_idx < decoded_runs.size() && sa_index < sa_values.size(); run_idx++) {
        pos_t graph_pos = decoded_runs[run_idx].first;
        uint16_t run_length = decoded_runs[run_idx].second;

        // Extract node_id, offset, strand from pos_t
        int64_t node_id = id(graph_pos);
        size_t offset = gbwtgraph::offset(graph_pos);
        bool is_reverse = is_rev(graph_pos);

        // Output for each BWT position covered by this run
        for (uint16_t pos_in_run = 0; pos_in_run < run_length && sa_index < sa_values.size(); pos_in_run++) {
            // seq_id from r_index.locate()
            size_type seq_id = sa_values[sa_index];

            // Statistics collection when debug is enabled
            if (debug_stats) {
                total_positions++;
                seq_id_occurrence_count[seq_id]++;

                if (seq_id_occurrence_count[seq_id] > 1) {
                    duplicate_occurrences++;
                }
            }

            // Only output if we haven't seen this seq_id before
            auto [it, inserted] = seen_seq_ids.insert(seq_id);
            if (inserted) {
                seq_id_counter[seq_id]++;
                unique_seq_ids++;
                mem_stats.entries_written++;

#if TIME
                auto file_write_start = chrono::high_resolution_clock::now();
#endif
                // Create output string
                std::string output_line = std::to_string(seq_id) + SEPARATOR +
                                        std::to_string(node_id) + SEPARATOR +
                                        std::to_string(offset) + SEPARATOR +
                                        std::to_string(is_reverse) + SEPARATOR +
                                        std::to_string(mem_length) + SEPARATOR +
                                        std::to_string(mem.start) + SEPARATOR +
                                        std::to_string(read_id);

                // Output to file if provided
                if (output_file && output_file->is_open()) {
                    *output_file << output_line << '\n';
                } else {
                    // Output to stdout
                    std::cout << output_line << std::endl;
                }

#if TIME
                auto file_write_end = chrono::high_resolution_clock::now();
                std::chrono::duration<double> file_write_duration = file_write_end - file_write_start;
                mem_stats.file_write_time += file_write_duration.count();
#endif
            }

            sa_index++;
        }
    }

    // Print statistics when debug is enabled
    if (debug_stats && (duplicate_occurrences > 0 || total_positions != unique_seq_ids)) {
        std::cerr << "=== MEM DUPLICATE STATISTICS ===" << std::endl;
        std::cerr << "Read ID: " << read_id << std::endl;
        std::cerr << "MEM: start=" << mem.start << ", end=" << mem.end
                  << ", length=" << (mem.end - mem.start) << std::endl;
        std::cerr << "BWT interval: [" << bwt_start << ", " << bwt_end
                  << "], size=" << mem.size << std::endl;
        std::cerr << "Total positions processed: " << total_positions << std::endl;
        std::cerr << "Unique sequence IDs: " << unique_seq_ids << std::endl;
        std::cerr << "Total duplicate occurrences: " << duplicate_occurrences << std::endl;

        std::cerr << "Sequences with duplicates:" << std::endl;
        for (const auto& pair : seq_id_occurrence_count) {
            if (pair.second > 1) {
                std::cerr << "  seq_id=" << pair.first << " appeared " << pair.second << " times" << std::endl;
            }
        }
        std::cerr << "=================================" << std::endl;
    }
}

// Helper function to dump MEM information (Optimized Algorithm)
// Only locates SA entries for the START of each unique tag run
// Uses amortized locate traversal but avoids allocating array for all positions
void dump_mem_info_unique_runs(const MEM& mem, const int read_id, TagArray& tag_array, FastLocate& r_index,
    vector<size_t> &seq_id_counter, std::vector<PackedEntry>& entries,
    MEMProfilingStats& mem_stats, bool debug_stats = false) {

    // Get BWT range for this MEM
    size_t bwt_start = mem.bwt_start;
    size_t bwt_end = mem.bwt_start + mem.size - 1;
    size_t mem_length = mem.end - mem.start;

#if TIME
    auto tag_query_start = chrono::high_resolution_clock::now();
#endif

    // Query tag array to get the runs for this BWT interval
    std::vector<std::pair<pos_t, uint16_t>> decoded_runs;
    tag_array.query_compressed_decoded_runs(mem.bwt_start, mem.bwt_start + mem.size - 1, decoded_runs);

#if TIME
    auto tag_query_end = chrono::high_resolution_clock::now();
    std::chrono::duration<double> tag_query_duration = tag_query_end - tag_query_start;
    mem_stats.tag_query_time = tag_query_duration.count();
    mem_stats.tag_runs = decoded_runs.size();
#endif

    if (decoded_runs.empty()) {
        return;
    }

#if TIME
    auto first_locate_start = chrono::high_resolution_clock::now();
#endif

    // Initialize locate traversal: get SA value at bwt_start
    size_type sa_value = locate_sa_value(r_index, bwt_start);
    mem_stats.first_locate_calls = 1;

#if TIME
    auto first_locate_end = chrono::high_resolution_clock::now();
    std::chrono::duration<double> first_locate_duration = first_locate_end - first_locate_start;
    mem_stats.first_locate_time = first_locate_duration.count();
#endif

    // Track unique graph positions to avoid duplicates
    std::unordered_set<pos_t> seen_graph_positions;
    
    // Statistics tracking when debug is enabled
    std::unordered_map<size_type, size_t> seq_id_occurrence_count;
    size_t total_runs_from_query = decoded_runs.size();
    size_t total_unique_positions = 0;
    size_t total_seq_ids_output = 0;
    size_t duplicate_positions_skipped = 0;
    
    // Track unique vs total node IDs from unique tag runs (for debug stats)
    std::unordered_set<int64_t> unique_node_ids;
    size_t total_node_ids_from_tag_runs = 0;
    
    // Track graph_pos -> total occurrence (sum of run lengths) when debug is enabled
    std::unordered_map<pos_t, size_t> graph_pos_occurrence_count;

    // Track node_id -> number of tag runs (in this MEM) whose tag starts at that node.
    // Counted once per physical tag-run-start (i.e. per pos_in_run == 0), regardless of
    // dedup. Lets us measure same-node repetition across non-adjacent runs within a MEM.
    std::unordered_map<int64_t, size_t> node_id_run_count;

    // Iterate through BWT interval, tracking position within current tag run
    size_t run_idx = 0;
    uint16_t pos_in_run = 0;
    
    for (size_t i = bwt_start; i <= bwt_end && run_idx < decoded_runs.size(); i++) {
        pos_t graph_pos = decoded_runs[run_idx].first;
        uint16_t run_length = decoded_runs[run_idx].second;

        // TODO: 05/01: calculate average run length
        
        // Process only at the START of each tag run (pos_in_run == 0)
        if (pos_in_run == 0) {
            if (debug_stats) {
                // std::cerr << "[DEBUG] run_idx=" << run_idx << " graph_pos=" << graph_pos << " run_length=" << run_length << std::endl;
                graph_pos_occurrence_count[graph_pos] += run_length;
                node_id_run_count[id(graph_pos)]++;
            }
            // Check if we've already seen this graph position
            auto [set_it, inserted] = seen_graph_positions.insert(graph_pos);
            if (!inserted) {
                // Duplicate graph position, skip it
                if (debug_stats) {
                    duplicate_positions_skipped++;
                }
            } else {
                total_unique_positions++;
                
                // Get seq_id and bp offset within that sequence for this position
                size_type seq_id = r_index.seqId(sa_value);
                size_type path_bp = r_index.seqOffset(sa_value);

                // Extract node_id, offset, strand from pos_t
                int64_t node_id = id(graph_pos);
                size_t offset = gbwtgraph::offset(graph_pos);
                bool is_reverse = is_rev(graph_pos);
                
                // Statistics collection when debug is enabled
                if (debug_stats) {
                    seq_id_occurrence_count[seq_id]++;
                    unique_node_ids.insert(node_id);
                    total_node_ids_from_tag_runs++;
                }

                // Output for this unique graph position
                seq_id_counter[seq_id]++;
                total_seq_ids_output++;
                mem_stats.entries_written++;

                // v2: only the fields gafpack actually needs. node_id/offset are
                // derived from path_bp by the walker; is_rev(graph_pos) is
                // independent of bucket parity and not used anywhere downstream.
                entries.push_back(PackedEntry{
                    static_cast<uint32_t>(seq_id),
                    static_cast<uint32_t>(path_bp),
                    static_cast<uint32_t>(mem_length),
                    static_cast<uint32_t>(mem.start),
                    static_cast<uint32_t>(read_id)
                });
            }

            // Early exit: we were at the start of the LAST run this iteration,
            // so we either just emitted or just skipped a duplicate at that
            // run's start. Either way, no more emissions are possible; the
            // remaining locateNext calls would advance the SA cursor through
            // the rest of the last run's BWT positions with no downstream
            // effect. Matches dump_mem_info_lightweight, which walks cur_bwt
            // only up to last_rid's start. Saves ~(last_run_length − 1)
            // locateNext calls per MEM.
            if (run_idx == decoded_runs.size() - 1) break;
        }

        // Advance position within current run
        pos_in_run++;
        
        // If we've reached the end of this run, move to the next run
        if (pos_in_run >= run_length) {
            run_idx++;
            pos_in_run = 0;
        }
        
#if TIME
        auto locate_next_start = chrono::high_resolution_clock::now();
#endif
        sa_value = r_index.locateNext(sa_value);
        mem_stats.locate_next_calls++;
#if TIME
        auto locate_next_end = chrono::high_resolution_clock::now();
        std::chrono::duration<double> locate_next_duration = locate_next_end - locate_next_start;
        mem_stats.locate_next_time += locate_next_duration.count();
#endif
    }
    
    // Compute total locate time and operations
    mem_stats.locate_time = mem_stats.first_locate_time + mem_stats.locate_next_time;
    mem_stats.locate_operations = mem_stats.first_locate_calls + mem_stats.locate_next_calls;

    // Print statistics when debug is enabled
    if (debug_stats) {
        std::cerr << "=== MEM STATISTICS (Unique Tag Runs Algorithm) ===" << std::endl;
        std::cerr << "Read ID: " << read_id << std::endl;
        std::cerr << "MEM: start=" << mem.start << ", end=" << mem.end
                  << ", size(no. of matches)=" << (mem.size)
        << ", length(len of mem)=" << (mem_length) << std::endl;
        std::cerr << "BWT interval: [" << bwt_start << ", " << bwt_end
                  << "], size=" << mem.size << std::endl;
        std::cerr << "Total tag runs from query: " << total_runs_from_query << std::endl;
        std::cerr << "Duplicate graph positions skipped: " << duplicate_positions_skipped << std::endl;
        std::cerr << "Unique graph positions processed: " << total_unique_positions << std::endl;
        std::cerr << "Total seq_ids output: " << total_seq_ids_output << std::endl;
        if (total_node_ids_from_tag_runs != unique_node_ids.size()) {
            std::cerr << "Total node IDs (from unique tag runs): " << total_node_ids_from_tag_runs << std::endl;
            std::cerr << "Unique node IDs: " << unique_node_ids.size() << std::endl;
        }

        // Count unique seq_ids
        std::cerr << "Unique sequence IDs: " << seq_id_occurrence_count.size() << std::endl;
        
        // Show seq_ids that appear multiple times (same seq_id with different graph positions)
        bool has_duplicates = false;
        for (const auto& pair : seq_id_occurrence_count) {
            if (pair.second > 1) {
                has_duplicates = true;
                break;
            }
        }
        
        if (has_duplicates) {
            std::cerr << "Sequences appearing in multiple unique graph positions:" << std::endl;
            for (const auto& pair : seq_id_occurrence_count) {
                if (pair.second > 1) {
                    std::cerr << "  seq_id=" << pair.first << " appears in " << pair.second << " unique positions" << std::endl;
                }
            }
        }
        
        // Output graph_pos -> total occurrence (sum of run lengths)
        std::cerr << "Graph position occurrence counts (sum of run lengths):" << std::endl;
        for (const auto& pair : graph_pos_occurrence_count) {
            std::cerr << "  graph_pos=" << pair.first << " total_occurrence=" << pair.second << std::endl;
        }

        // Per-MEM node repetition: how many tag runs start at each node_id.
        // run_count > 1 means the same node appears as the start of multiple
        // tag runs within this MEM (non-adjacent same-node tag-run repetition,
        // since adjacent same-position runs are already merged into one run).
        size_t nodes_with_repeats = 0;
        size_t total_node_runs = 0;
        size_t max_runs_per_node = 0;
        for (const auto& pr : node_id_run_count) {
            total_node_runs += pr.second;
            if (pr.second > 1) nodes_with_repeats++;
            if (pr.second > max_runs_per_node) max_runs_per_node = pr.second;
        }
        std::cerr << "Distinct node_ids: " << node_id_run_count.size()
                  << ", nodes with >1 run: " << nodes_with_repeats
                  << ", max runs/node: " << max_runs_per_node
                  << ", total node-run starts: " << total_node_runs << std::endl;
        std::cerr << "Node run-count histogram (runs_per_node count):" << std::endl;
        std::map<size_t, size_t> hist;
        for (const auto& pr : node_id_run_count) hist[pr.second]++;
        for (const auto& pr : hist) {
            std::cerr << "  " << pr.first << " run(s)/node: " << pr.second << " node(s)" << std::endl;
        }
        // Per-node detail for nodes appearing more than once (sorted descending).
        if (nodes_with_repeats > 0) {
            std::vector<std::pair<int64_t, size_t>> repeats;
            for (const auto& pr : node_id_run_count) if (pr.second > 1) repeats.push_back(pr);
            std::sort(repeats.begin(), repeats.end(),
                      [](const std::pair<int64_t,size_t>& a, const std::pair<int64_t,size_t>& b){
                          return a.second > b.second;
                      });
            std::cerr << "Node-run repeats (node_id: runs):" << std::endl;
            for (const auto& pr : repeats) {
                std::cerr << "  node_id=" << pr.first << " runs=" << pr.second << std::endl;
            }
        }
        std::cerr << "=================================" << std::endl;
    }
}

// Helper function to dump MEM information (Lightweight Tag Index Algorithm).
// Consumes a LightTagIndex which contains ONLY the bwt_intervals bitvector
// (no tag values). Emits one PackedEntry per tag run intersecting the MEM
// interval, using exactly one (seq_id, path_bp) per run. No graph-position
// dedup is performed (since tag values are unavailable); the per-MEM record
// count therefore exceeds dump_mem_info_unique_runs by the per-MEM duplicate
// fraction (~24% on yeast235 chrII; dataset-dependent).
//
// Algorithm:
//   1. Convert [bwt_start, bwt_end] -> [first_rid, last_rid] via rank.
//   2. Seed SA value at bwt_start with locate_sa_value.
//   3. For each run id in [first_rid+1 .. last_rid], walk locateNext from the
//      previous run start to the next run start.
//   4. At each run start, read seq_id/path_bp from the current SA value and
//      push a PackedEntry.
//
// Total locateNext calls per MEM = (last_run_start_bwt - bwt_start). The
// dump_mem_info_unique_runs path also stops at the last run start (after the
// matching early-exit fix); both paths now perform the same locate work.
// The remaining savings come from skipping graph-position decoding (no
// encoded_runs_iv reads) and the unordered_set<pos_t> allocation.
void dump_mem_info_lightweight(const MEM& mem, const int read_id,
                               const LightTagIndex& ltag, FastLocate& r_index,
                               vector<size_t>& seq_id_counter,
                               std::vector<PackedEntry>& entries,
                               MEMProfilingStats& mem_stats,
                               bool debug_stats = false) {
    const size_t bwt_start  = mem.bwt_start;
    const size_t bwt_end    = mem.bwt_start + mem.size - 1;
    const size_t mem_length = mem.end - mem.start;

#if TIME
    auto tag_query_start = chrono::high_resolution_clock::now();
#endif

    // Map [bwt_start, bwt_end] -> run id range via rank on bwt_intervals.
    const size_t first_rid = ltag.run_id_at(bwt_start);
    const size_t last_rid  = ltag.run_id_at(bwt_end);
    const size_t n_runs    = last_rid - first_rid + 1;
    mem_stats.tag_runs = n_runs;

#if TIME
    auto tag_query_end = chrono::high_resolution_clock::now();
    std::chrono::duration<double> tag_query_duration = tag_query_end - tag_query_start;
    mem_stats.tag_query_time = tag_query_duration.count();
#endif

    if (n_runs == 0) return;

    // Seed SA at the first run-start position covered by the MEM.
    // The first run may have started BEFORE bwt_start; we clip to bwt_start.
#if TIME
    auto first_locate_start = chrono::high_resolution_clock::now();
#endif
    size_type sa_value = locate_sa_value(r_index, bwt_start);
    mem_stats.first_locate_calls = 1;
#if TIME
    auto first_locate_end = chrono::high_resolution_clock::now();
    std::chrono::duration<double> first_locate_duration = first_locate_end - first_locate_start;
    mem_stats.first_locate_time = first_locate_duration.count();
#endif

    size_t cur_bwt = bwt_start;

    // Stats accumulators when debug_stats is enabled.
    size_t total_entries_in_mem = 0;
    std::unordered_map<size_type, size_t> seq_id_count;

    for (size_t rid = first_rid; rid <= last_rid; ++rid) {
        // BWT start of this run, clipped to bwt_start for the very first run
        // (since the first run may extend backward beyond the MEM).
        const size_t run_bwt = (rid == first_rid) ? bwt_start : ltag.run_start_bwt(rid);

        // Advance SA cursor from cur_bwt up to run_bwt by locateNext steps.
        while (cur_bwt < run_bwt) {
#if TIME
            auto locate_next_start = chrono::high_resolution_clock::now();
#endif
            sa_value = r_index.locateNext(sa_value);
            mem_stats.locate_next_calls++;
#if TIME
            auto locate_next_end = chrono::high_resolution_clock::now();
            std::chrono::duration<double> locate_next_duration = locate_next_end - locate_next_start;
            mem_stats.locate_next_time += locate_next_duration.count();
#endif
            ++cur_bwt;
        }

        // Emit one entry for this run-start position.
        const size_type seq_id  = r_index.seqId(sa_value);
        const size_type path_bp = r_index.seqOffset(sa_value);

        seq_id_counter[seq_id]++;
        mem_stats.entries_written++;
        ++total_entries_in_mem;
        if (debug_stats) seq_id_count[seq_id]++;

        entries.push_back(PackedEntry{
            static_cast<uint32_t>(seq_id),
            static_cast<uint32_t>(path_bp),
            static_cast<uint32_t>(mem_length),
            static_cast<uint32_t>(mem.start),
            static_cast<uint32_t>(read_id)
        });
    }

    mem_stats.locate_time = mem_stats.first_locate_time + mem_stats.locate_next_time;
    mem_stats.locate_operations = mem_stats.first_locate_calls + mem_stats.locate_next_calls;

    if (debug_stats) {
        std::cerr << "=== MEM STATISTICS (Lightweight Tag Index Algorithm) ===" << std::endl;
        std::cerr << "Read ID: " << read_id << std::endl;
        std::cerr << "MEM: start=" << mem.start << ", end=" << mem.end
                  << ", size(no. of matches)=" << mem.size
                  << ", length(len of mem)=" << mem_length << std::endl;
        std::cerr << "BWT interval: [" << bwt_start << ", " << bwt_end
                  << "], size=" << mem.size << std::endl;
        std::cerr << "Tag runs intersecting interval: " << n_runs << std::endl;
        std::cerr << "Entries emitted (no dedup): " << total_entries_in_mem << std::endl;
        std::cerr << "Distinct seq_ids in MEM: " << seq_id_count.size() << std::endl;
        size_t dup_seq = 0;
        for (const auto& p : seq_id_count) if (p.second > 1) ++dup_seq;
        std::cerr << "Seq_ids appearing in >1 entry: " << dup_seq << std::endl;
        std::cerr << "=================================" << std::endl;
    }
}

// Helper function to dump MEM information (All-Positions Algorithm).
// Verification / POC path: bypasses the tag array entirely and emits ONE
// PackedEntry per BWT position in [bwt_start, bwt_end] using only the
// r-index (locate_sa_value + locateNext). No run-start filtering, no
// graph-position dedup. seq_id and path_bp come straight from the SA value
// at each position.
//
// Output volume per MEM = mem.size (= total_mem_occurrences in profiling).
// At HPRC-noisy-alt scale this can be hundreds of millions of records;
// intended for small-scale A/B verification against the lightweight and
// full-tag dedup modes, not production use.
//
// MEMProfilingStats accounting:
//   tag_runs            = 0          (tag array not consulted)
//   first_locate_calls  = 1          (initial locate_sa_value)
//   locate_next_calls   = mem.size   (one per BWT position, including a final
//                                     unused step after the last emission;
//                                     trade a per-iteration branch for one
//                                     wasted locateNext per MEM)
//   entries_written     = mem.size
void dump_mem_info_all_positions(const MEM& mem, const int read_id, FastLocate& r_index,
                                 vector<size_t>& seq_id_counter,
                                 std::vector<PackedEntry>& entries,
                                 MEMProfilingStats& mem_stats,
                                 bool /*debug_stats*/ = false) {
    const size_t bwt_start  = mem.bwt_start;
    const size_t bwt_end    = mem.bwt_start + mem.size - 1;
    const size_t mem_length = mem.end - mem.start;

    // Tag array intentionally not consulted; tag_runs stays 0.
    mem_stats.tag_runs = 0;

    if (mem.size == 0) return;

    // Seed SA cursor at bwt_start.
#if TIME
    auto first_locate_start = chrono::high_resolution_clock::now();
#endif
    size_type sa_value = locate_sa_value(r_index, bwt_start);
    mem_stats.first_locate_calls = 1;
#if TIME
    auto first_locate_end = chrono::high_resolution_clock::now();
    std::chrono::duration<double> first_locate_duration = first_locate_end - first_locate_start;
    mem_stats.first_locate_time = first_locate_duration.count();
#endif

    // Walk every BWT position in the MEM interval; emit one entry each.
    for (size_t pos = bwt_start; pos <= bwt_end; ++pos) {
        const size_type seq_id  = r_index.seqId(sa_value);
        const size_type path_bp = r_index.seqOffset(sa_value);

        seq_id_counter[seq_id]++;
        mem_stats.entries_written++;

        entries.push_back(PackedEntry{
            static_cast<uint32_t>(seq_id),
            static_cast<uint32_t>(path_bp),
            static_cast<uint32_t>(mem_length),
            static_cast<uint32_t>(mem.start),
            static_cast<uint32_t>(read_id)
        });

        // Advance to the next BWT position. We unconditionally locateNext even
        // on the last iteration (result discarded next loop exit) -- trading
        // one wasted call per MEM for a tighter inner loop.
#if TIME
        auto locate_next_start = chrono::high_resolution_clock::now();
#endif
        sa_value = r_index.locateNext(sa_value);
        mem_stats.locate_next_calls++;
#if TIME
        auto locate_next_end = chrono::high_resolution_clock::now();
        std::chrono::duration<double> locate_next_duration = locate_next_end - locate_next_start;
        mem_stats.locate_next_time += locate_next_duration.count();
#endif
    }

    mem_stats.locate_time = mem_stats.first_locate_time + mem_stats.locate_next_time;
    mem_stats.locate_operations = mem_stats.first_locate_calls + mem_stats.locate_next_calls;
}

// Bucket-sort entries by seq_id (using prefix sums of seq_id_counter), then sort each
// bucket by path_bp ascending, and write _path_pos_v2.bin (+ optional .tsv) +
// _seq_id_starts.out. See PLAN_find_mems_binary_io_v2.md.
void write_sorted_entries(std::vector<PackedEntry>& entries, const std::string& output_prefix,
                          const vector<size_t>& seq_id_counter, bool emit_tsv) {
    const size_t num_seq = seq_id_counter.size();
    const size_t n = entries.size();
    std::cerr << "Sorting " << n << " entries across " << num_seq << " seq_id buckets" << std::endl;

    // Prefix sum -> bucket boundaries [0..num_seq]
    std::vector<size_t> bucket_start(num_seq + 1, 0);
    for (size_t s = 0; s < num_seq; s++) {
        bucket_start[s + 1] = bucket_start[s] + seq_id_counter[s];
    }
    assert(bucket_start[num_seq] == n);

    // Scatter into per-seq_id buckets (counting sort, O(n))
    std::vector<PackedEntry> out(n);
    {
        std::vector<size_t> cursor(bucket_start);
        for (const auto& e : entries) {
            out[cursor[e.seq_id]++] = e;
        }
    }
    // Reclaim input buffer before per-bucket sort / write
    entries.clear();
    entries.shrink_to_fit();

    // Per-bucket sort by path_bp ascending — v2 walker requires this for the
    // linear-merge cursor invariant (cum_bp[i] <= r.path_bp < cum_bp[i+1]).
    // Tie-break is unstable; if byte-stable .bin md5s are needed later, extend
    // the key to (path_bp, read_id, read_st). See PLAN_v2 §"Open questions".
    for (size_t s = 0; s < num_seq; s++) {
        auto first = out.begin() + bucket_start[s];
        auto last  = out.begin() + bucket_start[s + 1];
        std::sort(first, last, [](const PackedEntry& a, const PackedEntry& b) {
            return a.path_bp < b.path_bp;
        });
    }

    // Write _path_pos_v2.bin: contiguous 16-byte BinRecordV2 (seq_id stripped).
    // Single contiguous write per bucket avoids per-record ofstream overhead.
    {
        const std::string bin_path = output_prefix + "_path_pos_v2.bin";
        std::ofstream bin(bin_path, std::ios::binary);
        if (!bin.is_open()) {
            throw std::runtime_error("Cannot open output file: " + bin_path);
        }
        // Use a small reusable on-stack buffer; PackedEntry has an extra
        // seq_id field that must be stripped before write.
        for (const auto& e : out) {
            BinRecordV2 r{e.path_bp, e.match_len, e.read_st, e.read_id};
            bin.write(reinterpret_cast<const char*>(&r), sizeof(r));
        }
    }

    if (emit_tsv) {
        // v2 columns: path_bp match_len read_st read_id
        // (node_id/offset/is_rev dropped; recover node_id/offset via cum_bp walk if needed.)
        std::ofstream tsv(output_prefix + "_path_pos.tsv");
        if (!tsv.is_open()) {
            throw std::runtime_error("Cannot open output file: " + output_prefix + "_path_pos.tsv");
        }
        tsv << "# path_bp\tmatch_len\tread_st\tread_id\n";
        for (const auto& e : out) {
            tsv << e.path_bp << '\t' << e.match_len << '\t'
                << e.read_st << '\t' << e.read_id << '\n';
        }
    }

    // Write _seq_id_starts.out: bucket_start[0..=num_seq], one per line.
    // Unchanged from v1; values are record indices (byte offset = idx * 16 now).
    std::ofstream starts(output_prefix + "_seq_id_starts.out");
    if (!starts.is_open()) {
        throw std::runtime_error("Cannot open seq_id_starts file: " + output_prefix + "_seq_id_starts.out");
    }
    for (size_t s = 0; s <= num_seq; s++) {
        starts << bucket_start[s] << '\n';
    }
    starts.close();

    std::cerr << "Successfully wrote " << n << " sorted entries to " << output_prefix
              << "_path_pos_v2.bin (" << (n * sizeof(BinRecordV2)) << " bytes)" << std::endl;
}


int main(int argc, char **argv) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <r_index_file> <tag_array_index> <reads_file> <mem_length> <min_occ> [output_prefix] [--tsv] [--debug-stats] [--verbose] [--lightweight-tags|--all-positions]" << std::endl;
        std::cerr << "  --lightweight-tags: treat <tag_array_index> as a .ltags file (LightTagIndex)" << std::endl;
        std::cerr << "                      instead of a full compact tags file. Emits one entry per" << std::endl;
        std::cerr << "                      tag run intersecting each MEM, with no graph-position dedup." << std::endl;
        std::cerr << "  --all-positions:    verification/POC mode. Bypasses the tag array (the" << std::endl;
        std::cerr << "                      <tag_array_index> argument is ignored / can be /dev/null)." << std::endl;
        std::cerr << "                      Emits one entry per BWT position in each MEM interval" << std::endl;
        std::cerr << "                      (mem.size entries per MEM). Output volume can be very" << std::endl;
        std::cerr << "                      large; intended only for small-scale A/B against the" << std::endl;
        std::cerr << "                      lightweight/full-tag dedup modes." << std::endl;
        return EXIT_FAILURE;
    }

    string r_index_file = argv[1];
    string tag_array_index = argv[2];
    string reads_file = argv[3];
    size_t mem_length = std::stoi(argv[4]);
    size_t min_occ = std::stoi(argv[5]);
    string output_file_name;
    bool debug_stats = false;
    bool verbose = false;
    bool emit_tsv = false;
    bool lightweight_tags = false;
    bool all_positions = false;

    for (int i = 6; i < argc; i++) {
        std::string arg(argv[i]);
        if (arg == "--debug-stats") {
            debug_stats = true;
        } else if (arg == "--verbose" || arg == "--debug") {
            verbose = true;
        } else if (arg == "--tsv") {
            emit_tsv = true;
        } else if (arg == "--lightweight-tags") {
            lightweight_tags = true;
        } else if (arg == "--all-positions") {
            all_positions = true;
        } else if (output_file_name.empty()) {
            output_file_name = arg;
        }
    }
    if (lightweight_tags && all_positions) {
        std::cerr << "ERROR: --lightweight-tags and --all-positions are mutually exclusive" << std::endl;
        return EXIT_FAILURE;
    }

    // Print run parameters for logging
    std::cerr << "========================================" << std::endl;
    std::cerr << "find_mems - Run Parameters" << std::endl;
    std::cerr << "========================================" << std::endl;
    std::cerr << "R-index file:      " << r_index_file << std::endl;
    std::cerr << "Tag array index:   " << tag_array_index << std::endl;
    std::cerr << "Reads file:        " << reads_file << std::endl;
    std::cerr << "MEM length:        " << mem_length << std::endl;
    std::cerr << "Min occurrences:   " << min_occ << std::endl;
    std::cerr << "Output file:       " << (output_file_name.empty() ? "(stdout)" : output_file_name) << std::endl;
    std::cerr << "Debug stats:       " << (debug_stats ? "enabled" : "disabled") << std::endl;
    std::cerr << "Verbose mode:      " << (verbose ? "enabled" : "disabled") << std::endl;
    std::cerr << "Tag mode:          "
              << (all_positions   ? "all-positions (tag array bypassed)"
                  : lightweight_tags ? "lightweight (.ltags)"
                                     : "full compact tags") << std::endl;
    std::cerr << "========================================" << std::endl;

    // Initialize profiling data
    ProfilingData profiling;
    double initial_memory_mb = get_memory_usage_mb();

#if TIME
    auto rindex_load_start = chrono::high_resolution_clock::now();
#endif

    cerr << "Reading the rindex file" << endl;
    FastLocate r_index;
    {
        std::ifstream rin(r_index_file, std::ios::binary);
        if (!rin) {
            std::cerr << "Cannot open the r-index file " << r_index_file << std::endl;
            std::exit(EXIT_FAILURE);
        }
        r_index.load_encoded(rin);
    }

#if TIME
    auto rindex_load_end = chrono::high_resolution_clock::now();
    std::chrono::duration<double> rindex_load_duration = rindex_load_end - rindex_load_start;
    profiling.rindex_load_time = rindex_load_duration.count();
    profiling.rindex_load_memory_mb = get_memory_usage_mb() - initial_memory_mb;
    std::cerr << "Loading r-index into memory took " << profiling.rindex_load_time << " seconds" << std::endl;
    std::cerr << "R-index memory usage: " << profiling.rindex_load_memory_mb << " MB" << std::endl;
#endif

#if TIME
    auto tag_load_start = chrono::high_resolution_clock::now();
    double memory_before_tag = get_memory_usage_mb();
#endif

    // Load either the full compact TagArray or the LightTagIndex, depending on
    // --lightweight-tags. Exactly one of `tag_array` / `light_tag_index` is
    // populated; the other stays empty and is never accessed below. In
    // --all-positions mode neither is loaded (the tag_array_index CLI arg is
    // ignored).
    TagArray tag_array;
    LightTagIndex light_tag_index;
    if (all_positions) {
        cerr << "Skipping tag array load (--all-positions mode)" << endl;
    } else if (lightweight_tags) {
        cerr << "Reading the lightweight tag index (.ltags)" << endl;
        std::ifstream in_ds(tag_array_index, std::ios::binary);
        if (!in_ds) {
            std::cerr << "Cannot open the .ltags file " << tag_array_index << std::endl;
            std::exit(EXIT_FAILURE);
        }
        try {
            light_tag_index.load(in_ds);
        } catch (const std::exception& e) {
            std::cerr << "Failed to load .ltags file " << tag_array_index
                      << ": " << e.what() << std::endl;
            std::exit(EXIT_FAILURE);
        }
        std::cerr << "  bwt_size: " << light_tag_index.bwt_size()
                  << ", num_runs: " << light_tag_index.num_runs()
                  << ", memory: " << light_tag_index.memory_bytes() << " bytes"
                  << std::endl;
    } else {
        cerr << "Reading the tag array index" << endl;
        std::ifstream in_ds(tag_array_index);
        tag_array.load_compressed_tags_compact(in_ds);
    }

#if TIME
    auto tag_load_end = chrono::high_resolution_clock::now();
    std::chrono::duration<double> tag_load_duration = tag_load_end - tag_load_start;
    profiling.tag_index_load_time = tag_load_duration.count();
    profiling.tag_index_load_memory_mb = get_memory_usage_mb() - memory_before_tag;
    std::cerr << "Loading tag arrays took " << profiling.tag_index_load_time << " seconds" << std::endl;
    std::cerr << "Tag index memory usage: " << profiling.tag_index_load_memory_mb << " MB" << std::endl;
#endif

    // Open reads file
    std::ifstream reads(reads_file);
    if (!reads) {
        std::cerr << "Cannot open reads file: " << reads_file << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::string read;
    int i = 0;
    vector<size_t> seq_id_counter(r_index.tot_strings(), 0);
    std::cerr << "Reserved seq_id_ctr vector of length: " << r_index.tot_strings() << std::endl;
    std::vector<PackedEntry> entries;
    size_t total_mem_matches = 0;

#if TIME
    auto reads_processing_start = chrono::high_resolution_clock::now();
#endif

    while (std::getline(reads, read)) {
        if (read.empty()) continue;
        i++;
        profiling.total_reads++;

#if TIME
        auto mem_finding_start = chrono::high_resolution_clock::now();
#endif

        // Find MEMs for this read
        auto mems = find_all_mems(read, mem_length, min_occ, r_index);

#if TIME
        auto mem_finding_end = chrono::high_resolution_clock::now();
        std::chrono::duration<double> mem_finding_duration = mem_finding_end - mem_finding_start;
        profiling.total_mem_finding_time += mem_finding_duration.count();
#endif

#if DEBUG
        if (verbose) {
            std::cerr << "Read: " << i << std::endl;
            std::cerr << "Read sequence: " << read << std::endl;
        }
#endif
        
        profiling.total_mems_outputted += mems.size();
        
#if TIME
        auto mem_processing_start = chrono::high_resolution_clock::now();
#endif

        for (const auto& mem : mems) {
#if DEBUG
            if (verbose) {
                std::cerr << "MEM START: " << mem.start << ", MEM END: " << mem.end 
                          << " BWT START: " << mem.bwt_start << " BWT SIZE: " << mem.size << std::endl;
            }
#endif

            total_mem_matches += mem.size;
            
            // Track MEM statistics
            profiling.total_mem_length += (mem.end - mem.start);
            profiling.total_mem_occurrences += mem.size;

            // Per-MEM profiling
            MEMProfilingStats mem_stats;
            
            // Choose algorithm:
            // Option 1: Original - locates ALL BWT positions (slower, more SA lookups)
            // dump_mem_info(mem, i, tag_array, r_index, seq_id_counter, mem_stats,
            //     output_file.is_open() ? &output_file : nullptr, debug_stats);
            //
            // // Option 2 (default): Optimized - locates only START of each unique tag run,
            // dedup by graph_pos via unordered_set<pos_t>.
            //
            // Option 3 (--lightweight-tags): consumes a LightTagIndex (.ltags) with
            // bwt_intervals only. Emits one entry per intersecting tag run, no dedup.
            //
            // Option 4 (--all-positions): verification/POC. Skips the tag array
            // entirely; emits one entry per BWT position in [bwt_start, bwt_end].
            if (all_positions) {
                dump_mem_info_all_positions(mem, i, r_index,
                                            seq_id_counter, entries, mem_stats,
                                            debug_stats);
            } else if (lightweight_tags) {
                dump_mem_info_lightweight(mem, i, light_tag_index, r_index,
                                          seq_id_counter, entries, mem_stats,
                                          debug_stats);
            } else {
                dump_mem_info_unique_runs(mem, i, tag_array, r_index,
                                          seq_id_counter, entries, mem_stats,
                                          debug_stats);
            }

            // Aggregate per-MEM stats
            profiling.total_tag_query_time += mem_stats.tag_query_time;
            profiling.total_locate_time += mem_stats.locate_time;
            profiling.total_first_locate_time += mem_stats.first_locate_time;
            profiling.total_locate_next_time += mem_stats.locate_next_time;
            profiling.total_locate_operations += mem_stats.locate_operations;
            profiling.total_first_locate_calls += mem_stats.first_locate_calls;
            profiling.total_locate_next_calls += mem_stats.locate_next_calls;
            profiling.total_file_write_time += mem_stats.file_write_time;
            profiling.total_tag_runs += mem_stats.tag_runs;
            profiling.total_entries_written += mem_stats.entries_written;
            
            // Calculate n/r ratio (n = mem.size, r = number of decoded runs)
            if (mem_stats.tag_runs > 0) {
                double n_over_r = (double)mem.size / mem_stats.tag_runs;
                profiling.total_n_over_r_ratio += n_over_r;
            }
        }
        
#if TIME
        auto mem_processing_end = chrono::high_resolution_clock::now();
        std::chrono::duration<double> mem_processing_duration = mem_processing_end - mem_processing_start;
        profiling.total_mem_processing_time += mem_processing_duration.count();
#endif
        
#if DEBUG
        if (verbose && !mems.empty()) {
            std::cerr << std::endl;  // Add blank line after each read in verbose mode
        }
#endif
        
        // Update peak memory
        double current_memory = get_memory_usage_mb();
        if (current_memory > profiling.peak_memory_mb) {
            profiling.peak_memory_mb = current_memory;
        }
    }

#if TIME
    auto reads_processing_end = chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_reads_processing_duration = reads_processing_end - reads_processing_start;
    profiling.total_reads_processing_time = total_reads_processing_duration.count();
#endif

    reads.close();
    std::cerr << "Accumulated " << entries.size() << " entries in memory ("
              << (entries.size() * sizeof(PackedEntry)) / (1024.0 * 1024.0) << " MB)" << std::endl;

    if (!output_file_name.empty()) {
#if TIME
        auto sort_start = chrono::high_resolution_clock::now();
#endif

        write_sorted_entries(entries, output_file_name, seq_id_counter, emit_tsv);

#if TIME
        auto sort_end = chrono::high_resolution_clock::now();
        std::chrono::duration<double> sort_duration = sort_end - sort_start;
        profiling.total_sort_time = sort_duration.count();
#endif
    }

#if TIME
    std::cout << "\n================================================" << std::endl;
    std::cout << "=== COMPREHENSIVE PROFILING RESULTS ===" << std::endl;
    std::cout << "================================================" << std::endl;
    
    std::cout << "\n1. R-INDEX:" << std::endl;
    std::cout << "   Load time: " << profiling.rindex_load_time << " seconds" << std::endl;
    std::cout << "   Memory usage: " << profiling.rindex_load_memory_mb << " MB" << std::endl;
    
    std::cout << "\n2. TAG INDEX:" << std::endl;
    std::cout << "   Load time: " << profiling.tag_index_load_time << " seconds" << std::endl;
    std::cout << "   Memory usage: " << profiling.tag_index_load_memory_mb << " MB" << std::endl;
    
    std::cout << "\n3. READS:" << std::endl;
    std::cout << "   Total number of reads: " << format_number(profiling.total_reads) << std::endl;
    std::cout << "   Time taken to process all reads: " << profiling.total_reads_processing_time << " seconds" << std::endl;
    std::cout << "   Total number of MEMs outputted: " << format_number(profiling.total_mems_outputted) << std::endl;
    if (profiling.total_reads > 0) {
        std::cout << "   Average MEMs per read: " << (double)profiling.total_mems_outputted / profiling.total_reads << std::endl;
    }
    if (profiling.total_mems_outputted > 0) {
        std::cout << "   Average MEM length: " << (double)profiling.total_mem_length / profiling.total_mems_outputted << " bp" << std::endl;
        std::cout << "   Average MEM occurrence (mem.size): " << (double)profiling.total_mem_occurrences / profiling.total_mems_outputted << std::endl;
        std::cout << "   Mean per-MEM n/r ratio (n=mem.size, r=tag runs): " << profiling.total_n_over_r_ratio / profiling.total_mems_outputted << " (mean of per-MEM ratios; see Global n/r below for volume-weighted aggregate)" << std::endl;
    }
    std::cout << "   Time for finding all MEMs: " << profiling.total_mem_finding_time << " seconds" << std::endl;
    std::cout << "   Time for processing all MEMs: " << profiling.total_mem_processing_time << " seconds" << std::endl;
    if (profiling.total_reads > 0) {
        // std::cout << "   Average time for finding MEMs per read: " << profiling.total_mem_finding_time / profiling.total_reads << " seconds" << std::endl;
        // std::cout << "   Average time for processing MEMs per read: " << profiling.total_mem_processing_time / profiling.total_reads << " seconds" << std::endl;
    }
    
    std::cout << "\n4. PER-MEM OPERATIONS (aggregated across all MEMs):" << std::endl;
    std::cout << "   Total time for tag queries: " << profiling.total_tag_query_time << " seconds" << std::endl;
    std::cout << "   Total number of tag runs: " << format_number(profiling.total_tag_runs) << std::endl;
    std::cout << "   Total occurrence/bwt positions scanned: " << format_number(profiling.total_mem_occurrences) << std::endl;
    if (profiling.total_mems_outputted > 0) {
        std::cout << "   Average tag runs per MEM: " << (double)profiling.total_tag_runs / profiling.total_mems_outputted << std::endl;
    }
    // Global n/r ratio: total BWT positions scanned / total tag runs decoded.
    // Unlike the per-MEM mean(n/r) above (which can be skewed by a long tail
    // of high-occurrence MEMs), this is the volume-weighted aggregate density
    // of MEM hits per decoded tag run. It is the right comparison metric
    // across datasets of different MEM-occurrence distributions.
    if (profiling.total_tag_runs > 0) {
        std::cout << "   Global n/r ratio (Σmem.size / Σtag_runs): "
                  << (double)profiling.total_mem_occurrences / profiling.total_tag_runs << std::endl;
    }
    std::cout << "   Total time for locate operations: " << profiling.total_locate_time << " seconds" << std::endl;
    std::cout << "     - First locate (locate_sa_value): " << profiling.total_first_locate_time << " seconds" << std::endl;
    std::cout << "     - Locate next (locateNext calls): " << profiling.total_locate_next_time << " seconds" << std::endl;
    std::cout << "   Total number of locate operations: " << format_number(profiling.total_locate_operations) << std::endl;
    std::cout << "     - First locate calls: " << format_number(profiling.total_first_locate_calls) << std::endl;
    std::cout << "     - Locate next calls: " << format_number(profiling.total_locate_next_calls) << std::endl;
    if (profiling.total_mems_outputted > 0) {
        std::cout << "   Average locate operations per MEM: " << (double)profiling.total_locate_operations / profiling.total_mems_outputted << std::endl;
    }
    std::cout << "   Total time for file writes: " << profiling.total_file_write_time << " seconds" << std::endl;
    std::cout << "   Total number of entries written: " << format_number(profiling.total_entries_written) << std::endl;
    if (profiling.total_mems_outputted > 0) {
        std::cout << "   Average entries written per MEM: " << (double)profiling.total_entries_written / profiling.total_mems_outputted << std::endl;
    }
    
    std::cout << "\n5. SORTING:" << std::endl;
    std::cout << "   Total time for sorting (not included in read processing): " << profiling.total_sort_time << " seconds" << std::endl;
    
    std::cout << "\n6. MEMORY:" << std::endl;
    std::cout << "   Peak memory usage: " << profiling.peak_memory_mb << " MB" << std::endl;
    
    std::cout << "\n================================================" << std::endl;
    std::cout << "=== TIMING BREAKDOWN ===" << std::endl;
    std::cout << "================================================" << std::endl;
    double total_time = profiling.rindex_load_time + profiling.tag_index_load_time + 
                        profiling.total_reads_processing_time + profiling.total_sort_time;
    std::cout << "Total execution time: " << total_time << " seconds" << std::endl;
    std::cout << "  - R-index loading: " << profiling.rindex_load_time << " s" << std::endl;
    std::cout << "  - Tag index loading: " << profiling.tag_index_load_time << " s" << std::endl;
    std::cout << "  - Read processing: " << profiling.total_reads_processing_time << " s" << std::endl;
    std::cout << "  - Sorting: " << profiling.total_sort_time << " s" << std::endl;
    
    std::cout << "\n=== READ PROCESSING BREAKDOWN (% of read processing time) ===" << std::endl;
    if (profiling.total_reads_processing_time > 0) {
        std::cout << "  - MEM finding: " << profiling.total_mem_finding_time << " s (" 
                  << (profiling.total_mem_finding_time / profiling.total_reads_processing_time * 100) << "%)" << std::endl;
        std::cout << "  - MEM processing: " << profiling.total_mem_processing_time << " s (" 
                  << (profiling.total_mem_processing_time / profiling.total_reads_processing_time * 100) << "%)" << std::endl;
        std::cout << "    - Tag queries: " << profiling.total_tag_query_time << " s (" 
                  << (profiling.total_tag_query_time / profiling.total_mem_processing_time * 100) << "%)" << std::endl;
        std::cout << "    - Locate operations: " << profiling.total_locate_time << " s (" 
                  << (profiling.total_locate_time / profiling.total_mem_processing_time * 100) << "%)" << std::endl;
        std::cout << "      - First locate: " << profiling.total_first_locate_time << " s (" 
                  << (profiling.total_first_locate_time / profiling.total_locate_time * 100) << "%)" << std::endl;
        std::cout << "      - Locate next: " << profiling.total_locate_next_time << " s (" 
                  << (profiling.total_locate_next_time / profiling.total_locate_time * 100) << "%)" << std::endl;
        std::cout << "    - File writes: " << profiling.total_file_write_time << " s (" 
                  << (profiling.total_file_write_time / profiling.total_mem_processing_time * 100) << "%)" << std::endl;
        
        double other_time = profiling.total_reads_processing_time - 
                           (profiling.total_mem_finding_time + profiling.total_mem_processing_time);
        std::cout << "  - Other overhead: " << other_time << " s (" 
                  << (other_time / profiling.total_reads_processing_time * 100) << "%)" << std::endl;
    }
    std::cout << "================================================" << std::endl;
#endif
}
