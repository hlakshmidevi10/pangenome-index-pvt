#include "pangenome_index/algorithm.hpp"
#include "pangenome_index/tag_arrays.hpp"
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

// Profiling structure for per-MEM statistics
struct MEMProfilingStats {
    double tag_query_time = 0.0;
    double locate_time = 0.0;
    double write_time = 0.0;
    double file_write_time = 0.0;  // Time spent specifically on file I/O
    size_t unique_runs = 0;
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
    
    // MEM statistics
    size_t total_mem_length = 0;      // Sum of all MEM lengths
    size_t total_mem_occurrences = 0;  // Sum of all MEM occurrences (mem.size)
    double total_n_over_r_ratio = 0.0; // Sum of n/r ratios (n=mem.size, r=decoded_runs.size)
    
    // Per-MEM operations (aggregated)
    double total_tag_query_time = 0.0;
    double total_locate_time = 0.0;
    double total_write_time = 0.0;
    double total_file_write_time = 0.0;  // Time spent specifically on file I/O
    size_t total_unique_runs = 0;
    
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

// Helper function to dump MEM information
void dump_mem_info(const MEM& mem, const int read_id, TagArray& tag_array, FastLocate& r_index,
    vector<size_t> &seq_id_counter, MEMProfilingStats& mem_stats,
    std::ofstream* output_file = nullptr, bool debug_stats = false) {

    // Get BWT range for this MEM
    size_t bwt_start = mem.bwt_start;
    size_t bwt_end = mem.bwt_start + mem.size - 1;

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
    mem_stats.unique_runs = decoded_runs.size();
#endif

#if TIME
    auto locate_start = chrono::high_resolution_clock::now();
#endif

    // Get sequence IDs from r_index.locate()
    gbwt::range_type range(bwt_start, bwt_end);
    std::vector<size_type> sa_values = r_index.locate(range);   // TODO: Address this expensive operation

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

#if TIME
                auto file_write_start = chrono::high_resolution_clock::now();
#endif
                // Create output string
                std::string output_line = std::to_string(seq_id) + SEPARATOR +
                                        std::to_string(node_id) + SEPARATOR +
                                        std::to_string(offset) + SEPARATOR +
                                        std::to_string(is_reverse) + SEPARATOR +
                                        std::to_string(mem.end - mem.start) + SEPARATOR +
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

#if TIME
    auto write_end = chrono::high_resolution_clock::now();
    std::chrono::duration<double> write_duration = write_end - write_start;
    mem_stats.write_time = write_duration.count();
#endif

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
void sort_mem_output_by_seq_node(const std::string& input_file, const std::string& output_file,
                                  const vector<size_t> &seq_id_counter, size_t total_mem_matches) {
    std::cerr << "Sorting MEM output file: " << input_file << std::endl;

    // Structure to hold parsed CSV data
    struct MEMData {
        size_type seq_id;
        size_type node_id;
        std::string line;

        // Comparison operator for sorting
        bool operator<(const MEMData& other) const {    // TODO: Avoid global sort; instead sort only by node_id - focusin on chunks
            if (seq_id != other.seq_id) {
                return seq_id < other.seq_id;
            }
            return node_id < other.node_id;
        }
    };

    std::ifstream input(input_file);
    if (!input.is_open()) {
        throw std::runtime_error("Cannot open input file: " + input_file);
    }

    std::vector<MEMData> mem_entries;

    mem_entries.reserve(total_mem_matches);     // TODO: change

    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;

        MEMData entry;

        // Parse seq_id and node_id from CSV (first two columns)
        size_t first_separator = line.find(SEPARATOR);
        if (first_separator == std::string::npos) continue;

        entry.line = line.substr(first_separator + 1);

        size_t second_separator = line.find(SEPARATOR, first_separator + 1);
        if (second_separator == std::string::npos) continue;

        try {
            entry.seq_id = std::stoull(line.substr(0, first_separator));
            entry.node_id = std::stoull(line.substr(first_separator + 1, second_separator - first_separator - 1));
            mem_entries.push_back(entry);
        } catch (const std::exception& e) {
            std::cerr << "Warning: Could not parse line: " << line << std::endl;
        }
    }

    input.close();
    std::cerr << "Read " << mem_entries.size() << " entries" << std::endl;

    // Sort the entries
    std::cerr << "Sorting entries..." << std::endl;
    std::sort(mem_entries.begin(), mem_entries.end());

    // Write sorted entries to output file
    std::ofstream output(output_file + "_path_pos.tsv");
    if (!output.is_open()) {
        throw std::runtime_error("Cannot open output file: " + output_file + "_path_pos.tsv");
    }

    for (const auto & entry : mem_entries) {
        output << entry.line << "\n";       // TODO: line might already contain \n?
    }

    // Create seq_id_starts files
    std::ofstream seq_id_starts_file(output_file + "_seq_id_starts.out");
    if (!seq_id_starts_file.is_open()) {
        throw std::runtime_error("Cannot open seq_id_starts file: " + output_file + "_seq_id_starts.out");
    }


    size_t cumulative_count = 0;

    for (size_t curr_seq_id = 0; curr_seq_id < seq_id_counter.size(); curr_seq_id++) {
        // seq_id_starts_file << curr_seq_id << SEPARATOR << cumulative_count << "\n";
        seq_id_starts_file << cumulative_count << "\n";
        cumulative_count += seq_id_counter[curr_seq_id];
    }
    // n + 1
    seq_id_starts_file << cumulative_count << "\n";

    assert(cumulative_count == mem_entries.size());

    output.close();
    seq_id_starts_file.close();
    std::cerr << "Successfully wrote " << mem_entries.size() << " sorted entries to " << output_file << std::endl;
}


int main(int argc, char **argv) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <r_index_file> <tag_array_index> <reads_file> <mem_length> <min_occ> [output_file] [--debug-stats] [--verbose]" << std::endl;
        return EXIT_FAILURE;
    }

    string r_index_file = argv[1];
    string tag_array_index = argv[2];
    string reads_file = argv[3];
    size_t mem_length = std::stoi(argv[4]);
    size_t min_occ = std::stoi(argv[5]);
    string output_file_name, tmp_file_name;
    bool debug_stats = false;
    bool verbose = false;

    // Check for debug flags
    for (int i = 6; i < argc; i++) {
        if (std::string(argv[i]) == "--debug-stats") {
            debug_stats = true;
        }
        if (std::string(argv[i]) == "--verbose" || std::string(argv[i]) == "--debug") {
            verbose = true;
        }
    }

    // Optional output file parameter
    std::ofstream output_file;
    if (argc > 6) {
        bool found_output_file = false;
        for (int i = 6; i < argc; i++) {
            if (std::string(argv[i]) != "--debug-stats" && std::string(argv[i]) != "--verbose" && std::string(argv[i]) != "--debug") {
                output_file_name = argv[i];
                tmp_file_name = output_file_name + "_tmp.tsv";
                output_file.open(tmp_file_name);
                if (!output_file.is_open()) {
                    std::cerr << "Cannot open output file: " << argv[i] << std::endl;
                    return EXIT_FAILURE;
                }
                found_output_file = true;
                break;
            }
        }
    }

    // Initialize profiling data
    ProfilingData profiling;
    double initial_memory_mb = get_memory_usage_mb();

#if TIME
    auto rindex_load_start = chrono::high_resolution_clock::now();
#endif

    cerr << "Reading the rindex file" << endl;
    FastLocate r_index;
    if (!sdsl::load_from_file(r_index, r_index_file)) {
        std::cerr << "Cannot load the r-index from " << r_index_file << std::endl;
        std::exit(EXIT_FAILURE);
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

    cerr << "Reading the tag array index" << endl;
    TagArray tag_array;
    std::ifstream in_ds(tag_array_index);
    tag_array.load_compressed_tags(in_ds);

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
            
            dump_mem_info(mem, i, tag_array, r_index, seq_id_counter, mem_stats,
                output_file.is_open() ? &output_file : nullptr, debug_stats);

            // Aggregate per-MEM stats
            profiling.total_tag_query_time += mem_stats.tag_query_time;
            profiling.total_locate_time += mem_stats.locate_time;
            profiling.total_write_time += mem_stats.write_time;
            profiling.total_file_write_time += mem_stats.file_write_time;
            profiling.total_unique_runs += mem_stats.unique_runs;
            
            // Calculate n/r ratio (n = mem.size, r = number of decoded runs)
            if (mem_stats.unique_runs > 0) {
                double n_over_r = (double)mem.size / mem_stats.unique_runs;
                profiling.total_n_over_r_ratio += n_over_r;
            }
        }
        
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
    if (output_file.is_open()) {
        output_file.close();
    }

    if (argc > 6) {
#if TIME
        auto sort_start = chrono::high_resolution_clock::now();
#endif

        sort_mem_output_by_seq_node(tmp_file_name, output_file_name, seq_id_counter, total_mem_matches);

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
    std::cout << "   Total number of reads: " << profiling.total_reads << std::endl;
    std::cout << "   Time taken to process all reads: " << profiling.total_reads_processing_time << " seconds" << std::endl;
    std::cout << "   Total number of MEMs outputted: " << profiling.total_mems_outputted << std::endl;
    if (profiling.total_reads > 0) {
        std::cout << "   Average MEMs per read: " << (double)profiling.total_mems_outputted / profiling.total_reads << std::endl;
    }
    if (profiling.total_mems_outputted > 0) {
        std::cout << "   Average MEM length: " << (double)profiling.total_mem_length / profiling.total_mems_outputted << " bp" << std::endl;
        std::cout << "   Average MEM occurrence (mem.size): " << (double)profiling.total_mem_occurrences / profiling.total_mems_outputted << std::endl;
        std::cout << "   Average n/r ratio (n=mem.size, r=tag runs): " << profiling.total_n_over_r_ratio / profiling.total_mems_outputted << std::endl;
    }
    std::cout << "   Time for finding all MEMs: " << profiling.total_mem_finding_time << " seconds" << std::endl;
    if (profiling.total_reads > 0) {
        std::cout << "   Average time for finding MEMs per read: " << profiling.total_mem_finding_time / profiling.total_reads << " seconds" << std::endl;
    }
    
    std::cout << "\n4. PER-MEM OPERATIONS (aggregated across all MEMs):" << std::endl;
    std::cout << "   Total time for tag queries: " << profiling.total_tag_query_time << " seconds" << std::endl;
    std::cout << "   Total number of unique runs: " << profiling.total_unique_runs << std::endl;
    if (profiling.total_mems_outputted > 0) {
        std::cout << "   Average unique runs per MEM: " << (double)profiling.total_unique_runs / profiling.total_mems_outputted << std::endl;
        std::cout << "   Average tag query time per MEM: " << profiling.total_tag_query_time / profiling.total_mems_outputted << " seconds" << std::endl;
    }
    std::cout << "   Total time for locate operations: " << profiling.total_locate_time << " seconds" << std::endl;
    if (profiling.total_mems_outputted > 0) {
        std::cout << "   Average locate time per MEM: " << profiling.total_locate_time / profiling.total_mems_outputted << " seconds" << std::endl;
    }
    std::cout << "   Total time for writing output (overall): " << profiling.total_write_time << " seconds" << std::endl;
    std::cout << "   Total time for file I/O only: " << profiling.total_file_write_time << " seconds" << std::endl;
    if (profiling.total_mems_outputted > 0) {
        std::cout << "   Average write time per MEM: " << profiling.total_write_time / profiling.total_mems_outputted << " seconds" << std::endl;
        std::cout << "   Average file I/O time per MEM: " << profiling.total_file_write_time / profiling.total_mems_outputted << " seconds" << std::endl;
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
        std::cout << "  - Tag queries: " << profiling.total_tag_query_time << " s (" 
                  << (profiling.total_tag_query_time / profiling.total_reads_processing_time * 100) << "%)" << std::endl;
        std::cout << "  - Locate operations: " << profiling.total_locate_time << " s (" 
                  << (profiling.total_locate_time / profiling.total_reads_processing_time * 100) << "%)" << std::endl;
        std::cout << "  - Writing output (overall): " << profiling.total_write_time << " s (" 
                  << (profiling.total_write_time / profiling.total_reads_processing_time * 100) << "%)" << std::endl;
        std::cout << "    - File I/O only: " << profiling.total_file_write_time << " s (" 
                  << (profiling.total_file_write_time / profiling.total_reads_processing_time * 100) << "%)" << std::endl;
        
        double other_time = profiling.total_reads_processing_time - 
                           (profiling.total_mem_finding_time + profiling.total_tag_query_time + 
                            profiling.total_locate_time + profiling.total_write_time);
        std::cout << "  - Other overhead: " << other_time << " s (" 
                  << (other_time / profiling.total_reads_processing_time * 100) << "%)" << std::endl;
    }
    std::cout << "================================================" << std::endl;
#endif
}
