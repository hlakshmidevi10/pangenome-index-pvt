#include "pangenome_index/algorithm.hpp"
#include "pangenome_index/tag_arrays.hpp"
#include <chrono>
#include <unordered_set>

using namespace std;
using namespace gbwtgraph;
using namespace panindexer;

#ifndef TIME
#define TIME 1
#endif

const std::string SEPARATOR = "\t";


// TODO Optimizations:
//  1. remove SeqID counter - instead compute running length when processing final output
//  2. Optimize for large no. of reads / streaming reads scenario - sorting n chunks and merging n sorted chunks
//  3. total entries at the start of the file to efficiently reserve space
//  4. Encode all entries in bytes instead of csv

// Helper function to dump MEM information
void dump_mem_info(const MEM& mem, const int read_id, TagArray& tag_array, FastLocate& r_index,
    vector<size_t> &seq_id_counter, std::ofstream* output_file = nullptr) {
    // Get BWT range for this MEM
    size_t bwt_start = mem.bwt_start;
    size_t bwt_end = mem.bwt_start + mem.size - 1;

    // Get sequence IDs from r_index.locate()
    gbwt::range_type range(bwt_start, bwt_end);
    std::vector<size_type> sa_values = r_index.locate(range);

    // Query tag array to get the runs for this BWT interval
    std::vector<std::pair<pos_t, uint16_t>> decoded_runs;
    tag_array.query_compressed_decoded_runs(mem.bwt_start, mem.bwt_start + mem.size - 1, decoded_runs);

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
        // std::string strand = is_reverse ? "-" : "+";

        // Output for each BWT position covered by this run
        for (uint16_t pos_in_run = 0; pos_in_run < run_length && sa_index < sa_values.size(); pos_in_run++) {
            // seq_id from r_index.locate()
            size_type seq_id = sa_values[sa_index];

            // Only output if we haven't seen this seq_id before // TODO: address this (choose node the longest run length)
            if (seen_seq_ids.find(seq_id) == seen_seq_ids.end()) {
                seen_seq_ids.insert(seq_id);
                seq_id_counter[seq_id]++;

                // Create output string
                std::string output_line = std::to_string(seq_id) + SEPARATOR +
                                        std::to_string(node_id) + SEPARATOR +
                                        std::to_string(offset) + SEPARATOR +
                                        std::to_string(is_reverse) + SEPARATOR +
                                        std::to_string(mem.end - mem.start) + SEPARATOR +
                                        std::to_string(mem.start) + SEPARATOR +
                                        std::to_string(read_id);


                // Output to stdout
                std::cout << output_line << std::endl;

                // Output to file if provided
                if (output_file && output_file->is_open()) {
                    *output_file << output_line << std::endl;
                }
            }

            sa_index++;
        }
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

    // Reserve space based on seq_id_counter if available
    // if (!seq_id_counter.empty()) {
    //     size_t estimated_entries = 0;
    //     for (const auto& pair : seq_id_counter) {
    //         estimated_entries += pair.second;
    //     }
    //     mem_entries.reserve(estimated_entries);
    // }

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


    // TODO: change this
    size_type current_seq_id = -1; // Initialize to invalid value
    size_t cumulative_count = 0;

    for (size_t curr_seq_id = 0; curr_seq_id < seq_id_counter.size(); curr_seq_id++) {
        // seq_id_starts_file << curr_seq_id << SEPARATOR << cumulative_count << "\n";
        seq_id_starts_file << cumulative_count << "\n";
        cumulative_count += seq_id_counter[curr_seq_id];
    }
    // n + 1
    seq_id_starts_file << cumulative_count << "\n";

    // for (const auto& entry : mem_entries) {
    //     output << entry.line << "\n";       // TODO: line might already contain \n?
    //
    //     // Write seq_id info only when we encounter a new seq_id (since entries are sorted)
    //     if (entry.seq_id != current_seq_id) {
    //         current_seq_id = entry.seq_id;
    //         // auto it = seq_id_counter.find(current_seq_id);
    //         // size_t count = (it != seq_id_counter.end()) ? it->second : 0;
    //         seq_id_starts_file << current_seq_id << SEPARATOR << cumulative_count << "\n";
    //     }
    //     cumulative_count++;
    // }

    // n + 1
    // seq_id_starts_file << current_seq_id << SEPARATOR << cumulative_count << "\n";

    output.close();
    seq_id_starts_file.close();
    std::cerr << "Successfully wrote " << mem_entries.size() << " sorted entries to " << output_file << std::endl;
}


int main(int argc, char **argv) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <r_index_file> <tag_array_index> <reads_file> <mem_length> <min_occ> [output_file]" << std::endl;
        return EXIT_FAILURE;
    }

    string r_index_file = argv[1];
    string tag_array_index = argv[2];
    string reads_file = argv[3];
    size_t mem_length = std::stoi(argv[4]);
    size_t min_occ = std::stoi(argv[5]);
    string output_file_name, tmp_file_name;


    // Optional output file parameter
    std::ofstream output_file;
    if (argc > 6) {
        output_file_name = argv[6];
        tmp_file_name = output_file_name + "_tmp.tsv";
        output_file.open(tmp_file_name);
        if (!output_file.is_open()) {
            std::cerr << "Cannot open output file: " << argv[6] << std::endl;
            return EXIT_FAILURE;
        }
    }

#if TIME
    auto time1 = chrono::high_resolution_clock::now();
    double total_mem_time = 0.0;
    double total_tag_time = 0.0;
#endif

    cerr << "Reading the rindex file" << endl;
    FastLocate r_index;
    if (!sdsl::load_from_file(r_index, r_index_file)) {
        std::cerr << "Cannot load the r-index from " << r_index_file << std::endl;
        std::exit(EXIT_FAILURE);
    }



//    std::cerr << "sym map " << (int) r_index.sym_map[NENDMARKER] << endl;
//    std::cerr << "sym map " << (int) r_index.sym_map['A'] << endl;
//    std::cerr << "sym map " << (int) r_index.sym_map['C'] << endl;
//    std::cerr << "sym map " << (int) r_index.sym_map['G'] << endl;
//    std::cerr << "sym map " << (int) r_index.sym_map['T'] << endl;
//    std::cerr << "sym map " << (int) r_index.sym_map['N'] << endl;
//
//    // r_index.initialize_complement_table();
//
//     FastLocate::bi_interval bint = {0, 0, r_index.bwt_size()};
//     auto forw = r_index.forward_extend(bint, 'A');
//     forw = r_index.forward_extend(forw, 'C');
//     forw = r_index.forward_extend(forw, 'A');
//     forw = r_index.forward_extend(forw, 'G');
//     forw = r_index.forward_extend(forw, 'G');
//     forw = r_index.forward_extend(forw, 'G');
//     forw = r_index.forward_extend(forw, 'G');
//     forw = r_index.forward_extend(forw, 'G');
//     forw = r_index.forward_extend(forw, 'G');
//     forw = r_index.forward_extend(forw, 'G');
//     forw = r_index.forward_extend(forw, 'G');
//
//
//     // print forw
//     std::cerr << "forward: " << forw.forward << " size: " << forw.size << " reverse: " << forw.reverse << std::endl;
//     auto back = r_index.backward_extend(bint, 'G');
//
//     back = r_index.backward_extend(back, 'A');
//     back = r_index.backward_extend(back, 'C');
//     back = r_index.backward_extend(back, 'A');
//     std::cerr << "backward: " << back.forward << " size: " << back.size << " reverse: " << back.reverse << std::endl;
//     exit(0);


#if TIME
    auto time2 = chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration1 = time2 - time1;
    std::cerr << "Loading r-index into memory took " << duration1.count() << " seconds" << std::endl;
#endif

    cerr << "Reading the tag array index" << endl;
    TagArray tag_array;
    std::ifstream in_ds(tag_array_index);
    tag_array.load_compressed_tags(in_ds);

#if TIME
    auto time3 = chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration2 = time3 - time2;
    std::cerr << "Loading tag arrays took " << duration2.count() << " seconds" << std::endl;
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

    while (std::getline(reads, read)) {
        if (read.empty()) continue;
        i++;

#if TIME
        auto time4 = chrono::high_resolution_clock::now();
#endif

        // Find MEMs for this read
        auto mems = find_all_mems(read, mem_length, min_occ, r_index);

#if TIME
        auto time5 = chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration3 = time5 - time4;
        total_mem_time += duration3.count();
//        std::cerr << "Finding MEMs took " << duration3.count() << " seconds" << std::endl;
#endif

        // std::vector<size_type> sa = r_index.decompressSA();
        // Print results for this read
        std::cout << "Read: " << i << std::endl;
        for (const auto& mem : mems) {
            std::cout << "Read: " << read << std::endl;
            std::cout << "MEM START: " << mem.start << ", MEM END: " << mem.end << " BWT START: " << mem.bwt_start << " BWT SIZE: " << mem.size << std::endl;
//            std::cout << "BWT interval start: " << mem.bi_interval.forward << ", size: " << mem.bi_interval.size << std::endl;

            total_mem_matches += mem.size;

            dump_mem_info(mem, i, tag_array, r_index, seq_id_counter,  output_file.is_open() ? &output_file : nullptr);
            size_t tag_nums = 0;

#if TIME
            auto time6 = chrono::high_resolution_clock::now();
#endif

            // Query the tag array for this MEM
            tag_array.query_compressed(mem.bwt_start, mem.bwt_start + mem.size - 1, tag_nums);

#if TIME
            auto time7 = chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration4 = time7 - time6;
            total_tag_time += duration4.count();
#endif
        }
        std::cout << std::endl;
    }

    // size_t tag_nums = 0;
    // tag_array.query_compressed(20, 24, tag_nums);
    // tag_array.query_compressed(1, 1 + 7, tag_nums);
    // tag_array.query_compressed(7717, 7717 + 3, tag_nums);
    // tag_array.query_compressed(7720, 7720 + 1, tag_nums);
    // tag_array.query_compressed(7720, 7720 + 2, tag_nums);


    reads.close();
    if (output_file.is_open()) {
        output_file.close();
    }

    if (argc > 6) {
        std::cerr << "seq_id_counter vector contents:" << std::endl;
        for (size_t k = 0; k < seq_id_counter.size(); k++) {
            if (seq_id_counter[k] > 0) {  // Only print non-zero entries
                std::cerr << "seq_id[" << k << "] = " << seq_id_counter[k] << std::endl;
            }
        }
        sort_mem_output_by_seq_node(tmp_file_name, output_file_name, seq_id_counter, total_mem_matches);
    }

#if TIME
    std::cout << "\nTotal time for finding all MEMs: " << total_mem_time << " seconds" << std::endl;
    std::cout << "Total time for all tag queries: " << total_tag_time << " seconds" << std::endl;
#endif
}
