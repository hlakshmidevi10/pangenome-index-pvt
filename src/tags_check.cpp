//
// Created by seeskand on 3/4/25.
//

#include "pangenome_index/r-index.hpp"
#include "pangenome_index/algorithm.hpp"
#include "pangenome_index/tag_arrays.hpp"
#include "pangenome_index/sampled_tag_array.hpp"
#include <gbwtgraph/utils.h>
#include <gbwt/gbwt.h>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/gbz.h>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <omp.h>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <utility>
#include <stdexcept>
#include <string>
#include <cstring>


#ifndef TIME
#define TIME 1
#endif






//using namespace gbwtgraph;
using namespace panindexer;
using namespace std;
using namespace gbwtgraph;
using handlegraph::pos_t;

namespace fs = std::filesystem;





class FileReader {
public:


    FileReader(const std::vector <std::string> &files, size_t n_threads, size_t batch_size)
            : files(files), n_threads(n_threads), current_thread_id(0), batch_size(batch_size) {
        if (files.empty()) {
            throw std::invalid_argument("File list cannot be empty.");
        }

        // Initialize per-file mutexes
        for (size_t i = 0; i < files.size(); ++i) {
            mutexes.emplace_back(std::make_unique<std::mutex>());
        }

        // Initialize file data
        file_positions.resize(files.size(), 0);
        file_end_reached.resize(files.size(), false);
        file_buffers.resize(files.size());
        tag_batches.resize(files.size());
        current_index_tag_batch.resize(files.size(), 0);
        current_length_tag_batch.resize(files.size(), 0);
        for (size_t i = 0; i < files.size(); ++i) {
            tag_batches[i].resize(batch_size);
        }

        initializeFiles();

    }

    void print_all_tags(){
        for (size_t i = 0; i < files.size(); i++){
            size_t tag_tot = 0;
            std::cerr << "File: " << files[i] << std::endl;
            for (size_t j = 0; j < current_length_tag_batch[i]; j++){
                std::cerr << "Tag: " << tag_batches[i][j].first << " Length: " << int(tag_batches[i][j].second) << std::endl;
                tag_tot += tag_batches[i][j].second;
            }
            std::cerr << "=============Total tags: " << tag_tot << std::endl;
        }
    }

    int get_file_num(){
        return files.size();
    }


    pos_t get_first_tag(int fileIndex) {
//        std::cerr << "Starting tag runs of file " << files[fileIndex] << std::endl;
//        for (size_t i = 0; i < 3; i++) {
//            std::cerr << tag_batches[fileIndex][i].first << " " << int(tag_batches[fileIndex][i].second) << std::endl;
//        }
        return tag_batches[fileIndex][0].first;
    }

    void print_first_n_item(int n, int fileIndex){
        for (int i = 0; i < n; i++){
            std::cerr << "Tag: " << tag_batches[fileIndex][i].first << " Length: " << int(tag_batches[fileIndex][i].second) << std::endl;
        }
    }



    void closeAllFiles() {
        std::cerr << "Closing all files" << std::endl;
        for (size_t i = 0; i < file_buffers.size(); i++) {
            if (file_buffers[i].is_open()) {
                file_buffers[i].close(); // Explicitly close the file buffer
                std::cerr << "Closed file: " << files[i] << std::endl;
            }
        }
        file_buffers.clear(); // Optionally clear the vector to release all resources
        file_positions.clear(); // Clear positions as files are closed
        file_end_reached.clear(); // Clear end-of-file flags
        tag_batches.clear(); // Clear tag buffersStarting creating the request list
    }

    pos_t get_next_tag(int fileIndex) {
        if (current_length_tag_batch[fileIndex] == 0) {
            std::cerr << "GET NEXT TAG: No more tags is possible to read from file: " << files[fileIndex] << std::endl;
        }

        pos_t res = tag_batches[fileIndex][current_index_tag_batch[fileIndex]].first;
        tag_batches[fileIndex][current_index_tag_batch[fileIndex]].second--;

        if (tag_batches[fileIndex][current_index_tag_batch[fileIndex]].second == 0){
            current_index_tag_batch[fileIndex] = (current_index_tag_batch[fileIndex] + 1) % batch_size;
            current_length_tag_batch[fileIndex]--;
        }

        return res;
    }





    std::vector<std::vector<std::pair<pos_t, uint8_t>>> extract_requested_tags(
            size_t thread_id, const std::vector<size_t>& requests) {

        waitForTurn(thread_id);

        std::vector<std::vector<std::pair<pos_t, uint8_t>>> extracted_tags(files.size());

        for (size_t i = 0; i < files.size(); i++) {
            size_t current_extracted = 0;


            while (current_extracted < requests[i]){
                if (current_length_tag_batch[i] == 0) {
                    if (file_end_reached[i]) {
                        std::cerr << "No more tags is possible to read from file: " << files[i] << std::endl;
                        std::cerr << "Needed " << requests[i] << " but only extracted " << current_extracted << std::endl;
                        break;
                    } else {
                        refill_tags();
                    }
                }
                if (current_extracted + tag_batches[i][current_index_tag_batch[i]].second <= requests[i]){
                    // add the whole tag run to the extracted tags and delete it from the tag_batches
                    extracted_tags[i].push_back(tag_batches[i][current_index_tag_batch[i]]);
                    current_extracted += tag_batches[i][current_index_tag_batch[i]].second;
                    current_index_tag_batch[i] = (current_index_tag_batch[i] + 1) % batch_size;
                    current_length_tag_batch[i]--;

                } else {
                    // add the first part of the tag run to the extracted tags and update the tag run in the tag_batches
                    extracted_tags[i].push_back(std::make_pair(tag_batches[i][current_index_tag_batch[i]].first, requests[i] - current_extracted));
                    tag_batches[i][current_index_tag_batch[i]].second -= (requests[i] - current_extracted);
                    current_extracted = requests[i];
                }

            }

        }


        refill_tags();

        notifyNext(thread_id);
        return extracted_tags;
    }



private:
    void initializeFiles() {
        std::cerr << "Initializing files" << std::endl;
        for (size_t i = 0; i < files.size(); i++) {
            sdsl::int_vector_buffer<8> in(files[i], std::ios::in);
            if (!in.is_open()) {
                throw std::runtime_error("Cannot open file: " + files[i]);
            }
            file_buffers[i] = std::move(in);


            // want to read batch_size of the tags from each file and store them in the tag_batches

            for (size_t j = 0; j < batch_size; j++) {
                if (file_positions[i] >= file_buffers[i].size()) {
                    std::cerr << "End of file reached for: " << files[i] << std::endl;
                    file_end_reached[i] = true;
                    break;
                }

                auto tag_block = panindexer::TagArray::decode_run(
                        gbwt::ByteCode::read(file_buffers[i], file_positions[i])
                );
                tag_batches[i][j] = tag_block;
                current_length_tag_batch[i]++;

            }

        }
    }



    // This function checks for all the tag files it has that if the current number of batch files are less than batch_size/3 it will refill the tags
    void refill_tags() {
        for (size_t i = 0; i < files.size(); i++) {
            if (current_length_tag_batch[i] < batch_size / 3 && !file_end_reached[i]) {
                size_t remaining_tags = current_length_tag_batch[i];
                std::vector <std::pair<pos_t, uint8_t>> new_tags;
                size_t new_tags_size = 0;
                while (new_tags_size + remaining_tags < batch_size) {
                    if (file_positions[i] >= file_buffers[i].size()) {
                        // End of file reached, stop reading
                        file_end_reached[i] = true;
                        std::cerr << "REFILL End of file reached for: " << files[i] << std::endl;
                        break;
                    }

                    auto tag_block = panindexer::TagArray::decode_run(
                            gbwt::ByteCode::read(file_buffers[i], file_positions[i])
                    );

                    tag_batches[i][(current_length_tag_batch[i] + current_index_tag_batch[i]) % batch_size] = tag_block;
                    current_length_tag_batch[i]++;
                    new_tags_size++;
                }
            }
        }
    }




    void waitForTurn(size_t thread_id) {
        std::unique_lock <std::mutex> lock(thread_mutex);
        thread_cv.wait(lock, [this, thread_id] {
            return thread_id == current_thread_id;
        });
    }

    void notifyNext(size_t thread_id) {
        std::lock_guard <std::mutex> lock(thread_mutex);

        current_thread_id++;
        if (current_thread_id == n_threads) {
            current_thread_id = 0;
        }

        thread_cv.notify_all();
    }

    std::vector <std::string> files;          // List of file paths
    std::vector <bool> file_end_reached;      // Flag to indicate end of file

    size_t n_threads;                        // Number of threads
    size_t current_thread_id;                // Current thread ID to execute
    std::vector <std::unique_ptr<std::mutex>> mutexes; // Mutex for each file
    std::mutex thread_mutex;                 // Mutex for thread coordination
    std::condition_variable thread_cv;       // Condition variable for thread coordination

    std::vector <gbwt::size_type> file_positions;      // Current file positions
    std::vector <sdsl::int_vector_buffer<8>> file_buffers; // File buffers for reading
    size_t batch_size;                      // Number of tags to read in a batch
    std::vector <size_t> current_index_tag_batch; // Current index of tag batch for each file
    std::vector <std::vector<std::pair < pos_t, uint8_t>>> tag_batches; // buffer of tags of each tag_file
    std::vector <size_t> current_length_tag_batch; // Remaining length of tag batch for each file

};








std::vector <std::string> get_files_in_dir(const std::string &directoryPath) {
    std::vector <std::string> files;
    if (!fs::is_directory(directoryPath)) {
        std::cerr << "Path is not a directory: " << directoryPath << std::endl;
        return files;
    }

    for (const auto &entry: fs::directory_iterator(directoryPath)) {
        if (fs::is_regular_file(entry.status())) {
            files.push_back(entry.path().string());
        }
    }

    return files;
}

// --- Verify sampled tag array against GBWT path ---
// Decode tag code to (node_id, is_rev) using same convention as SampledTagArray::encode_value
static std::pair<int64_t, bool> decode_tag(uint64_t tag_code) {
    uint64_t decoded = tag_code - 1;
    int64_t node_id = (decoded >> 1) + 1;
    bool is_rev = (decoded & 1) != 0;
    return {node_id, is_rev};
}

// Recover packed text position at a given BWT position using r-index (locate).
static size_t recover_text_pos_from_bwt(panindexer::FastLocate& r_index, size_t bwt_pos) {
    size_t run_id = 0, bwt_run_start = 0;
    r_index.run_id_and_offset_at(bwt_pos, run_id, bwt_run_start);
    size_t packed = r_index.getSample(run_id);
    for (size_t k = bwt_run_start; k < bwt_pos; ++k) {
        packed = r_index.locateNext(packed);
    }
    return packed;
}

// Get BWT position for a given packed text position using last_successor + LF backwards
static bool text_pos_to_bwt_pos(panindexer::FastLocate& r_index, size_t text_pos,
                                size_t seq_end_text_pos, size_t& out_bwt_pos) {
    r_index.ensure_last_rank();
    r_index.ensure_last_select();
    auto successor_result = r_index.last_successor(text_pos);
    size_t text_pos_x = successor_result.first;
    size_t rank_x = successor_result.second;
    if (rank_x >= r_index.last_to_run.size()) {
        return false;
    }
    if (text_pos_x > seq_end_text_pos) {
        return false;  // beyond sequence end
    }
    size_t run_id = r_index.last_to_run[rank_x];
    size_t current_bwt_pos = r_index.bwt_end_position_of_run(run_id);
    size_t current_text_pos = text_pos_x;
    while (current_text_pos > text_pos) {
        current_text_pos--;
        current_bwt_pos = r_index.LF(current_bwt_pos);
    }
    out_bwt_pos = current_bwt_pos;
    return true;
}

// Verify sampled tag array by traversing a GBWT path from the start and comparing
// the tag at each (seq_id, base_offset) with the node on the path at that offset.
static int verify_sampled_against_gbwt(const std::string& r_index_file,
                                       const std::string& sampled_tags_file,
                                       const std::string& gbwt_index_file,
                                       const std::string& gbz_file,
                                       size_t path_id,
                                       size_t sample_every) {
    using namespace panindexer;
    std::cerr << "Loading RLBWT r-index: " << r_index_file << std::endl;
    FastLocate r_index;
    {
        std::ifstream rin(r_index_file, std::ios::binary);
        if (!rin) {
            std::cerr << "Cannot open r-index file: " << r_index_file << std::endl;
            return 1;
        }
        r_index.load_encoded(rin);
    }
    std::cerr << "Loading sampled tag array: " << sampled_tags_file << std::endl;
    SampledTagArray sampled;
    {
        std::ifstream sin(sampled_tags_file, std::ios::binary);
        if (!sin) {
            std::cerr << "Cannot open sampled tags file: " << sampled_tags_file << std::endl;
            return 1;
        }
        sampled.load(sin);
    }
    size_t rindex_bwt_size = r_index.bwt_size();
    size_t sampled_bwt_size = (sampled.run_starts().size() > 0) ? (sampled.run_starts().size() - 1) : 0;
    if (rindex_bwt_size != sampled_bwt_size) {
        std::cerr << "WARNING: RLBWT BWT size (" << rindex_bwt_size
                  << ") != sampled tag array BWT size (" << sampled_bwt_size
                  << "). Verification may be wrong." << std::endl;
    }
    std::cerr << "Loading GBWT index: " << gbwt_index_file << std::endl;
    gbwt::GBWT gbwt_index;
    {
        std::ifstream gin(gbwt_index_file, std::ios::binary);
        if (!gin) {
            std::cerr << "Cannot open GBWT index: " << gbwt_index_file << std::endl;
            return 1;
        }
        sdsl::load_from_file(gbwt_index, gbwt_index_file);
    }
    std::cerr << "Loading GBZ (graph): " << gbz_file << std::endl;
    gbwtgraph::GBZ gbz;
    sdsl::simple_sds::load_from(gbz, gbz_file);
    const gbwtgraph::GBWTGraph& graph = gbz.graph;

    gbwt::vector_type path = gbwt_index.extract(gbwt::Path::encode(path_id, false));
    if (path.empty()) {
        std::cerr << "Path " << path_id << " not found in GBWT." << std::endl;
        return 1;
    }
    size_t path_length_bases = 0;
    for (gbwt::node_type node : path) {
        if (node == gbwt::ENDMARKER) break;
        path_length_bases += graph.get_length(
            graph.get_handle(gbwt::Node::id(node), gbwt::Node::is_reverse(node)));
    }
    std::cerr << "Path " << path_id << " has " << path.size() << " nodes, "
              << path_length_bases << " bases." << std::endl;

    // Sequence in RLBWT: assume path_id in GBWT corresponds to seq_id in RLBWT
    size_t seq_id = path_id;
    size_t seq_end_text_pos = r_index.pack(seq_id, path_length_bases);
    if (seq_end_text_pos > r_index.bwt_size()) {
        std::cerr << "Path length exceeds RLBWT sequence length; is seq_id correct?" << std::endl;
    }
    sampled.ensure_run_rank();
    sampled.ensure_run_select();
    const auto& wm = sampled.values();
    bool first_run_is_gap = sampled.is_first_run_gap();

    size_t cumulative_bases = 0;
    size_t node_index = 0;
    size_t checked = 0;
    size_t mismatches = 0;
    for (gbwt::node_type node : path) {
        if (node == gbwt::ENDMARKER) break;
        int64_t node_id = gbwt::Node::id(node);
        bool node_rev = gbwt::Node::is_reverse(node);
        gbwtgraph::handle_t handle = graph.get_handle(node_id, node_rev);
        size_t node_len = graph.get_length(handle);
        uint64_t expected_tag = SampledTagArray::encode_value(node_id, node_rev);

        // For this GBWT node, find all runs in the sampled tag array with this tag
        {
            size_t total_occ = wm.rank(wm.size(), expected_tag);
            std::cerr << "  [node_index=" << node_index << " node_id=" << node_id << " rev=" << node_rev
                      << " tag=" << expected_tag << "] In sampled tag array: " << total_occ << " run(s) with this tag:";
            if (total_occ == 0) {
                std::cerr << " (none)" << std::endl;
            } else {
                std::cerr << std::endl;
                for (size_t j = 1; j <= total_occ; ++j) {
                    size_t tag_run_id = wm.select(j, expected_tag);  // 0-based index in non-gap runs
                    size_t bwt_run_id = 2 * tag_run_id + 1 - static_cast<size_t>(first_run_is_gap);
                    auto span = sampled.run_span(bwt_run_id);
                    std::cerr << "      run " << j << "/" << total_occ << ": bwt_run_id=" << bwt_run_id
                              << "  BWT [" << span.first << ", " << span.second << "]"
                              << " (length " << (span.second - span.first + 1) << ")";
                    // Locate BWT start and end to get text positions (seq_id, offset)
                    if (span.first < r_index.bwt_size()) {
                        size_t packed_start = recover_text_pos_from_bwt(r_index, span.first);
                        auto [seq_start, offset_start] = r_index.unpack(packed_start);
                        std::cerr << "  text@start: packed=" << packed_start << " seq_id=" << seq_start << " offset=" << offset_start;
                        if (span.second < r_index.bwt_size() && span.second != span.first) {
                            size_t packed_end = recover_text_pos_from_bwt(r_index, span.second);
                            auto [seq_end, offset_end] = r_index.unpack(packed_end);
                            std::cerr << "  text@end: packed=" << packed_end << " seq_id=" << seq_end << " offset=" << offset_end;
                        }
                    }
                    std::cerr << std::endl;
                }
            }
        }

        for (size_t off_in_node = 0; off_in_node < node_len; off_in_node += sample_every) {
            size_t base_offset = cumulative_bases + off_in_node;
            size_t text_pos = r_index.pack(seq_id, base_offset);
            size_t bwt_pos = 0;
            if (!text_pos_to_bwt_pos(r_index, text_pos, seq_end_text_pos, bwt_pos)) {
                std::cerr << "  Skip base_offset=" << base_offset << " (could not get BWT position)" << std::endl;
                continue;
            }
            size_t run_id = sampled.run_id_at(bwt_pos);
            uint64_t tag_val = sampled.run_value(run_id);
            checked++;
            bool match = (tag_val == expected_tag);
            if (!match) mismatches++;

            // Print only at the beginning of each node (first base of the node)
            if (off_in_node == 0) {
                size_t packed_at = recover_text_pos_from_bwt(r_index, bwt_pos);
                auto [seq_id_at, offset_at] = r_index.unpack(packed_at);
                std::cerr << "  base_offset=" << base_offset
                          << " bwt_pos=" << bwt_pos
                          << " text: packed=" << packed_at << " seq_id=" << seq_id_at << " offset=" << offset_at;
                if (seq_id_at != path_id || offset_at != base_offset) {
                    std::cerr << " [expected seq_id=" << path_id << " offset=" << base_offset << "]";
                }
                std::cerr << "  node_index=" << node_index
                          << " expected: node_id=" << node_id << " rev=" << node_rev << " tag=" << expected_tag;
                if (tag_val == 0) {
                    std::cerr << "  sampled: gap (tag=0)";
                } else {
                    auto [got_id, got_rev] = decode_tag(tag_val);
                    std::cerr << "  sampled: tag=" << tag_val << " node_id=" << got_id << " rev=" << got_rev;
                }
                std::cerr << "  " << (match ? "OK" : "MISMATCH") << std::endl;
            }
        }
        cumulative_bases += node_len;
        node_index++;
    }
    std::cerr << "Verification done: " << checked << " positions checked, "
              << mismatches << " mismatches." << std::endl;
    return mismatches > 0 ? 1 : 0;
}

static void usage_verify(const char* prog) {
    std::cerr << "Usage: " << prog << " --verify-sampled <r_index.ri> <sampled.tags> --gbwt-index <graph.gbwt> --gbz <graph.gbz> [options]\n"
              << "  Verify that the sampled tag array matches the GBWT path by traversing from the start.\n"
              << "Options:\n"
              << "  --path-id N       GBWT path/sequence ID to traverse (default: 0)\n"
              << "  --sample-every N  Check every N bases along the path (default: 1)\n" << std::endl;
}

int main(int argc, char **argv) {
    if (argc >= 2 && std::strcmp(argv[1], "--verify-sampled") == 0) {
        std::string r_index_file;
        std::string sampled_tags_file;
        std::string gbwt_index_file;
        std::string gbz_file;
        size_t path_id = 0;
        size_t sample_every = 1;
        for (int i = 2; i < argc; i++) {
            if (std::strcmp(argv[i], "--gbwt-index") == 0 && i + 1 < argc) {
                gbwt_index_file = argv[++i];
            } else if (std::strcmp(argv[i], "--gbz") == 0 && i + 1 < argc) {
                gbz_file = argv[++i];
            } else if (std::strcmp(argv[i], "--path-id") == 0 && i + 1 < argc) {
                path_id = static_cast<size_t>(std::stoull(argv[++i]));
            } else if (std::strcmp(argv[i], "--sample-every") == 0 && i + 1 < argc) {
                sample_every = static_cast<size_t>(std::stoull(argv[++i]));
                if (sample_every == 0) sample_every = 1;
            } else if (argv[i][0] != '-') {
                if (r_index_file.empty()) r_index_file = argv[i];
                else if (sampled_tags_file.empty()) sampled_tags_file = argv[i];
            }
        }
        if (r_index_file.empty() || sampled_tags_file.empty() || gbwt_index_file.empty() || gbz_file.empty()) {
            usage_verify(argv[0]);
            return 1;
        }
        return verify_sampled_against_gbwt(r_index_file, sampled_tags_file,
                                          gbwt_index_file, gbz_file, path_id, sample_every);
    }

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <gbz_graph> <r_index.ri> <tag_array_index_dir>\n"
                  << "   Or: " << argv[0] << " --verify-sampled <r_index.ri> <sampled.tags> --gbwt-index <graph.gbwt> --gbz <graph.gbz> [--path-id N] [--sample-every N]" << std::endl;
        return 1;
    }

    std::string gbz_graph = std::string(argv[1]);
    std::string r_index_file = std::string(argv[2]);
    std::string tag_array_index_dir = std::string(argv[3]);
    int threads = 8;

    // GBZ gbz;
    // cerr << "Loading the graph file" << endl;
    // sdsl::simple_sds::load_from(gbz, gbz_graph);

    std::cerr << "Getting the lists of tag files" << std::endl;
    // get the list of files in the directory
    std::vector <std::string> files = get_files_in_dir(tag_array_index_dir);

    int number_of_file = files.size();

    std::cerr << "The list of files are: " << std::endl;
    for (auto &file: files) {
        std::cerr << file << std::endl;
    }


    // cerr << "Reading the whole genome r-index file" << endl;
    // FastLocate r_index;
    // if (!sdsl::load_from_file(r_index, r_index_file)) {
    //     std::cerr << "Cannot load the r-index from " << r_index_file << std::endl;
    //     std::exit(EXIT_FAILURE);
    // }


    vector<size_t> tag_count_per_file;
    tag_count_per_file.resize(number_of_file, 0);


    for (size_t i = 0; i < files.size(); i++) {
        sdsl::int_vector_buffer<8> in(files[i], std::ios::in);
        gbwt::size_type file_pos = 0;
        if (!in.is_open()) {
            throw std::runtime_error("Cannot open file: " + files[i]);
        }
        while (file_pos < in.size()) {
            auto tag_block = panindexer::TagArray::decode_run(
                    gbwt::ByteCode::read(in, file_pos)
            );
            tag_count_per_file[i]++;
            // tag_count_per_file[i] += tag_block.second;
        }


    }



    for (size_t i = 0; i < number_of_file; i++) {
        std::cerr << "File: " << files[i] << " Tags runs: " << tag_count_per_file[i] << std::endl;
    }
    exit(0);


//     std::cerr << "Finding the node to component mapping" << std::endl;
//     // finding the components
//     std::unordered_map<nid_t, size_t> node_to_comp_map = node_to_component(gbz);

//     std::vector<int> file_to_comp(number_of_file);
//     std::vector<int> comp_to_file(number_of_file);

//     std::cerr << "Initializing the reader" << std::endl;
//     FileReader reader(files, threads, 1000000);
//     std::cerr << "Creating the mapping from comp to tag files" << std::endl;
//     // for each tag block files, we read the first block and read the first node
//     for (auto i = 0; i < number_of_file; i++) {
//         // get the component of the node of the first block
//         size_t comp = node_to_comp_map[id(reader.get_first_tag(i))];
//         std::cerr << "The component of the first block of file " << files[i] << " is " << comp << " first tag " << reader.get_first_tag(i)  << " node id is " << id(reader.get_first_tag(i)) << std::endl;
//         file_to_comp[i] = comp;
//         comp_to_file[comp] = i;

//     }

//     std::cerr << "The mapping from comp to tag files is done" << std::endl;


//     auto total_strings = r_index.tot_strings();
//     std::vector<size_t> seq_id_to_comp_id;
//     seq_id_to_comp_id.resize(total_strings);



//     // get the first node of each path and get the component id of the node
// #pragma omp parallel for
//     for (size_t i = 0; i < total_strings; i++) {
//         auto seq_graph_nodes = gbz.index.extract(i * 2);
//         if (!seq_graph_nodes.empty()) {
//             size_t node_id = gbwt::Node::id(seq_graph_nodes[0]);
//             seq_id_to_comp_id[i] = node_to_comp_map.at(node_id);
//         }
//     }
//     std::cerr << "The mapping from seq id to comp id is done" << std::endl;


//     // note that the first #num_seq tags are correspond to the ENDMARKERs
//     auto total_tags_count = r_index.get_sequence_size() - total_strings;
//     std::cerr << "Total tags count " << total_tags_count << std::endl;



// //    auto total_strings = r_index.tot_strings();
//     auto first = r_index.locateFirst();

//     for (int i = 0; i < total_strings - 1 ; i++){
//         first = r_index.locateNext(first);
//     }

//     auto total_bwt_size = r_index.get_sequence_size();

//     vector <size_t> rindex_per_file;
//     rindex_per_file.resize(number_of_file, 0);

//     std::cerr << "Total string " << total_strings << " Total bwt size " << total_bwt_size << std::endl;

//     for (size_t i = total_strings; i < total_bwt_size; i++){
//         auto seq_id = r_index.seqId(first);
// //        std::cerr << "Seq id " << seq_id << std::endl;
//         // want to get the file number that is associated with the seq id
//         auto current_file = comp_to_file[seq_id_to_comp_id[seq_id]];
//         rindex_per_file[current_file] += 1;
//         first = r_index.locateNext(first);
//     }


//     for (size_t i = 0; i < number_of_file; i++) {
//         std::cerr << "File: " << files[i] << " R-index: " << rindex_per_file[i] << " Tag count from file " << tag_count_per_file[i] << std::endl;
//     }


}