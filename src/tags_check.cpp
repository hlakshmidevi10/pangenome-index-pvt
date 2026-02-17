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
#include <random>
#include <limits>
#include <tuple>
#include <vector>


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

// Print the sampled tag array completely: every run with run_id, BWT span, length, and tag (node_id, is_rev).
static int print_sampled_tag_array(const std::string& sampled_tags_file) {
    using namespace panindexer;
    std::cerr << "Loading sampled tag array: " << sampled_tags_file << std::endl;
    SampledTagArray arr;
    {
        std::ifstream sin(sampled_tags_file, std::ios::binary);
        if (!sin) {
            std::cerr << "Cannot open sampled tags file: " << sampled_tags_file << std::endl;
            return 1;
        }
        arr.load(sin);
    }
    arr.ensure_run_rank();
    arr.ensure_run_select();

    size_t n_runs = arr.total_runs();
    size_t bwt_size = (arr.run_starts().size() > 0) ? (arr.run_starts().size() - 1) : 0;
    std::cout << "# Sampled tag array: " << n_runs << " runs, BWT size " << bwt_size
              << " first_run_is_gap=" << arr.is_first_run_gap() << "\n";
    std::cout << "# run_id\tbwt_start\tbwt_end\tlen\ttag\tnode_id\tis_rev\n";

    for (size_t run_id = 0; run_id < n_runs; ++run_id) {
        auto [bwt_start, bwt_end] = arr.run_span(run_id);
        size_t len = (bwt_end >= bwt_start) ? (bwt_end - bwt_start + 1) : 0;
        uint64_t tag = arr.run_value(run_id);
        std::cout << run_id << "\t" << bwt_start << "\t" << bwt_end << "\t" << len << "\t" << tag;
        if (tag == 0) {
            std::cout << "\t-\t-\n";
        } else {
            auto [node_id, is_rev] = decode_tag(tag);
            std::cout << "\t" << node_id << "\t" << (is_rev ? 1 : 0) << "\n";
        }
    }
    std::cerr << "Printed " << n_runs << " runs." << std::endl;
    return 0;
}

// ======================= Verify non-encoded vs encoded r-index =======================
// Complete comparison: sizes, samples, last/last_to_run, LF, rankAt, run_id_and_offset_at,
// single-step LF, bwt_end_position_of_run, full decompressDA; reports every mismatch.
static int verify_rindex_encoded_vs_legacy(const std::string& legacy_ri_file,
                                           const std::string& encoded_ri_file,
                                           size_t random_trials) {
    using namespace panindexer;
    std::cerr << "Loading legacy (non-encoded) r-index: " << legacy_ri_file << std::endl;
    FastLocate legacy;
    {
        std::ifstream in(legacy_ri_file, std::ios::binary);
        if (!in) {
            std::cerr << "Cannot open legacy r-index: " << legacy_ri_file << std::endl;
            return 1;
        }
        legacy.load_encoded(in);
    }
    if (legacy.is_encoded()) {
        std::cerr << "ERROR: First file is encoded; expected non-encoded (legacy). Pass legacy .ri first, then encoded .ri." << std::endl;
        return 1;
    }

    std::cerr << "Loading encoded r-index: " << encoded_ri_file << std::endl;
    FastLocate encoded;
    {
        std::ifstream in(encoded_ri_file, std::ios::binary);
        if (!in) {
            std::cerr << "Cannot open encoded r-index: " << encoded_ri_file << std::endl;
            return 1;
        }
        encoded.load_encoded(in);
    }
    if (!encoded.is_encoded()) {
        std::cerr << "ERROR: Second file is not encoded; expected encoded .ri. Pass legacy .ri first, then encoded .ri." << std::endl;
        return 1;
    }

    int exit_code = 0;
    const size_t n = legacy.bwt_size();
    const size_t num_runs = legacy.tot_runs();
    std::vector<std::string> passed, failed;

    auto ok = [&](int c, const std::string& name, bool cond) {
        if (cond) { passed.push_back(name); std::cerr << "  Check " << c << " (" << name << ") passed." << std::endl; }
        else { failed.push_back(name); exit_code = 1; }
    };

    std::cerr << "--- Running checks ---" << std::endl;

    // 1) Global sizes
    bool c1 = (encoded.bwt_size() == n && encoded.tot_runs() == num_runs &&
               encoded.get_sequence_size() == legacy.get_sequence_size());
    if (!c1) {
        if (encoded.bwt_size() != n) std::cerr << "MISMATCH bwt_size(): legacy=" << n << " encoded=" << encoded.bwt_size() << std::endl;
        if (encoded.tot_runs() != num_runs) std::cerr << "MISMATCH tot_runs(): legacy=" << num_runs << " encoded=" << encoded.tot_runs() << std::endl;
        if (encoded.get_sequence_size() != legacy.get_sequence_size())
            std::cerr << "MISMATCH get_sequence_size(): legacy=" << legacy.get_sequence_size() << " encoded=" << encoded.get_sequence_size() << std::endl;
    }
    ok(1, "bwt_size, tot_runs, sequence_size", c1);

    // 2) C array (same size and values)
    bool c2 = (legacy.C.size() == encoded.C.size());
    if (c2) {
        for (size_t i = 0; i < legacy.C.size() && c2; ++i)
            if (legacy.C[i] != encoded.C[i]) {
                std::cerr << "MISMATCH C[" << i << "]: legacy=" << legacy.C[i] << " encoded=" << encoded.C[i] << std::endl;
                c2 = false;
            }
    } else {
        std::cerr << "MISMATCH C.size(): legacy=" << legacy.C.size() << " encoded=" << encoded.C.size() << std::endl;
    }
    ok(2, "C array", c2);

    // 3) sym_map (same)
    bool c3 = (legacy.sym_map.size() == encoded.sym_map.size());
    if (c3) {
        for (size_t i = 0; i < legacy.sym_map.size() && c3; ++i)
            if (legacy.sym_map[i] != encoded.sym_map[i]) {
                std::cerr << "MISMATCH sym_map[" << (int)i << "]: legacy=" << (int)legacy.sym_map[i] << " encoded=" << (int)encoded.sym_map[i] << std::endl;
                c3 = false;
            }
    } else {
        std::cerr << "MISMATCH sym_map.size(): legacy=" << legacy.sym_map.size() << " encoded=" << encoded.sym_map.size() << std::endl;
    }
    ok(3, "sym_map", c3);

    // 4) samples
    bool c4 = true;
    for (size_t run_id = 0; run_id < num_runs && c4; ++run_id) {
        if (legacy.getSample(run_id) != encoded.getSample(run_id)) {
            std::cerr << "MISMATCH getSample(run_id=" << run_id << "): legacy=" << legacy.getSample(run_id) << " encoded=" << encoded.getSample(run_id) << std::endl;
            c4 = false;
        }
    }
    ok(4, "samples", c4);

    // 5) last and last_to_run (check only a prefix and a tail to keep runtime low)
    const size_t last_check_prefix = 100000;
    const size_t last_check_tail = 5000;
    bool c5 = (legacy.last.size() == encoded.last.size() && legacy.last_to_run.size() == encoded.last_to_run.size());
    if (c5) {
        size_t n_last = legacy.last.size();
        size_t n_ltr = legacy.last_to_run.size();
        for (size_t i = 0; i < n_last && i < last_check_prefix && c5; ++i)
            if (legacy.last[i] != encoded.last[i]) {
                std::cerr << "MISMATCH last[" << i << "]" << std::endl;
                c5 = false;
            }
        for (size_t i = (n_last > last_check_prefix + last_check_tail ? n_last - last_check_tail : n_last); i < n_last && c5; ++i)
            if (legacy.last[i] != encoded.last[i]) {
                std::cerr << "MISMATCH last[" << i << "]" << std::endl;
                c5 = false;
            }
        for (size_t i = 0; i < n_ltr && i < last_check_prefix && c5; ++i)
            if (legacy.last_to_run[i] != encoded.last_to_run[i]) {
                std::cerr << "MISMATCH last_to_run[" << i << "]" << std::endl;
                c5 = false;
            }
        for (size_t i = (n_ltr > last_check_prefix + last_check_tail ? n_ltr - last_check_tail : n_ltr); i < n_ltr && c5; ++i)
            if (legacy.last_to_run[i] != encoded.last_to_run[i]) {
                std::cerr << "MISMATCH last_to_run[" << i << "]" << std::endl;
                c5 = false;
            }
    } else {
        std::cerr << "MISMATCH last.size() or last_to_run.size()" << std::endl;
    }
    ok(5, "last, last_to_run", c5);

    // 6) LF(range, sym) on random ranges and symbols
    std::mt19937 rng(12345);
    const unsigned char syms[] = {'A', 'C', 'G', 'T'};
    bool c6 = true;
    for (size_t t = 0; t < random_trials && n > 0 && c6; ++t) {
        size_t first = rng() % n;
        size_t second = rng() % n;
        if (first > second) std::swap(first, second);
        for (unsigned char sym : syms) {
            gbwt::range_type r_leg = legacy.LF({first, second}, sym);
            gbwt::range_type r_enc = encoded.LF_encoded({first, second}, sym);
            if (r_leg.first != r_enc.first || r_leg.second != r_enc.second) {
                std::cerr << "MISMATCH LF(range=[" << first << "," << second << "], sym='" << (char)sym << "')" << std::endl;
                c6 = false;
                break;
            }
        }
    }
    ok(6, "LF(range,sym)", c6);

    // 7) rankAt(pos, sym)
    bool c7 = true;
    for (size_t t = 0; t < random_trials && n > 0 && c7; ++t) {
        size_t pos = rng() % (n + 1);
        for (unsigned char sym : syms) {
            if (legacy.rankAt(pos, sym) != encoded.rankAt(pos, sym)) {
                std::cerr << "MISMATCH rankAt(pos=" << pos << ", sym='" << (char)sym << "')" << std::endl;
                c7 = false;
                break;
            }
        }
    }
    ok(7, "rankAt", c7);

    // 8) run_id_and_offset_at(pos)
    bool c8 = true;
    for (size_t t = 0; t < random_trials && n > 0 && c8; ++t) {
        size_t pos = rng() % n;
        size_t run_leg, run_enc, off_leg, off_enc;
        legacy.run_id_and_offset_at(pos, run_leg, off_leg);
        encoded.run_id_and_offset_at(pos, run_enc, off_enc);
        if (run_leg != run_enc || off_leg != off_enc) {
            std::cerr << "MISMATCH run_id_and_offset_at(pos=" << pos << ")" << std::endl;
            c8 = false;
        }
    }
    ok(8, "run_id_and_offset_at", c8);

    // 9) Single-step LF(pos)
    bool c9 = true;
    for (size_t t = 0; t < random_trials && n > 0 && c9; ++t) {
        size_t pos = rng() % n;
        if (legacy.LF(pos) != encoded.LF(pos)) {
            std::cerr << "MISMATCH LF(pos=" << pos << ")" << std::endl;
            c9 = false;
        }
    }
    ok(9, "LF(pos)", c9);

    // 10) bwt_end_position_of_run(run_id)
    bool c10 = true;
    for (size_t run_id = 0; run_id < num_runs && c10; ++run_id) {
        if (legacy.bwt_end_position_of_run(run_id) != encoded.bwt_end_position_of_run(run_id)) {
            std::cerr << "MISMATCH bwt_end_position_of_run(run_id=" << run_id << ")" << std::endl;
            c10 = false;
        }
    }
    ok(10, "bwt_end_position_of_run", c10);

    // 11) Full decompressDA()
    auto da_leg = legacy.decompressDA();
    auto da_enc = encoded.decompressDA();
    bool c11 = (da_leg.size() == da_enc.size());
    if (c11) {
        for (size_t i = 0; i < da_leg.size() && c11; ++i)
            if (da_leg[i] != da_enc[i]) {
                std::cerr << "MISMATCH decompressDA()[" << i << "]" << std::endl;
                c11 = false;
            }
    } else {
        std::cerr << "MISMATCH decompressDA().size(): legacy=" << da_leg.size() << " encoded=" << da_enc.size() << std::endl;
    }
    ok(11, "decompressDA", c11);

    // 12) Block runs: Blocks (legacy) vs EncodedBlock (encoded) — same runs per block
    size_t num_blocks_legacy = legacy.num_blocks();
    size_t num_blocks_encoded = encoded.num_blocks();
    bool c12 = (num_blocks_legacy == num_blocks_encoded);
    if (!c12) {
        std::cerr << "MISMATCH num_blocks(): legacy=" << num_blocks_legacy << " encoded=" << num_blocks_encoded << std::endl;
        failed.push_back("block_runs");
        exit_code = 1;
    } else {
        std::vector<std::pair<size_t, size_t>> runs_legacy, runs_encoded;
        size_t first_bad_block = static_cast<size_t>(-1);
        for (size_t block_id = 0; block_id < num_blocks_legacy; ++block_id) {
            legacy.get_block_runs(block_id, runs_legacy);
            encoded.get_block_runs(block_id, runs_encoded);
            if (runs_legacy.size() != runs_encoded.size()) {
                c12 = false;
                std::cerr << "MISMATCH block " << block_id << " run count: legacy=" << runs_legacy.size()
                          << " encoded=" << runs_encoded.size() << std::endl;
                // Print actual runs for this block (sym as byte; show char if printable)
                std::cerr << "  Legacy runs (" << runs_legacy.size() << "):" << std::endl;
                for (size_t r = 0; r < runs_legacy.size(); ++r) {
                    size_t sym = runs_legacy[r].first;
                    size_t len = runs_legacy[r].second;
                    char sym_char = (sym < 128 && sym >= 32) ? static_cast<char>(sym) : '?';
                    std::cerr << "    run " << r << ": sym=" << sym << " ('" << sym_char << "') len=" << len << std::endl;
                }
                std::cerr << "  Encoded runs (" << runs_encoded.size() << "):" << std::endl;
                for (size_t r = 0; r < runs_encoded.size(); ++r) {
                    size_t sym = runs_encoded[r].first;
                    size_t len = runs_encoded[r].second;
                    char sym_char = (sym < 128 && sym >= 32) ? static_cast<char>(sym) : '?';
                    std::cerr << "    run " << r << ": sym=" << sym << " ('" << sym_char << "') len=" << len << std::endl;
                }
                // Extra context only for first block run-count mismatch
                if (first_bad_block == static_cast<size_t>(-1)) {
                    first_bad_block = block_id;
                    std::cerr << "  Context: num_blocks(legacy)=" << num_blocks_legacy
                              << " num_blocks(encoded)=" << num_blocks_encoded
                              << " bwt_size(legacy)=" << legacy.bwt_size()
                              << " bwt_size(encoded)=" << encoded.bwt_size()
                              << " tot_runs(legacy)=" << legacy.tot_runs()
                              << " tot_runs(encoded)=" << encoded.tot_runs() << std::endl;
                    std::cerr << "  Legacy block_size=" << legacy.block_size
                              << " encoded_block_size(encoded)=" << encoded.encoded_block_size << std::endl;
                    size_t lo = (block_id >= 3) ? block_id - 3 : 0;
                    size_t hi = (block_id + 4 <= num_blocks_legacy) ? block_id + 4 : num_blocks_legacy;
                    std::cerr << "  Run counts per block (block_id -> legacy, encoded) for blocks [" << lo << ".." << (hi-1) << "]:" << std::endl;
                    for (size_t b = lo; b < hi; ++b) {
                        legacy.get_block_runs(b, runs_legacy);
                        encoded.get_block_runs(b, runs_encoded);
                        std::cerr << "    block " << b << ": legacy=" << runs_legacy.size() << " encoded=" << runs_encoded.size()
                                  << (runs_legacy.size() != runs_encoded.size() ? "  <-- MISMATCH" : "") << std::endl;
                    }
                }
                exit_code = 1;
            } else {
                for (size_t r = 0; r < runs_legacy.size(); ++r) {
                    if (runs_legacy[r].first != runs_encoded[r].first || runs_legacy[r].second != runs_encoded[r].second) {
                        c12 = false;
                        std::cerr << "MISMATCH block " << block_id << " run " << r << ": legacy (sym=" << runs_legacy[r].first
                                  << ", len=" << runs_legacy[r].second << ") encoded (sym=" << runs_encoded[r].first
                                  << ", len=" << runs_encoded[r].second << ")" << std::endl;
                        exit_code = 1;
                    }
                }
            }
        }
        if (c12)
            std::cerr << "  Check 12 (block_runs) passed." << std::endl;
        else
            failed.push_back("block_runs");
        if (c12) passed.push_back("block_runs");
    }

    // Summary
    std::cerr << "--- Summary ---" << std::endl;
    std::cerr << "Passed (" << passed.size() << "): ";
    for (size_t i = 0; i < passed.size(); ++i) std::cerr << (i ? ", " : "") << passed[i];
    std::cerr << std::endl;
    if (!failed.empty()) {
        std::cerr << "Failed (" << failed.size() << "): ";
        for (size_t i = 0; i < failed.size(); ++i) std::cerr << (i ? ", " : "") << failed[i];
        std::cerr << std::endl;
        std::cerr << "One or more mismatches found; see above for details." << std::endl;
    } else {
        std::cerr << "All checks passed: legacy and encoded r-index agree (bwt_size=" << n << ", tot_runs=" << num_runs
                  << ", random trials=" << random_trials << ")." << std::endl;
    }
    return exit_code;
}

// Encode pos_t to the same integer code used by SampledTagArray (gaps -> 0).
static uint64_t encode_pos_for_sampled(handlegraph::pos_t p) {
    if (gbwtgraph::offset(p) != 0) return 0;
    if (gbwtgraph::id(p) == 0) return 0;
    return panindexer::SampledTagArray::encode_value(gbwtgraph::id(p), gbwtgraph::is_rev(p));
}

// Verify sampled tag array against the full (compressed) tag array by traversing from the beginning
// run-by-run and comparing the tag value at each run start.
static int verify_sampled_vs_tag_array(const std::string& compressed_tags_file,
                                       const std::string& sampled_tags_file) {
    using namespace panindexer;
    std::cerr << "Loading compressed tag array: " << compressed_tags_file << std::endl;
    TagArray tag_array;
    {
        std::ifstream tin(compressed_tags_file, std::ios::binary);
        if (!tin) {
            std::cerr << "Cannot open compressed tags file: " << compressed_tags_file << std::endl;
            return 1;
        }
        tag_array.load_compressed_tags_compact(tin);
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

    size_t tag_bwt_size = tag_array.bwt_size();
    size_t sampled_bwt_size = (sampled.run_starts().size() > 0) ? (sampled.run_starts().size() - 1) : 0;
    if (tag_bwt_size != sampled_bwt_size) {
        std::cerr << "WARNING: Tag array BWT size (" << tag_bwt_size
                  << ") != sampled tag array BWT size (" << sampled_bwt_size
                  << "). Verification may be wrong." << std::endl;
    }

    sampled.ensure_run_rank();
    sampled.ensure_run_select();

    size_t run_index = 0;
    size_t mismatches = 0;
    const size_t max_mismatch_print = 20;

    tag_array.for_each_run_compact_with_bwt([&](handlegraph::pos_t p, uint64_t len,
                                                 size_t bwt_start, size_t bwt_end) {
        uint64_t expected = encode_pos_for_sampled(p);
        if (bwt_start >= sampled.run_starts().size()) {
            std::cerr << "Run " << run_index << ": bwt_start=" << bwt_start
                      << " >= sampled size (" << sampled.run_starts().size() << "), skipping." << std::endl;
            run_index++;
            return;
        }
        size_t run_id = sampled.run_id_at(bwt_start);
        uint64_t got = sampled.run_value(run_id);
        if (expected != got) {
            mismatches++;
            if (mismatches <= max_mismatch_print) {
                std::cerr << "MISMATCH run " << run_index
                          << " bwt_start=" << bwt_start << " bwt_end=" << bwt_end
                          << " tag_array: node_id=" << gbwtgraph::id(p)
                          << " is_rev=" << gbwtgraph::is_rev(p)
                          << " offset=" << gbwtgraph::offset(p)
                          << " len=" << len << " expected_tag=" << expected
                          << " sampled: run_id=" << run_id << " got_tag=" << got;
                if (got != 0) {
                    auto [got_node_id, got_is_rev] = decode_tag(got);
                    std::cerr << " sampled_decoded: node_id=" << got_node_id << " is_rev=" << got_is_rev;
                } else {
                    std::cerr << " sampled_decoded: gap";
                }
                std::cerr << std::endl;
            }
        }
        run_index++;
    });

    std::cerr << "Verification done: " << run_index << " runs checked, "
              << mismatches << " mismatch(es)." << std::endl;
    if (mismatches > max_mismatch_print) {
        std::cerr << "  (first " << max_mismatch_print << " mismatches printed)" << std::endl;
    }
    return mismatches > 0 ? 1 : 0;
}

// Verify sampled tag array by traversing a GBWT path from the start and comparing
// the tag at each (seq_id, base_offset) with the node on the path at that offset.
// Uses only .gbz (GBZ contains both the GBWT index and the graph).
static int verify_sampled_against_gbwt(const std::string& r_index_file,
                                       const std::string& sampled_tags_file,
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
    std::cerr << "Encoded r-index BWT character counts:" << std::endl;
    r_index.print_character_counts(std::cerr);
    std::cerr << "Loading GBZ (GBWT index + graph): " << gbz_file << std::endl;
    gbwtgraph::GBZ gbz;
    sdsl::simple_sds::load_from(gbz, gbz_file);
    const gbwt::GBWT& gbwt_index = gbz.index;
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
    const size_t max_mismatch_print = 30;
    for (gbwt::node_type node : path) {
        if (node == gbwt::ENDMARKER) break;
        int64_t node_id = gbwt::Node::id(node);
        bool node_rev = gbwt::Node::is_reverse(node);
        gbwtgraph::handle_t handle = graph.get_handle(node_id, node_rev);
        size_t node_len = graph.get_length(handle);
        uint64_t expected_tag = SampledTagArray::encode_value(node_id, node_rev);

        // Check only at the start of each node (offset 0 within node)
        size_t base_offset = cumulative_bases;
        size_t text_pos = r_index.pack(seq_id, base_offset);
        size_t bwt_pos = 0;
        if (text_pos_to_bwt_pos(r_index, text_pos, seq_end_text_pos, bwt_pos)) {
            size_t run_id = sampled.run_id_at(bwt_pos);
            uint64_t tag_val = sampled.run_value(run_id);
            checked++;
            bool match = (tag_val == expected_tag);
            if (!match) {
                mismatches++;
                if (mismatches <= max_mismatch_print) {
                    size_t packed_at = recover_text_pos_from_bwt(r_index, bwt_pos);
                    auto [seq_id_at, offset_at] = r_index.unpack(packed_at);
                    std::cerr << "MISMATCH node_index=" << node_index << " base_offset=" << base_offset
                              << " bwt_pos=" << bwt_pos
                              << " text: packed=" << packed_at << " seq_id=" << seq_id_at << " offset=" << offset_at;
                    if (seq_id_at != path_id || offset_at != base_offset) {
                        std::cerr << " [expected seq_id=" << path_id << " offset=" << base_offset << "]";
                    }
                    std::cerr << "  expected: node_id=" << node_id << " rev=" << node_rev << " tag=" << expected_tag;
                    if (tag_val == 0) {
                        std::cerr << "  sampled: gap (tag=0)";
                    } else {
                        auto [got_id, got_rev] = decode_tag(tag_val);
                        std::cerr << "  sampled: tag=" << tag_val << " node_id=" << got_id << " rev=" << got_rev;
                    }
                    std::cerr << std::endl;
                }
            }
        }
        cumulative_bases += node_len;
        node_index++;
    }
    std::cerr << "Verification done: " << checked << " node starts checked, "
              << mismatches << " mismatch(es)." << std::endl;
    if (mismatches > max_mismatch_print) {
        std::cerr << "  (first " << max_mismatch_print << " mismatches printed above)" << std::endl;
    }
    return mismatches > 0 ? 1 : 0;
}

// Verify full compact tag array against GBWT: for each run, get expected (node_id, is_rev, offset)
// from the GBWT path at the text position corresponding to the run start, and compare with the tag.
// Uses only .gbz (GBZ contains both the GBWT index and the graph).
static int verify_compact_tags_against_gbwt(const std::string& r_index_file,
                                            const std::string& compact_tags_file,
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
    std::cerr << "Loading compact tag array: " << compact_tags_file << std::endl;
    TagArray tag_array;
    {
        std::ifstream tin(compact_tags_file, std::ios::binary);
        if (!tin) {
            std::cerr << "Cannot open compact tags file: " << compact_tags_file << std::endl;
            return 1;
        }
        tag_array.load_compressed_tags_compact(tin);
    }
    std::cerr << "Loading GBZ (contains GBWT index + graph): " << gbz_file << std::endl;
    gbwtgraph::GBZ gbz;
    sdsl::simple_sds::load_from(gbz, gbz_file);
    const gbwtgraph::GBWTGraph& graph = gbz.graph;
    const gbwt::GBWT& gbwt_index = gbz.index;

    gbwt::vector_type path = gbwt_index.extract(gbwt::Path::encode(path_id, false));
    if (path.empty()) {
        std::cerr << "Path " << path_id << " not found in GBWT." << std::endl;
        return 1;
    }
    // Build prefix sum of node lengths so we can map base_offset -> (node_index, off_in_node)
    std::vector<size_t> node_prefix;
    node_prefix.push_back(0);
    size_t path_length_bases = 0;
    for (gbwt::node_type node : path) {
        if (node == gbwt::ENDMARKER) break;
        path_length_bases += graph.get_length(
            graph.get_handle(gbwt::Node::id(node), gbwt::Node::is_reverse(node)));
        node_prefix.push_back(path_length_bases);
    }
    std::cerr << "Path " << path_id << " has " << (node_prefix.size() - 1) << " nodes, "
              << path_length_bases << " bases." << std::endl;

    auto offset_to_pos = [&](size_t base_offset) -> handlegraph::pos_t {
        if (base_offset >= path_length_bases) {
            return make_pos_t(0, false, 0);
        }
        size_t i = 0;
        while (i + 1 < node_prefix.size() && node_prefix[i + 1] <= base_offset) ++i;
        gbwt::node_type node = path[i];
        if (node == gbwt::ENDMARKER) return make_pos_t(0, false, 0);
        int64_t node_id = gbwt::Node::id(node);
        bool node_rev = gbwt::Node::is_reverse(node);
        size_t off_in_node = base_offset - node_prefix[i];
        return make_pos_t(node_id, node_rev, static_cast<size_t>(off_in_node));
    };

    size_t run_index = 0;
    size_t checked = 0;
    size_t mismatches = 0;
    size_t skipped_other_seq = 0;
    const size_t max_mismatch_print = 30;

    tag_array.for_each_run_compact_with_bwt([&](handlegraph::pos_t p, uint64_t len,
                                                  size_t bwt_start, size_t bwt_end) {
        if (run_index % sample_every != 0) {
            run_index++;
            return;
        }
        size_t text_pos = recover_text_pos_from_bwt(r_index, bwt_start);
        auto [seq_id, offset] = r_index.unpack(text_pos);
        if (seq_id != path_id) {
            skipped_other_seq++;
            run_index++;
            return;
        }
        if (offset >= path_length_bases) {
            run_index++;
            return;
        }
        handlegraph::pos_t expected = offset_to_pos(offset);
        bool match = (gbwtgraph::id(p) == gbwtgraph::id(expected) &&
                      gbwtgraph::is_rev(p) == gbwtgraph::is_rev(expected) &&
                      gbwtgraph::offset(p) == gbwtgraph::offset(expected));
        checked++;
        if (!match) mismatches++;

        if (!match && mismatches <= max_mismatch_print) {
            std::cerr << "MISMATCH run " << run_index << " bwt_start=" << bwt_start << " bwt_end=" << bwt_end
                      << " text: seq_id=" << seq_id << " offset=" << offset
                      << " compact: node_id=" << gbwtgraph::id(p) << " is_rev=" << gbwtgraph::is_rev(p)
                      << " offset=" << gbwtgraph::offset(p)
                      << " expected: node_id=" << gbwtgraph::id(expected) << " is_rev=" << gbwtgraph::is_rev(expected)
                      << " offset=" << gbwtgraph::offset(expected) << std::endl;
        }
        run_index++;
    });

    std::cerr << "Verification done: " << checked << " runs checked (path_id=" << path_id << "), "
              << skipped_other_seq << " runs skipped (other seq), "
              << mismatches << " mismatch(es)." << std::endl;
    if (mismatches > max_mismatch_print) {
        std::cerr << "  (first " << max_mismatch_print << " mismatches printed)" << std::endl;
    }
    return mismatches > 0 ? 1 : 0;
}

static void usage_verify(const char* prog) {
    std::cerr << "Usage: " << prog << " --verify-sampled <r_index.ri> <sampled.tags> --gbz <graph.gbz> [options]\n"
              << "  Verify that the sampled tag array matches the GBWT path at node starts only.\n"
              << "  Only positions at the start of each node are checked; only incorrect ones are printed.\n"
              << "  Uses only .gbz (GBZ contains both the GBWT index and the graph).\n"
              << "Options:\n"
              << "  --gbz FILE        GBZ file (required; contains GBWT + graph)\n"
              << "  --path-id N       GBWT path/sequence ID to traverse (default: 0)\n" << std::endl;
}

static void usage_verify_compact_vs_gbwt(const char* prog) {
    std::cerr << "Usage: " << prog << " --verify-compact-vs-gbwt <r_index.ri> <compact_tags.tags> --gbz <graph.gbz> [options]\n"
              << "  Verify that the full compact tag array matches the GBWT path.\n"
              << "  Uses only .gbz (GBZ contains both the GBWT index and the graph).\n"
              << "  For each run, locates the run start in text space and compares (node_id, is_rev, offset) with the path.\n"
              << "Options:\n"
              << "  --gbz FILE        GBZ file (required; contains GBWT + graph)\n"
              << "  --path-id N       GBWT path/sequence ID to compare (default: 0)\n"
              << "  --sample-every N  Check every N-th run only (default: 1 = every run)\n" << std::endl;
}

static void usage_verify_sampled_vs_tags(const char* prog) {
    std::cerr << "Usage: " << prog << " --verify-sampled-vs-tags <compressed_tags.tags> <sampled.tags>\n"
              << "  Verify that the sampled tag array was built correctly from the full (compressed) tag array.\n"
              << "  Traverses from the beginning run-by-run and compares the tag value at each run start.\n"
              << std::endl;
}

static void usage_verify_rindex_encoded_vs_legacy(const char* prog) {
    std::cerr << "Usage: " << prog << " --verify-rindex-encoded-vs-legacy <legacy.ri> <encoded.ri> [--trials N]\n"
              << "  Complete comparison of non-encoded (legacy) and encoded r-index.\n"
              << "  Checks: bwt_size, tot_runs, C, sym_map, samples, last, last_to_run,\n"
              << "  LF(range,sym), rankAt, run_id_and_offset_at, LF(pos), bwt_end_position_of_run, decompressDA,\n"
              << "  and block runs (Blocks vs EncodedBlock: same runs per block).\n"
              << "  --trials N  number of random positions/ranges for LF/rank/run_id checks (default: 500)\n"
              << std::endl;
}

// Print statistics for compact tag array: total runs, BWT size, first few tags, and node_id range
static int print_tags_stats(const std::string& compact_tags_file, size_t num_first_tags = 10) {
    using namespace panindexer;
    std::cerr << "Loading compact tag array: " << compact_tags_file << std::endl;
    TagArray tag_array;
    {
        std::ifstream tin(compact_tags_file, std::ios::binary);
        if (!tin) {
            std::cerr << "Cannot open compact tags file: " << compact_tags_file << std::endl;
            return 1;
        }
        tag_array.load_compressed_tags_compact(tin);
    }

    size_t bwt_size = tag_array.bwt_size();
    size_t total_runs = 0;
    size_t total_tag_length = 0;
    size_t gap_runs = 0;
    size_t gap_length = 0;
    
    int64_t min_node_id = std::numeric_limits<int64_t>::max();
    int64_t max_node_id = std::numeric_limits<int64_t>::min();
    bool has_valid_node = false;

    std::vector<std::tuple<handlegraph::pos_t, uint64_t, size_t, size_t>> first_tags;

    tag_array.for_each_run_compact_with_bwt([&](handlegraph::pos_t p, uint64_t len,
                                                 size_t bwt_start, size_t bwt_end) {
        total_runs++;
        total_tag_length += len;
        
        int64_t node_id = gbwtgraph::id(p);
        bool is_rev = gbwtgraph::is_rev(p);
        size_t offset = gbwtgraph::offset(p);
        
        // Track node_id range for all tags with valid node_id
        if (node_id > 0) {
            has_valid_node = true;
            if (node_id < min_node_id) min_node_id = node_id;
            if (node_id > max_node_id) max_node_id = node_id;
        }
        
        if (node_id == 0 || offset != 0) {
            // Gap or invalid tag
            gap_runs++;
            gap_length += len;
        } else {
            // Valid tag (node_id > 0 and offset == 0)
            if (first_tags.size() < num_first_tags) {
                first_tags.push_back({p, len, bwt_start, bwt_end});
            }
        }
    });

    // Print statistics
    std::cout << "# Compact Tag Array Statistics\n";
    std::cout << "# File: " << compact_tags_file << "\n";
    std::cout << "# BWT size: " << bwt_size << "\n";
    std::cout << "# Total runs: " << total_runs << "\n";
    std::cout << "# Total tag length: " << total_tag_length << "\n";
    std::cout << "# Gap runs: " << gap_runs << " (length: " << gap_length << ")\n";
    std::cout << "# Non-gap runs: " << (total_runs - gap_runs) << " (length: " << (total_tag_length - gap_length) << ")\n";
    
    if (has_valid_node) {
        std::cout << "# Node ID range: [" << min_node_id << ", " << max_node_id << "]\n";
    } else {
        std::cout << "# Node ID range: [no valid nodes found]\n";
    }
    
    std::cout << "\n# First " << first_tags.size() << " tags:\n";
    std::cout << "# run_index\tbwt_start\tbwt_end\tlen\tnode_id\tis_rev\toffset\n";
    
    for (size_t i = 0; i < first_tags.size(); ++i) {
        auto [p, len, bwt_start, bwt_end] = first_tags[i];
        int64_t node_id = gbwtgraph::id(p);
        bool is_rev = gbwtgraph::is_rev(p);
        size_t offset = gbwtgraph::offset(p);
        
        std::cout << i << "\t" << bwt_start << "\t" << bwt_end << "\t" << len << "\t";
        if (node_id == 0 || offset != 0) {
            std::cout << "-\t-\t-\n";
        } else {
            std::cout << node_id << "\t" << (is_rev ? 1 : 0) << "\t" << offset << "\n";
        }
    }
    
    std::cerr << "Statistics printed successfully." << std::endl;
    return 0;
}

int main(int argc, char **argv) {
    if (argc >= 2 && std::strcmp(argv[1], "--tags-stats") == 0) {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " --tags-stats <compact_tags.tags> [--num-first N]\n"
                      << "  Print statistics for compact tag array:\n"
                      << "    - Total runs, BWT size, gap statistics\n"
                      << "    - First N tags (default: 10)\n"
                      << "    - Node ID range (min, max)\n"
                      << "Options:\n"
                      << "  --num-first N  Number of first tags to print (default: 10)\n";
            return 1;
        }
        std::string compact_tags_file = argv[2];
        size_t num_first = 10;
        for (int i = 3; i < argc; ++i) {
            if (std::strcmp(argv[i], "--num-first") == 0 && i + 1 < argc) {
                num_first = static_cast<size_t>(std::stoull(argv[++i]));
            }
        }
        return print_tags_stats(compact_tags_file, num_first);
    }

    if (argc >= 3 && std::strcmp(argv[1], "--print-sampled-tags") == 0) {
        std::string sampled_tags_file = argv[2];
        return print_sampled_tag_array(sampled_tags_file);
    }

    if (argc >= 4 && std::strcmp(argv[1], "--verify-sampled-vs-tags") == 0) {
        std::string compressed_tags_file = argv[2];
        std::string sampled_tags_file = argv[3];
        return verify_sampled_vs_tag_array(compressed_tags_file, sampled_tags_file);
    }

    if (argc >= 2 && std::strcmp(argv[1], "--verify-rindex-encoded-vs-legacy") == 0) {
        std::string legacy_ri, encoded_ri;
        size_t trials = 500;
        for (int i = 2; i < argc; ++i) {
            if (std::strcmp(argv[i], "--trials") == 0 && i + 1 < argc) {
                trials = static_cast<size_t>(std::stoull(argv[++i]));
                if (trials == 0) trials = 500;
            } else if (argv[i][0] != '-') {
                if (legacy_ri.empty()) legacy_ri = argv[i];
                else if (encoded_ri.empty()) encoded_ri = argv[i];
            }
        }
        if (legacy_ri.empty() || encoded_ri.empty()) {
            usage_verify_rindex_encoded_vs_legacy(argv[0]);
            return 1;
        }
        return verify_rindex_encoded_vs_legacy(legacy_ri, encoded_ri, trials);
    }

    if (argc >= 2 && std::strcmp(argv[1], "--verify-sampled") == 0) {
        std::string r_index_file;
        std::string sampled_tags_file;
        std::string gbz_file;
        size_t path_id = 0;
        size_t sample_every = 1;
        for (int i = 2; i < argc; i++) {
            if (std::strcmp(argv[i], "--gbz") == 0 && i + 1 < argc) {
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
        if (r_index_file.empty() || sampled_tags_file.empty() || gbz_file.empty()) {
            usage_verify(argv[0]);
            return 1;
        }
        return verify_sampled_against_gbwt(r_index_file, sampled_tags_file,
                                          gbz_file, path_id, sample_every);
    }

    if (argc >= 2 && std::strcmp(argv[1], "--verify-compact-vs-gbwt") == 0) {
        std::string r_index_file;
        std::string compact_tags_file;
        std::string gbz_file;
        size_t path_id = 0;
        size_t sample_every = 1;
        for (int i = 2; i < argc; i++) {
            if (std::strcmp(argv[i], "--gbz") == 0 && i + 1 < argc) {
                gbz_file = argv[++i];
            } else if (std::strcmp(argv[i], "--path-id") == 0 && i + 1 < argc) {
                path_id = static_cast<size_t>(std::stoull(argv[++i]));
            } else if (std::strcmp(argv[i], "--sample-every") == 0 && i + 1 < argc) {
                sample_every = static_cast<size_t>(std::stoull(argv[++i]));
            } else if (argv[i][0] != '-') {
                if (r_index_file.empty()) r_index_file = argv[i];
                else if (compact_tags_file.empty()) compact_tags_file = argv[i];
            }
        }
        if (r_index_file.empty() || compact_tags_file.empty() || gbz_file.empty()) {
            usage_verify_compact_vs_gbwt(argv[0]);
            return 1;
        }
        return verify_compact_tags_against_gbwt(r_index_file, compact_tags_file,
                                                gbz_file, path_id, sample_every);
    }

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <gbz_graph> <r_index.ri> <tag_array_index_dir>\n"
                  << "   Or: " << argv[0] << " --tags-stats <compact_tags.tags> [--num-first N]\n"
                  << "   Or: " << argv[0] << " --print-sampled-tags <sampled.tags>\n"
                  << "   Or: " << argv[0] << " --verify-sampled-vs-tags <compressed_tags.tags> <sampled.tags>\n"
                  << "   Or: " << argv[0] << " --verify-rindex-encoded-vs-legacy <legacy.ri> <encoded.ri> [--trials N]\n"
                  << "   Or: " << argv[0] << " --verify-sampled <r_index.ri> <sampled.tags> --gbz <graph.gbz> [--path-id N] [--sample-every N]\n"
                  << "   Or: " << argv[0] << " --verify-compact-vs-gbwt <r_index.ri> <compact_tags.tags> --gbz <graph.gbz> [--path-id N] [--sample-every N]\n";
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