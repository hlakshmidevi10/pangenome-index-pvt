//
// Created by Harieasswar Lakshmidevi on 8/30/25.
//

#include "pangenome_index/tag_arrays.hpp"
#include "pangenome_index/r-index.hpp"
#include <iostream>
#include <filesystem>

using namespace std;
using namespace panindexer;

void compress_tags_file(const std::string& r_index_file, const std::string& input_tags_file, const std::string& output_prefix) {
    std::cerr << "Loading tags from file: " << input_tags_file << std::endl;

    // Open input file for reading
    sdsl::int_vector_buffer<8> in(input_tags_file, std::ios::in);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open input file: " + input_tags_file);
    }
    std::cerr << "DEBUG: Input file opened, size: " << in.size() << " bytes" << std::endl;

    cerr << "Reading the whole genome r-index file" << endl;
    FastLocate r_index;
    {
        std::ifstream rin(r_index_file, std::ios::binary);
        if (!rin) {
            std::cerr << "Cannot open the r-index file " << r_index_file << std::endl;
            std::exit(EXIT_FAILURE);
        }
        r_index.load_encoded(rin);
    }


    // Create output file names
    const std::string compressed_tags_file = output_prefix + "_compressed.tags";
    const std::string encoded_starts_file = output_prefix + "_encoded_starts.bin";
    const std::string bwt_intervals_file = output_prefix + "_bwt_intervals.bin";

    // Remove existing files if they exist
    if (std::filesystem::exists(compressed_tags_file)) {
        std::remove(compressed_tags_file.c_str());
        std::cerr << "Existing compressed tags file deleted.\n";
    }

    if (std::filesystem::exists(encoded_starts_file)) {
        std::remove(encoded_starts_file.c_str());
        std::cerr << "Existing encoded starts file deleted.\n";
    }

    if (std::filesystem::exists(bwt_intervals_file)) {
        std::remove(bwt_intervals_file.c_str());
        std::cerr << "Existing bwt intervals file deleted.\n";
    }

    // Open output files
    std::ofstream out(compressed_tags_file, std::ios::binary);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot open compressed tags file for writing.\n";
        return;
    }

    std::ofstream out_encoded_starts(encoded_starts_file, std::ios::binary);
    if (!out_encoded_starts.is_open()) {
        std::cerr << "Error: Cannot open encoded starts file for writing.\n";
        return;
    }

    std::ofstream out_bwt_intervals(bwt_intervals_file, std::ios::binary);
    if (!out_bwt_intervals.is_open()) {
        std::cerr << "Error: Cannot open bwt intervals file for writing.\n";
        return;
    }
    std::cerr << "DEBUG: Output files opened successfully" << std::endl;

    // Write placeholder for size
    size_t placeholder = 0;
    out.write(reinterpret_cast<const char*>(&placeholder), sizeof(size_t));

    TagArray tag_array;
    gbwt::size_type file_position = 0;
    size_t batch_size = 100000; // Process tags in batches
    size_t total_runs_processed = 0;
    size_t total_tags_processed = 0;


    std::vector<std::pair<pos_t, uint16_t>> temp_tag_runs;
    auto num_endmarkers = r_index.tot_strings();
    temp_tag_runs.push_back(std::make_pair(pos_t{0, 1, 0}, num_endmarkers));
    total_tags_processed += num_endmarkers;
    total_runs_processed += 1;
    tag_array.compressed_serialize(out, out_encoded_starts, out_bwt_intervals, temp_tag_runs);
    //    tag_array.serialize_run_by_run(out, endmarkers);


    cerr << "Num Endmarkers: " << num_endmarkers << endl;
    size_t bwt_index = 0;

    std::cerr << "Processing tags in batches of " << batch_size << std::endl;

    while (file_position < in.size()) {
        std::vector<std::pair<pos_t, uint16_t>> tag_batch;
        tag_batch.reserve(batch_size);

        // Read a batch of tags
        size_t batch_count = 0;
        while (batch_count < batch_size && file_position < in.size()) {
            auto tag_block = panindexer::TagArray::decode_run(
                gbwt::ByteCode::read(in, file_position)
            );

            // // Print decoded position and iteration count
            // cerr << "BWT index: " << bwt_index << ": "
            //      << "pos_t{node=" << gbwtgraph::id(tag_block.first)
            //      << ", offset=" << gbwtgraph::offset(tag_block.first)
            //      << ", is_rev=" << gbwtgraph::is_rev(tag_block.first)
            //      << "}, run_len=" << tag_block.second << endl;

            bwt_index += tag_block.second;

            tag_batch.push_back(tag_block);
            batch_count++;
        }

        if (!tag_batch.empty()) {
            // Compress and serialize the batch
            tag_array.compressed_serialize(out, out_encoded_starts, out_bwt_intervals, tag_batch);
            total_runs_processed += tag_batch.size();
            std:cerr << "Tag Batch Size processed: " << tag_batch.size() << std::endl;

            std::cerr << "Iterating through runs_to_add vector (size: " << tag_batch.size() << ")" << std::endl;
            for (size_t i = 0; i < tag_batch.size(); ++i) {
                const auto& run = tag_batch[i];
                // std::cerr << "BWT Start Pos: " << total_tags_processed << " Run " << i << ": Position=" << run.first
                //           << ", Length=" << run.second << std::endl;
                total_tags_processed += run.second;
            }

            if (total_runs_processed % 1000000 == 0) {
                std::cerr << "Processed " << total_runs_processed << " tag runs..." << std::endl;
            }
        }
    }
    std::cerr << "DEBUG: Batch processing complete" << std::endl;

    // Close individual files
    out.close();
    out_encoded_starts.close();
    out_bwt_intervals.close();

    std::cerr << "Total tag runs processed: " << total_runs_processed << std::endl;
    std::cerr << "Merging compressed files..." << std::endl;

    // Merge the compressed files into final format
    tag_array.merge_compressed_files(compressed_tags_file, encoded_starts_file, bwt_intervals_file);

    std::cerr << "Compression complete! Output files:" << std::endl;
    std::cerr << "  - " << compressed_tags_file << std::endl;
    std::cerr << "Successfully compressed tags file." << std::endl;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <r_index_file> <input_tags_file> <output_prefix>" << std::endl;
        std::cerr << "Example: " << argv[0] << "input.ri input.tags output" << std::endl;
        return 1;
    }

    std::string r_index_file = argv[1];
    std::string input_tags_file = argv[2];
    std::string output_prefix = argv[3];

    std::cerr << "DEBUG: Starting compression of " << input_tags_file << std::endl;

    try {
        compress_tags_file(r_index_file, input_tags_file, output_prefix);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    std::cerr << "DEBUG: Compression completed successfully" << std::endl;
    return 0;
}