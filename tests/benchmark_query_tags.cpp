//
// Created by seeskand on 9/18/24.
//
// Benchmark test for query_tags.cpp performance
//

#include "../include/pangenome_index/r-index.hpp"
#include "../include/pangenome_index/tag_arrays.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>

using namespace std;
using namespace panindexer;
using namespace std::chrono;

#ifndef TIME
#define TIME 1
#endif

std::vector<std::string> readSequencesFromFile(const std::string &filename) {
    std::ifstream file(filename);
    std::vector<std::string> sequences;
    std::string line;
    while (std::getline(file, line)) {
        if (!line.empty()) {
            sequences.push_back(line);
        }
    }
    return sequences;
}

void sampleKmersFromFile(const std::string &filename, int k, int sampleCount, std::vector<std::string> &kmers) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Error opening file: " + filename);
    }
    std::default_random_engine generator(static_cast<long unsigned int>(std::time(0)));
    std::string line;
    int totalKmersSampled = 0;
    int lineIndex = 0;
    while (std::getline(file, line)) {
        if (line.empty() || line.size() < k) continue; // Skip empty or short sequences
        std::uniform_int_distribution<int> positionDist(0, line.size() - k);
        int kmersToSample = 50;
        for (int i = 0; i < kmersToSample && totalKmersSampled < sampleCount; ++i) {
            int pos = positionDist(generator);
            kmers.push_back(line.substr(pos, k));
            totalKmersSampled++;
        }
        if (totalKmersSampled >= sampleCount) break; // Stop if enough k-mers are collected
        lineIndex++;
    }
}

int main(int argc, char **argv) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <r_index.ri> <compressed_tags.tags> <k> <sequence_file>" << std::endl;
        return 1;
    }

    std::string r_index_file = std::string(argv[1]);
    std::string tag_array_index = std::string(argv[2]);
    int k = std::stoi(argv[3]);
    std::string sequence_file = std::string(argv[4]);

    cerr << "Reading the rindex file" << endl;

#if TIME
    auto time1 = chrono::high_resolution_clock::now();
#endif

    FastLocate r_index;
    {
        std::ifstream rin(r_index_file, std::ios::binary);
        if (!rin) {
            std::cerr << "Cannot open r-index: " << r_index_file << std::endl;
            return 1;
        }
        r_index.load_encoded(rin);
    }

#if TIME
    auto time2 = chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration1 = time2 - time1;
    std::cerr << "Loading r-index into memory took " << duration1.count() << " seconds" << std::endl;
#endif

    cerr << "Reading the tag array index" << endl;
    TagArray tag_array;
    std::ifstream in_ds(tag_array_index, std::ios::binary);
    tag_array.load_compressed_tags_compact(in_ds);

#if TIME
    auto time3 = chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration2 = time3 - time2;
    std::cerr << "Loading tag arrays took " << duration2.count() << " seconds" << std::endl;
#endif

    cerr << "Sample kmers created" << endl;

    std::vector<int> k_values = {10, 15, 20, 30, 50, 100, 200, 500, 1000, 2000};

    for (int k_val : k_values) {
        cerr << "k is: " << k_val << endl;
        std::vector<std::string> kmers;
        sampleKmersFromFile(sequence_file, k_val, 10000, kmers);
        cerr << "Finished creating the test kmers" << endl;

        int count = 0;

        // Shuffling the kmers
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(kmers.begin(), kmers.end(), g);

        int size_kmers = kmers.size();
        size_t total_number_of_tag_runs = 0;
        size_t total_query_intervals = 0;

#if TIME
        auto time10 = chrono::high_resolution_clock::now();
#endif

        double time_r_index = 0.0;
        double time_tag_query = 0.0;

        // For the first 10000 kmers
        for (int i = 0; i < 10000 && i < size_kmers; i++) {
            size_t tag_nums = 0;
            count++;

            string kmer = kmers[i];

            auto start_rindex = chrono::high_resolution_clock::now();
            gbwt::range_type range;
            if (r_index.is_encoded()) {
                range = r_index.count_encoded(kmer);
            } else {
                range = r_index.count(kmer);
            }
            auto end_rindex = chrono::high_resolution_clock::now();
            time_r_index += chrono::duration<double>(end_rindex - start_rindex).count();

            if (range.first > range.second) {
                continue; // Skip if no matches
            }

            auto start_tagquery = chrono::high_resolution_clock::now();
            tag_array.query_compressed_compact(range.first, range.second, tag_nums);
            auto end_tagquery = chrono::high_resolution_clock::now();
            time_tag_query += chrono::duration<double>(end_tagquery - start_tagquery).count();

            total_number_of_tag_runs += tag_nums;
            total_query_intervals += (range.second - range.first + 1);
        }

        cerr << "All queries done" << endl;

#if TIME
        auto time11 = chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration10 = time11 - time10;
        std::cerr << "Querying " << count << " kmers with size " << k_val << " took " << duration10.count() << " seconds" << std::endl;
        std::cerr << "r_index took " << time_r_index << " seconds in total" << std::endl;
        std::cerr << "tag_array took " << time_tag_query << " seconds in total" << std::endl;
        std::cerr << "Total number of tag runs: " << total_number_of_tag_runs << std::endl;
        std::cerr << "Total length of query intervals: " << total_query_intervals << std::endl;
        if (total_query_intervals > 0) {
            std::cerr << "Tag runs / query interval: " << (double) total_number_of_tag_runs / total_query_intervals << std::endl;
        }
#endif
    }

    return 0;
}

