//
// Created by Parsa Eskandar on 10/10/24.
//


#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <random>

#include "../include/pangenome_index/r-index.hpp"
#include "../deps/grlBWT/scripts/fm_index.h"
#include "../include/pangenome_index/algorithm.hpp"


using namespace panindexer;

namespace {

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

struct Rotation {
    std::string rotation;
    int seq_index;
    size_t start_pos;  // Starting position in original string (SA value)
};

// Function to create BWT from a given string and track the sequence index
std::pair <std::string, std::vector<panindexer::FastLocate::size_type>>
createBWTWithSequenceInfo(const std::string &input, const std::vector<int> &seq_indices) {
    int n = input.size();
    std::vector <Rotation> rotations;

    // Generate all rotations of the input string and associate each rotation with a sequence index
    for (int i = 0; i < n; ++i) {
        rotations.push_back({input.substr(i) + input.substr(0, i), seq_indices[i], static_cast<size_t>(i)});
    }

    std::sort(rotations.begin(), rotations.end(), [](const Rotation &a, const Rotation &b) {
        return a.rotation < b.rotation;
    });

    // Create BWT by taking the last column of sorted rotations
    std::string bwt;
    std::vector<panindexer::FastLocate::size_type> result_indices;

    for (const auto &rotation: rotations) {
        bwt += rotation.rotation.back();             // Collect BWT
        result_indices.push_back(rotation.seq_index); // Collect corresponding sequence index
    }


    return {bwt, result_indices};
}

// Create BWT and return (bwt, result_indices, sa_values) where sa_values[i] = text position for BWT position i
std::tuple<std::string, std::vector<panindexer::FastLocate::size_type>, std::vector<panindexer::FastLocate::size_type>>
createBWTWithSA(const std::string &input, const std::vector<int> &seq_indices) {
    int n = input.size();
    std::vector<Rotation> rotations;

    for (int i = 0; i < n; ++i) {
        rotations.push_back({input.substr(i) + input.substr(0, i), seq_indices[i], static_cast<size_t>(i)});
    }

    std::sort(rotations.begin(), rotations.end(), [](const Rotation &a, const Rotation &b) {
        return a.rotation < b.rotation;
    });

    std::string bwt;
    std::vector<panindexer::FastLocate::size_type> result_indices;
    std::vector<panindexer::FastLocate::size_type> sa_values;

    for (const auto &rotation : rotations) {
        bwt += rotation.rotation.back();
        result_indices.push_back(rotation.seq_index);
        sa_values.push_back(rotation.start_pos);
    }

    return {bwt, result_indices, sa_values};
}

// ======================= Additional encoded consistency tests =======================

TEST(RINDEX_Test, Count_Consistency_Encoded) {
    std::string text_file = "../test_data/big_test/merged_info";
    std::string rlbwt_file = "../test_data/big_test/merged_info.rl_bwt";

    FastLocate legacy(rlbwt_file);
    FastLocate encoded;
    {
        std::string enc_path = "./tmp_big_test.ri";
        std::ofstream out(enc_path, std::ios::binary | std::ios::trunc);
        ASSERT_TRUE(out.good());
        legacy.serialize_encoded(out);
        out.close();
        std::ifstream in(enc_path, std::ios::binary);
        ASSERT_TRUE(in.good());
        encoded.load_encoded(in);
    }

    std::ifstream in_txt(text_file);
    ASSERT_TRUE(in_txt.is_open()) << "Could not open input text file";
    std::string full_text, line;
    while (std::getline(in_txt, line)) { if (!line.empty()) full_text += line; }
    ASSERT_GE(full_text.size(), 100) << "Input text too short";

    std::mt19937 rng(12345);
    std::uniform_int_distribution<size_t> dist_pos(0, full_text.size() - 21);
    std::uniform_int_distribution<size_t> dist_len(5, 20);

    for (size_t i = 0; i < 100; ++i) {
        size_t pos = dist_pos(rng);
        size_t len = dist_len(rng);
        std::string pat = full_text.substr(pos, len);
        auto r1 = legacy.count(pat);
        auto r2 = encoded.count_encoded(pat);
        ASSERT_EQ(r1.first, r2.first) << "Mismatch in range start for pattern '" << pat << "'";
        ASSERT_EQ(r1.second, r2.second) << "Mismatch in range end for pattern '" << pat << "'";
    }
}

TEST(RINDEX_Test, BackwardExtend_Consistency_Encoded) {
    std::string rlbwt_file = "../test_data/big_test/merged_info.rl_bwt";
    FastLocate legacy(rlbwt_file);
    FastLocate encoded;
    {
        std::string enc_path = "./tmp_big_test.ri";
        std::ofstream out(enc_path, std::ios::binary | std::ios::trunc);
        ASSERT_TRUE(out.good());
        legacy.serialize_encoded(out);
        out.close();
        std::ifstream in(enc_path, std::ios::binary);
        ASSERT_TRUE(in.good());
        encoded.load_encoded(in);
    }

    std::ifstream in_txt("../test_data/big_test/merged_info");
    ASSERT_TRUE(in_txt.is_open()) << "Could not open input text file";
    std::string full_text, line;
    while (std::getline(in_txt, line)) { if (!line.empty()) full_text += line; }
    ASSERT_GE(full_text.size(), 50) << "Input text too short";

    std::mt19937 rng(6789);
    std::uniform_int_distribution<size_t> dist_pos(0, full_text.size() - 21);
    const size_t k = 20;

    for (size_t i = 0; i < 50; ++i) {
        size_t pos = dist_pos(rng);
        std::string kmer = full_text.substr(pos, k);

        panindexer::FastLocate::bi_interval b1 = {0, 0, legacy.bwt_size()};
        panindexer::FastLocate::bi_interval b2 = {0, 0, encoded.bwt_size()};
        for (int j = static_cast<int>(kmer.size()) - 1; j >= 0; --j) {
            b1 = legacy.backward_extend(b1, kmer[j]);
            b2 = encoded.backward_extend_encoded(b2, kmer[j]);
            ASSERT_EQ(b1.size, b2.size) << "Step mismatch at j=" << j << " for kmer '" << kmer << "'";
            if (b1.size == 0) break;
            ASSERT_EQ(b1.forward, b2.forward) << "Forward start mismatch at j=" << j;
        }
    }
}

}

// Function to read strings from a file, concatenate them using $_i, and create BWT with sequence index tracking
std::vector<panindexer::FastLocate::size_type> readFileAndCreateBWTWithIndices(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the file!" << std::endl;
        return {"", {}};
    }

    std::string line;
    std::string concatenatedString;
    std::vector<int> seq_indices;
    int sequence_number = 0;


    char added_char = '$';
    while (std::getline(file, line)) {


        // Mark the characters of the current string as belonging to this sequence
        for (size_t i = 0; i < line.size(); ++i) {
            seq_indices.push_back(sequence_number);
        }

        // Add $ between strings

        seq_indices.push_back(sequence_number);
        concatenatedString += line;
        concatenatedString += added_char;


        added_char += 1;
        ++sequence_number;
    }

    file.close();

    return createBWTWithSequenceInfo(concatenatedString, seq_indices).second;
}



TEST(RINDEX_Test, Locate_small_test) {
    std::string filename = "../test_data/small_test_nl.txt";
    auto sequence_indices = readFileAndCreateBWTWithIndices(filename);


    std::string rlbwt_file = "../test_data/small_test_nl.rl_bwt";
    FastLocate r_index(rlbwt_file);
    std::cerr << "r_index size: " << r_index.size() << std::endl;

    // print header max length
    std::cerr << "header max length: " << r_index.header.max_length << std::endl;
    // print the samples 
    std::cerr << "samples: " << std::endl;
    for (size_t i = 0; i < r_index.samples.size(); i++) {
        std::cerr << r_index.samples[i] << " ";
    }
    std::cerr << std::endl;

    auto x = r_index.decompressDA();
    std::cerr << "x: " << std::endl;
    for (size_t i = 0; i < x.size(); i++) {
        std::cerr << x[i] << " ";
    }
    std::cerr << std::endl;
    std::cerr << "sequence_indices: " << std::endl;
    for (size_t i = 0; i < sequence_indices.size(); i++) {
        std::cerr << sequence_indices[i] << " ";
    }
    std::cerr << std::endl;
    ASSERT_EQ(x, sequence_indices) << "Invalid Locate results from the r-index";

}


TEST(RINDEX_Test, Locate_small_test_encoded) {
    std::string filename = "../test_data/small_test_nl.txt";
    auto sequence_indices = readFileAndCreateBWTWithIndices(filename);

    std::string rlbwt_file = "../test_data/small_test_nl.rl_bwt";
    // Build legacy, serialize encoded, load encoded
    FastLocate built(rlbwt_file);
    std::string enc_path = "./tmp_small_test.ri";
    {
        std::ofstream out(enc_path, std::ios::binary | std::ios::trunc);
        ASSERT_TRUE(out.good());
        built.serialize_encoded(out);
    }
    FastLocate r_index;
    {
        std::ifstream in(enc_path, std::ios::binary);
        ASSERT_TRUE(in.good());
        r_index.load_encoded(in);
    }

    auto x = r_index.decompressDA_encoded();
    ASSERT_EQ(x, sequence_indices) << "Invalid Locate results from the encoded r-index";
}

TEST(RINDEX_Test, Locate_medium_test) {
    std::string filename = "../test_data/med_test.txt";
    auto sequence_indices = readFileAndCreateBWTWithIndices(filename);

    std::string rlbwt_file = "../test_data/med_test.rl_bwt";
    FastLocate r_index(rlbwt_file);

    auto x = r_index.decompressDA();
    ASSERT_EQ(x, sequence_indices) << "Invalid Locate results from the r-index";

}

TEST(RINDEX_Test, Locate_medium_test_encoded) {
    std::string filename = "../test_data/med_test.txt";
    auto sequence_indices = readFileAndCreateBWTWithIndices(filename);

    std::string rlbwt_file = "../test_data/med_test.rl_bwt";
    FastLocate built(rlbwt_file);
    std::string enc_path = "./tmp_med_test.ri";
    {
        std::ofstream out(enc_path, std::ios::binary | std::ios::trunc);
        ASSERT_TRUE(out.good());
        built.serialize_encoded(out);
    }
    FastLocate r_index;
    {
        std::ifstream in(enc_path, std::ios::binary);
        ASSERT_TRUE(in.good());
        r_index.load_encoded(in);
    }

    auto x = r_index.decompressDA_encoded();
    ASSERT_EQ(x, sequence_indices) << "Invalid Locate results from the encoded r-index";
}

TEST(RINDEX_Test, Locate_big_test) {
    std::string filename = "../test_data/x.newline_separated";
    auto sequence_indices = readFileAndCreateBWTWithIndices(filename);

    std::string rlbwt_file = "../test_data/x.rl_bwt";
    FastLocate r_index(rlbwt_file);

    auto x = r_index.decompressDA();
    ASSERT_EQ(x, sequence_indices) << "Invalid Locate results from the r-index";

}

TEST(RINDEX_Test, Locate_big_test_encoded) {
    std::string filename = "../test_data/x.newline_separated";
    auto sequence_indices = readFileAndCreateBWTWithIndices(filename);

    std::string rlbwt_file = "../test_data/x.rl_bwt";
    FastLocate built(rlbwt_file);
    std::string enc_path = "./tmp_big_test.ri";
    {
        std::ofstream out(enc_path, std::ios::binary | std::ios::trunc);
        ASSERT_TRUE(out.good());
        built.serialize_encoded(out);
    }
    FastLocate r_index;
    {
        std::ifstream in(enc_path, std::ios::binary);
        ASSERT_TRUE(in.good());
        r_index.load_encoded(in);
    }

    auto x = r_index.decompressDA_encoded();
    ASSERT_EQ(x, sequence_indices) << "Invalid Locate results from the encoded r-index";
}

TEST(RINDEX_Test, Locate_N_test) {
    std::string filename = "../test_data/N_test.txt";
    auto sequence_indices = readFileAndCreateBWTWithIndices(filename);

    std::string rlbwt_file = "../test_data/N_test.rl_bwt";
    FastLocate r_index(rlbwt_file);

    auto x = r_index.decompressDA();
    ASSERT_EQ(x, sequence_indices) << "Invalid Locate results from the r-index";

}

TEST(RINDEX_Test, Locate_N_test_encoded) {
    std::string filename = "../test_data/N_test.txt";
    auto sequence_indices = readFileAndCreateBWTWithIndices(filename);

    std::string rlbwt_file = "../test_data/N_test.rl_bwt";
    FastLocate built(rlbwt_file);
    std::string enc_path = "./tmp_N_test.ri";
    {
        std::ofstream out(enc_path, std::ios::binary | std::ios::trunc);
        ASSERT_TRUE(out.good());
        built.serialize_encoded(out);
    }
    FastLocate r_index;
    {
        std::ifstream in(enc_path, std::ios::binary);
        ASSERT_TRUE(in.good());
        r_index.load_encoded(in);
    }

    auto x = r_index.decompressDA_encoded();
    ASSERT_EQ(x, sequence_indices) << "Invalid Locate results from the encoded r-index";
}

// TEST(FMINDEX_Test, small_test){
//     std::string filename = "../big_test/merged_info.rl_bwt";
//     fm_index index(filename);
//     FastLocate r_index(filename);

//     for (int j = 0; j < 3; j++){
//         auto start_fm = j;
//         auto start_r = j;
//         auto a = index.lf(start_fm);
//         auto b = r_index.psi(start_r);

//         start_fm = a.second;
//         auto fm_char = a.first;
//         start_r = b.second;
//         auto r_char = b.first;

//         int i = 0;
//         while (r_char != NENDMARKER && fm_char != NENDMARKER){
// //                std::cerr << i << " " << start_fm << " " << start_r << std::endl;
//             ASSERT_EQ(fm_char, r_char) << "Invalid LF results from the FM-index and the r-index";
//             ASSERT_EQ(start_fm, start_r) << "Invalid LF results from the FM-index and the r-index";
//             a = index.lf(start_fm);
//             b = r_index.psi(start_r);
//             start_fm = a.second;
//             fm_char = a.first;
//             start_r = b.second;
//             r_char = b.first;
//             i++;
//         }



//     }



// }

// ======================= Encoded vs legacy parity (hard tests) =======================
// Load legacy from .rl_bwt, serialize to encoded, load encoded; then compare core operations.

static void load_legacy_and_encoded(const std::string& rlbwt_file, const std::string& temp_ri,
                                     FastLocate& legacy, FastLocate& encoded) {
    legacy = FastLocate(rlbwt_file);
    {
        std::ofstream out(temp_ri, std::ios::binary | std::ios::trunc);
        ASSERT_TRUE(out.good());
        legacy.serialize_encoded(out);
    }
    {
        std::ifstream in(temp_ri, std::ios::binary);
        ASSERT_TRUE(in.good());
        encoded.load_encoded(in);
    }
    ASSERT_FALSE(legacy.is_encoded());
    ASSERT_TRUE(encoded.is_encoded());
    ASSERT_EQ(legacy.bwt_size(), encoded.bwt_size());
    ASSERT_EQ(legacy.tot_runs(), encoded.tot_runs());
}

TEST(RINDEX_Test, EncodedVsLegacy_LF_Range_Small) {
    const char* rlbwt = "../test_data/small_test_nl.rl_bwt";
    FastLocate legacy, encoded;
    load_legacy_and_encoded(rlbwt, "./tmp_parity_small.ri", legacy, encoded);
    const size_t n = legacy.bwt_size();
    std::mt19937 rng(1001);
    const size_t num_trials = 200;
    for (size_t t = 0; t < num_trials; ++t) {
        size_t first = n ? (rng() % n) : 0;
        size_t second = n ? (rng() % n) : 0;
        if (first > second) std::swap(first, second);
        for (unsigned char sym : {'A', 'C', 'G', 'T'}) {
            gbwt::range_type leg = legacy.LF({first, second}, sym);
            gbwt::range_type enc = encoded.LF_encoded({first, second}, sym);
            ASSERT_EQ(leg.first, enc.first) << "LF range first mismatch pos=[" << first << "," << second << "] sym=" << (char)sym;
            ASSERT_EQ(leg.second, enc.second) << "LF range second mismatch pos=[" << first << "," << second << "] sym=" << (char)sym;
        }
    }
}

TEST(RINDEX_Test, EncodedVsLegacy_LF_Range_Medium) {
    const char* rlbwt = "../test_data/med_test.rl_bwt";
    FastLocate legacy, encoded;
    load_legacy_and_encoded(rlbwt, "./tmp_parity_med.ri", legacy, encoded);
    const size_t n = legacy.bwt_size();
    std::mt19937 rng(1002);
    for (size_t t = 0; t < 300; ++t) {
        size_t first = n ? (rng() % n) : 0;
        size_t second = n ? (rng() % n) : 0;
        if (first > second) std::swap(first, second);
        for (unsigned char sym : {'A', 'C', 'G', 'T'}) {
            gbwt::range_type leg = legacy.LF({first, second}, sym);
            gbwt::range_type enc = encoded.LF_encoded({first, second}, sym);
            ASSERT_EQ(leg.first, enc.first) << "LF range first mismatch sym=" << (char)sym;
            ASSERT_EQ(leg.second, enc.second) << "LF range second mismatch sym=" << (char)sym;
        }
    }
}

TEST(RINDEX_Test, EncodedVsLegacy_LF_Range_NTest) {
    const char* rlbwt = "../test_data/N_test.rl_bwt";
    FastLocate legacy, encoded;
    load_legacy_and_encoded(rlbwt, "./tmp_parity_N.ri", legacy, encoded);
    const size_t n = legacy.bwt_size();
    std::mt19937 rng(1003);
    for (size_t t = 0; t < 200; ++t) {
        size_t first = n ? (rng() % n) : 0;
        size_t second = n ? (rng() % n) : 0;
        if (first > second) std::swap(first, second);
        for (unsigned char sym : {'A', 'C', 'G', 'T', 'N'}) {
            gbwt::range_type leg = legacy.LF({first, second}, sym);
            gbwt::range_type enc = encoded.LF_encoded({first, second}, sym);
            ASSERT_EQ(leg.first, enc.first) << "LF range first mismatch sym=" << (char)sym;
            ASSERT_EQ(leg.second, enc.second) << "LF range second mismatch sym=" << (char)sym;
        }
    }
}

TEST(RINDEX_Test, EncodedVsLegacy_RankAt_Small) {
    const char* rlbwt = "../test_data/small_test_nl.rl_bwt";
    FastLocate legacy, encoded;
    load_legacy_and_encoded(rlbwt, "./tmp_parity_small.ri", legacy, encoded);
    const size_t n = legacy.bwt_size();
    std::mt19937 rng(2001);
    for (size_t t = 0; t < 300; ++t) {
        size_t pos = n ? (rng() % (n + 1)) : 0;  // 0..n inclusive (rank at n = total count)
        for (unsigned char sym : {'A', 'C', 'G', 'T'}) {
            size_t r_leg = legacy.rankAt(pos, sym);
            size_t r_enc = encoded.rankAt(pos, sym);
            ASSERT_EQ(r_leg, r_enc) << "rankAt mismatch pos=" << pos << " sym=" << (char)sym;
        }
    }
}

TEST(RINDEX_Test, EncodedVsLegacy_RankAt_Medium) {
    const char* rlbwt = "../test_data/med_test.rl_bwt";
    FastLocate legacy, encoded;
    load_legacy_and_encoded(rlbwt, "./tmp_parity_med.ri", legacy, encoded);
    const size_t n = legacy.bwt_size();
    std::mt19937 rng(2002);
    for (size_t t = 0; t < 400; ++t) {
        size_t pos = n ? (rng() % (n + 1)) : 0;
        for (unsigned char sym : {'A', 'C', 'G', 'T'}) {
            ASSERT_EQ(legacy.rankAt(pos, sym), encoded.rankAt(pos, sym))
                << "rankAt pos=" << pos << " sym=" << (char)sym;
        }
    }
}

TEST(RINDEX_Test, EncodedVsLegacy_RunIdAndGetSample_Small) {
    const char* rlbwt = "../test_data/small_test_nl.rl_bwt";
    FastLocate legacy, encoded;
    load_legacy_and_encoded(rlbwt, "./tmp_parity_small.ri", legacy, encoded);
    const size_t n = legacy.bwt_size();
    const size_t num_runs = legacy.tot_runs();
    for (size_t run_id = 0; run_id < num_runs; ++run_id) {
        ASSERT_EQ(legacy.getSample(run_id), encoded.getSample(run_id))
            << "getSample mismatch run_id=" << run_id;
    }
    std::mt19937 rng(3001);
    for (size_t t = 0; t < 200; ++t) {
        size_t pos = n ? (rng() % n) : 0;
        size_t run_leg, run_enc, off_leg, off_enc;
        legacy.run_id_and_offset_at(pos, run_leg, off_leg);
        encoded.run_id_and_offset_at(pos, run_enc, off_enc);
        ASSERT_EQ(run_leg, run_enc) << "run_id_and_offset_at run_id mismatch pos=" << pos;
        ASSERT_EQ(off_leg, off_enc) << "run_id_and_offset_at offset mismatch pos=" << pos;
        ASSERT_EQ(legacy.getSample(run_leg), encoded.getSample(run_enc)) << "getSample at pos=" << pos;
    }
}

TEST(RINDEX_Test, EncodedVsLegacy_RunIdAndGetSample_Medium) {
    const char* rlbwt = "../test_data/med_test.rl_bwt";
    FastLocate legacy, encoded;
    load_legacy_and_encoded(rlbwt, "./tmp_parity_med.ri", legacy, encoded);
    const size_t n = legacy.bwt_size();
    const size_t num_runs = legacy.tot_runs();
    for (size_t run_id = 0; run_id < num_runs; ++run_id) {
        ASSERT_EQ(legacy.getSample(run_id), encoded.getSample(run_id)) << "run_id=" << run_id;
    }
    std::mt19937 rng(3002);
    for (size_t t = 0; t < 400; ++t) {
        size_t pos = n ? (rng() % n) : 0;
        size_t run_leg, run_enc, off_leg, off_enc;
        legacy.run_id_and_offset_at(pos, run_leg, off_leg);
        encoded.run_id_and_offset_at(pos, run_enc, off_enc);
        ASSERT_EQ(run_leg, run_enc);
        ASSERT_EQ(off_leg, off_enc);
    }
}

TEST(RINDEX_Test, EncodedVsLegacy_SingleLF_Small) {
    const char* rlbwt = "../test_data/small_test_nl.rl_bwt";
    FastLocate legacy, encoded;
    load_legacy_and_encoded(rlbwt, "./tmp_parity_small.ri", legacy, encoded);
    const size_t n = legacy.bwt_size();
    std::mt19937 rng(4001);
    for (size_t t = 0; t < 200; ++t) {
        size_t pos = n ? (rng() % n) : 0;
        size_t next_leg = legacy.LF(pos);
        size_t next_enc = encoded.LF(pos);
        ASSERT_EQ(next_leg, next_enc) << "LF(idx) mismatch at pos=" << pos;
    }
}

TEST(RINDEX_Test, EncodedVsLegacy_SingleLF_Medium) {
    const char* rlbwt = "../test_data/med_test.rl_bwt";
    FastLocate legacy, encoded;
    load_legacy_and_encoded(rlbwt, "./tmp_parity_med.ri", legacy, encoded);
    const size_t n = legacy.bwt_size();
    std::mt19937 rng(4002);
    for (size_t t = 0; t < 500; ++t) {
        size_t pos = n ? (rng() % n) : 0;
        ASSERT_EQ(legacy.LF(pos), encoded.LF(pos)) << "pos=" << pos;
    }
}

TEST(RINDEX_Test, EncodedVsLegacy_DecompressDA_Small) {
    const char* rlbwt = "../test_data/small_test_nl.rl_bwt";
    FastLocate legacy, encoded;
    load_legacy_and_encoded(rlbwt, "./tmp_parity_small.ri", legacy, encoded);
    auto da_leg = legacy.decompressDA();
    auto da_enc = encoded.decompressDA();
    ASSERT_EQ(da_leg.size(), da_enc.size());
    for (size_t i = 0; i < da_leg.size(); ++i) {
        ASSERT_EQ(da_leg[i], da_enc[i]) << "decompressDA mismatch at i=" << i;
    }
}

TEST(RINDEX_Test, EncodedVsLegacy_DecompressDA_Medium) {
    const char* rlbwt = "../test_data/med_test.rl_bwt";
    FastLocate legacy, encoded;
    load_legacy_and_encoded(rlbwt, "./tmp_parity_med.ri", legacy, encoded);
    auto da_leg = legacy.decompressDA();
    auto da_enc = encoded.decompressDA();
    ASSERT_EQ(da_leg.size(), da_enc.size());
    for (size_t i = 0; i < da_leg.size(); ++i) {
        ASSERT_EQ(da_leg[i], da_enc[i]) << "i=" << i;
    }
}

TEST(RINDEX_Test, EncodedVsLegacy_DecompressDA_NTest) {
    const char* rlbwt = "../test_data/N_test.rl_bwt";
    FastLocate legacy, encoded;
    load_legacy_and_encoded(rlbwt, "./tmp_parity_N.ri", legacy, encoded);
    auto da_leg = legacy.decompressDA();
    auto da_enc = encoded.decompressDA();
    ASSERT_EQ(da_leg.size(), da_enc.size());
    for (size_t i = 0; i < da_leg.size(); ++i) {
        ASSERT_EQ(da_leg[i], da_enc[i]) << "i=" << i;
    }
}

TEST(RINDEX_Test, EncodedVsLegacy_BackwardExtend_StepByStep_Small) {
    const char* rlbwt = "../test_data/small_test_nl.rl_bwt";
    FastLocate legacy, encoded;
    load_legacy_and_encoded(rlbwt, "./tmp_parity_small.ri", legacy, encoded);
    std::string pattern = "ACGTACGT";
    panindexer::FastLocate::bi_interval b_leg = {0, 0, legacy.bwt_size()};
    panindexer::FastLocate::bi_interval b_enc = {0, 0, encoded.bwt_size()};
    for (char c : pattern) {
        b_leg = legacy.backward_extend(b_leg, static_cast<size_t>(c));
        b_enc = encoded.backward_extend_encoded(b_enc, static_cast<size_t>(c));
        ASSERT_EQ(b_leg.forward, b_enc.forward) << "backward_extend forward after '" << c << "'";
        ASSERT_EQ(b_leg.reverse, b_enc.reverse) << "backward_extend reverse after '" << c << "'";
        ASSERT_EQ(b_leg.size, b_enc.size) << "backward_extend size after '" << c << "'";
    }
}

TEST(RINDEX_Test, EncodedVsLegacy_BackwardExtend_StepByStep_Medium) {
    const char* rlbwt = "../test_data/med_test.rl_bwt";
    FastLocate legacy, encoded;
    load_legacy_and_encoded(rlbwt, "./tmp_parity_med.ri", legacy, encoded);
    std::mt19937 rng(5001);
    const char alphabet[] = "ACGT";
    for (size_t trial = 0; trial < 100; ++trial) {
        std::string pattern;
        for (int i = 0; i < 15; ++i) pattern += alphabet[rng() % 4];
        panindexer::FastLocate::bi_interval b_leg = {0, 0, legacy.bwt_size()};
        panindexer::FastLocate::bi_interval b_enc = {0, 0, encoded.bwt_size()};
        for (char c : pattern) {
            b_leg = legacy.backward_extend(b_leg, static_cast<size_t>(c));
            b_enc = encoded.backward_extend_encoded(b_enc, static_cast<size_t>(c));
            if (b_leg.size == 0) break;
            ASSERT_EQ(b_leg.forward, b_enc.forward);
            ASSERT_EQ(b_leg.reverse, b_enc.reverse);
            ASSERT_EQ(b_leg.size, b_enc.size);
        }
    }
}

TEST(RINDEX_Test, EncodedVsLegacy_Count_BigTest) {
    const char* rlbwt = "../test_data/big_test/merged_info.rl_bwt";
    FastLocate legacy, encoded;
    load_legacy_and_encoded(rlbwt, "./tmp_parity_big.ri", legacy, encoded);
    std::ifstream in("../test_data/big_test/merged_info");
    ASSERT_TRUE(in.is_open());
    std::string full_text, line;
    while (std::getline(in, line)) { if (!line.empty()) full_text += line; }
    in.close();
    ASSERT_GE(full_text.size(), 50u);
    std::mt19937 rng(6001);
    for (size_t t = 0; t < 100; ++t) {
        size_t pos = rng() % (full_text.size() - 20);
        std::string pat = full_text.substr(pos, 10 + (rng() % 10));
        auto r_leg = legacy.count(pat);
        auto r_enc = encoded.count_encoded(pat);
        ASSERT_EQ(r_leg.first, r_enc.first) << "count start pattern '" << pat << "'";
        ASSERT_EQ(r_leg.second, r_enc.second) << "count end pattern '" << pat << "'";
    }
}

TEST(FMDINDEX_Test, BackwardExtensionMatchesLF) {
    std::string rlbwt_file = "../test_data/big_test/merged_info.rl_bwt";
    FastLocate r_index(rlbwt_file);
//    r_index.initialize_complement_table();

    std::string kmer = "ATCAAAGAAAAAAGCCCAACATATCCATTACCATTACTAGTTACACATAGCATCAGGAACCAGAGAGTTGGA";
    std::string revcomp = kmer;
//    std::reverse(revcomp.begin(), revcomp.end());
//    for (char& c : revcomp) {
//        c = r_index.complement(c);
//    }

//    std::cerr << "Testing kmer: " << kmer << " and its reverse complement: " << revcomp << std::endl;

    // Initial full-range interval
    panindexer::FastLocate::bi_interval bint = {0, 0, r_index.bwt_size()};

    for (int i = kmer.size() - 1; i >= 0; --i) {
        char fwd_char = kmer[i];
        char rev_char = revcomp[i];

        // Compute manually via LF
        range_type fwd_expected = r_index.LF({bint.forward, bint.forward + bint.size - 1}, fwd_char);
        range_type rev_expected = r_index.LF({bint.reverse, bint.reverse + bint.size}, rev_char);

        size_t expected_size = 0;
        if (fwd_expected.first <= fwd_expected.second) {
            expected_size = fwd_expected.second - fwd_expected.first + 1;
        }

        panindexer::FastLocate::bi_interval extended = r_index.backward_extend(bint, fwd_char);


//        std::cerr << "char fwd: " << fwd_char << " char rev: " << rev_char << ", Interval size: " << extended.size << " Interval reverse " <<  extended.reverse << " Interval forward " << extended.forward << "\n";
        // Debug
//        std::cerr << "Char: " << fwd_char
//        << " | Expected size: " << expected_size
//        << " | Actual size: " << extended.size << std::endl;

        if (expected_size > 0){
            ASSERT_EQ(extended.size, expected_size) << "Mismatch in size for char " << fwd_char;
            ASSERT_EQ(extended.forward, fwd_expected.first);
        }

//        ASSERT_EQ(extended.reverse, rev_expected.first);

        // Move to next iteration
        bint = extended;
    }
}

TEST(FMDINDEX_Test, BackwardExtensionMatchesLF_Encoded) {
    std::string rlbwt_file = "../test_data/big_test/merged_info.rl_bwt";
    // Build and load encoded
    FastLocate built(rlbwt_file);
    std::string enc_path = "./tmp_bwd_ext.ri";
    {
        std::ofstream out(enc_path, std::ios::binary | std::ios::trunc);
        ASSERT_TRUE(out.good());
        built.serialize_encoded(out);
    }
    FastLocate r_index;
    {
        std::ifstream in(enc_path, std::ios::binary);
        ASSERT_TRUE(in.good());
        r_index.load_encoded(in);
    }

    std::string kmer = "ATCAAAGAAAAAAGCCCAACATATCCATTACCATTACTAGTTACACATAGCATCAGGAACCAGAGAGTTGGA";
    panindexer::FastLocate::bi_interval bint = {0, 0, r_index.bwt_size()};

    for (int i = kmer.size() - 1; i >= 0; --i) {
        char fwd_char = kmer[i];
        range_type fwd_expected = r_index.LF_encoded({bint.forward, bint.forward + bint.size - 1}, fwd_char);
        size_t expected_size = 0;
        if (fwd_expected.first <= fwd_expected.second) {
            expected_size = fwd_expected.second - fwd_expected.first + 1;
        }
        panindexer::FastLocate::bi_interval extended = r_index.backward_extend_encoded(bint, fwd_char);
        if (expected_size > 0){
            ASSERT_EQ(extended.size, expected_size) << "Mismatch in size for char " << fwd_char;
            ASSERT_EQ(extended.forward, fwd_expected.first);
        }
        bint = extended;
    }
}


TEST(FMDINDEX_Test, CompareSampledKmersWithReverseComplementsBIGTEST) {
    std::string rlbwt_file = "../test_data/big_test/merged_info.rl_bwt";
    std::string text_file = "../test_data/big_test/merged_info";

    FastLocate r_index(rlbwt_file);
    r_index.initialize_complement_table();

    std::ifstream in(text_file);
    ASSERT_TRUE(in.is_open()) << "Could not open input text file";

    std::string full_text;
    std::string line;
    while (std::getline(in, line)) {
    if (!line.empty()) full_text += line;
    }

    ASSERT_GE(full_text.size(), 20) << "Input text is too short to sample from";

    const size_t k = 12;
    const size_t num_samples = 100;
    std::mt19937 rng(42); // fixed seed for reproducibility
    std::uniform_int_distribution<size_t> dist(0, full_text.size() - k);

    for (size_t i = 0; i < num_samples; ++i) {
    size_t pos = dist(rng);
    std::string kmer = full_text.substr(pos, k);

    // Skip invalid kmers with non-ACGTN symbols
    if (kmer.find_first_not_of("ACGTN") != std::string::npos) {
    i--;
    continue;
    }

    std::string revcomp = kmer;
    std::reverse(revcomp.begin(), revcomp.end());
    for (char& c : revcomp) {
    c = r_index.complement(c);
    }

    // Backward extend original
    panindexer::FastLocate::bi_interval int_kmer = {0, 0, r_index.bwt_size()};
    for (int j = kmer.size() - 1; j >= 0; --j) {
    int_kmer = r_index.backward_extend(int_kmer, kmer[j]);
    }

    // Backward extend reverse complement
    panindexer::FastLocate::bi_interval int_rc = {0, 0, r_index.bwt_size()};
    for (int j = revcomp.size() - 1; j >= 0; --j) {
    int_rc = r_index.backward_extend(int_rc, revcomp[j]);
    }

//    std::cerr << "k-mer       : " << kmer << ", Interval size: " << int_kmer.size << " Interval reverse " <<  int_kmer.reverse << " Interval forward " << int_kmer.forward << "\n";
//    std::cerr << "RevComp     : " << revcomp << ", Interval size: " << int_rc.size << " Interval reverse " <<  int_rc.reverse << " Interval forward " << int_rc.forward << "\n";

    ASSERT_EQ(int_kmer.forward, int_rc.reverse) << "FMD symmetry violated between '" << kmer << "' and '" << revcomp << "'";
    ASSERT_EQ(int_kmer.reverse, int_rc.forward) << "FMD symmetry violated between '" << kmer << "' and '" << revcomp << "'";
    ASSERT_EQ(int_kmer.size, int_rc.size) << "FMD symmetry violated between '" << kmer << "' and '" << revcomp << "'";

    }
}

TEST(FMDINDEX_Test, CompareSampledKmersWithReverseComplementsBIGTEST_Encoded) {
    std::string rlbwt_file = "../test_data/big_test/merged_info.rl_bwt";
    FastLocate built(rlbwt_file);
    std::string enc_path = "./tmp_big_rc.ri";
    {
        std::ofstream out(enc_path, std::ios::binary | std::ios::trunc);
        ASSERT_TRUE(out.good());
        built.serialize_encoded(out);
    }
    FastLocate r_index;
    {
        std::ifstream in(enc_path, std::ios::binary);
        ASSERT_TRUE(in.good());
        r_index.load_encoded(in);
    }

    r_index.initialize_complement_table();

    std::ifstream in_txt("../test_data/big_test/merged_info");
    ASSERT_TRUE(in_txt.is_open()) << "Could not open input text file";

    std::string full_text, line;
    while (std::getline(in_txt, line)) { if (!line.empty()) full_text += line; }

    ASSERT_GE(full_text.size(), 20) << "Input text is too short to sample from";

    const size_t k = 12;
    const size_t num_samples = 100;
    std::mt19937 rng(42);
    std::uniform_int_distribution<size_t> dist(0, full_text.size() - k);

    for (size_t i = 0; i < num_samples; ++i) {
        size_t pos = dist(rng);
        std::string kmer = full_text.substr(pos, k);
        if (kmer.find_first_not_of("ACGTN") != std::string::npos) { i--; continue; }
        std::string revcomp = kmer;
        std::reverse(revcomp.begin(), revcomp.end());
        for (char& c : revcomp) { c = r_index.complement(c); }

        panindexer::FastLocate::bi_interval int_kmer = {0, 0, r_index.bwt_size()};
        for (int j = kmer.size() - 1; j >= 0; --j) { int_kmer = r_index.backward_extend_encoded(int_kmer, kmer[j]); }

        panindexer::FastLocate::bi_interval int_rc = {0, 0, r_index.bwt_size()};
        for (int j = revcomp.size() - 1; j >= 0; --j) { int_rc = r_index.backward_extend_encoded(int_rc, revcomp[j]); }

        ASSERT_EQ(int_kmer.forward, int_rc.reverse) << "FMD symmetry violated between '" << kmer << "' and '" << revcomp << "'";
        ASSERT_EQ(int_kmer.reverse, int_rc.forward) << "FMD symmetry violated between '" << kmer << "' and '" << revcomp << "'";
        ASSERT_EQ(int_kmer.size, int_rc.size) << "FMD symmetry violated between '" << kmer << "' and '" << revcomp << "'";
    }
}

TEST(FMDINDEX_Test, CompareSampledKmersNoN_Bidirectional) {
    std::string rlbwt_file = "../test_data/bidirectional_test/contigs_xy.rl_bwt";
    std::string text_file = "../test_data/bidirectional_test/contigs_xy";

    FastLocate r_index(rlbwt_file);
    r_index.initialize_complement_table();

    std::ifstream in(text_file);
    ASSERT_TRUE(in.is_open()) << "Could not open input text file";

    std::string full_text;
    std::string line;
    while (std::getline(in, line)) { if (!line.empty()) full_text += line; }

    ASSERT_GE(full_text.size(), 20) << "Input text is too short to sample from";

    const size_t k = 12;
    const size_t num_samples = 100;
    std::mt19937 rng(7);
    std::uniform_int_distribution<size_t> dist(0, full_text.size() - k);

    for (size_t i = 0; i < num_samples; ++i) {
        size_t pos = dist(rng);
        std::string kmer = full_text.substr(pos, k);
        // Ensure no 'N' in sampled kmer
        if (kmer.find('N') != std::string::npos) { i--; continue; }

        std::string revcomp = kmer;
        std::reverse(revcomp.begin(), revcomp.end());
        for (char& c : revcomp) { c = r_index.complement(c); }

        panindexer::FastLocate::bi_interval int_kmer = {0, 0, r_index.bwt_size()};
        for (int j = kmer.size() - 1; j >= 0; --j) { int_kmer = r_index.backward_extend(int_kmer, kmer[j]); }

        panindexer::FastLocate::bi_interval int_rc = {0, 0, r_index.bwt_size()};
        for (int j = revcomp.size() - 1; j >= 0; --j) { int_rc = r_index.backward_extend(int_rc, revcomp[j]); }

        ASSERT_EQ(int_kmer.forward, int_rc.reverse) << "FMD symmetry violated between '" << kmer << "' and '" << revcomp << "'";
        ASSERT_EQ(int_kmer.reverse, int_rc.forward) << "FMD symmetry violated between '" << kmer << "' and '" << revcomp << "'";
        ASSERT_EQ(int_kmer.size, int_rc.size) << "FMD symmetry violated between '" << kmer << "' and '" << revcomp << "'";
    }
}

TEST(RINDEX_Test, Samples_Offset_From_End) {
    std::string filename = "../test_data/small_test_nl.txt";
    std::ifstream file(filename);
    ASSERT_TRUE(file.is_open()) << "Could not open " << filename;

    std::string line;
    std::string concatenatedString;
    std::vector<int> seq_indices;
    int sequence_number = 0;
    char added_char = '$';

    while (std::getline(file, line)) {
        for (size_t i = 0; i < line.size(); ++i) {
            seq_indices.push_back(sequence_number);
        }
        seq_indices.push_back(sequence_number);
        concatenatedString += line;
        concatenatedString += added_char;
        added_char += 1;
        ++sequence_number;
    }
    file.close();

    std::cout << "[Samples_Offset_From_End] Input file: " << filename << std::endl;
    std::cout << "[Samples_Offset_From_End] Concatenated string length: " << concatenatedString.size()
              << " (\"" << concatenatedString << "\")" << std::endl;
    std::cout << "[Samples_Offset_From_End] Number of sequences: " << sequence_number << std::endl;

    ASSERT_GE(concatenatedString.size(), 1u) << "Empty or invalid test file";

    auto [bwt, result_indices, expected_sa] = createBWTWithSA(concatenatedString, seq_indices);
    ASSERT_EQ(expected_sa.size(), concatenatedString.size()) << "SA size mismatch";

    int n_seq = sequence_number;
    std::vector<size_t> seq_start(n_seq), seq_length(n_seq);
    for (int s = 0; s < n_seq; ++s) {
        seq_start[s] = 0;
        for (int t = 0; t < s; ++t) {
            size_t count = 0;
            for (int idx : seq_indices) { if (idx == t) ++count; }
            seq_start[s] += count;
        }
        seq_length[s] = 0;
        for (int idx : seq_indices) { if (idx == s) ++seq_length[s]; }
    }

    std::cout << "[Samples_Offset_From_End] Sequence boundaries:" << std::endl;
    for (int s = 0; s < n_seq; ++s) {
        std::cout << "  seq_id=" << s << " start=" << seq_start[s] << " length=" << seq_length[s] << std::endl;
    }

    std::string rlbwt_file = "../test_data/small_test_nl.rl_bwt";
    std::cout << "[Samples_Offset_From_End] Loading r-index: " << rlbwt_file << std::endl;
    FastLocate r_index(rlbwt_file);

    std::vector<panindexer::FastLocate::size_type> r_index_sa = r_index.decompressSA();
    ASSERT_EQ(r_index_sa.size(), expected_sa.size()) << "r-index SA size mismatch with expected";

    std::cout << "[Samples_Offset_From_End] BWT size: " << r_index_sa.size() << std::endl;
    std::cout << "[Samples_Offset_From_End] Comparing stored samples vs expected (offset-from-start vs offset-from-end):" << std::endl;

    // Detect whether samples use offset-from-START or offset-from-END
    bool matches_offset_from_end = true;
    bool matches_offset_from_start = true;
    for (size_t i = 0; i < r_index_sa.size(); ++i) {
        size_t packed = r_index_sa[i];
        size_t seq_id = r_index.seqId(packed);
        size_t seq_offset = r_index.seqOffset(packed);

        size_t text_pos = expected_sa[i];
        size_t expected_seq_id = static_cast<size_t>(seq_indices[text_pos]);
        size_t offset_in_seq = text_pos - seq_start[expected_seq_id];
        size_t seq_len = seq_length[expected_seq_id];
        size_t expected_offset_from_end = (seq_len - 1) - offset_in_seq;

        std::cout << "  i=" << i << " text_pos=" << text_pos << " char='" << (text_pos < concatenatedString.size() ? concatenatedString[text_pos] : '?')
                  << "' stored(seq_id=" << seq_id << ", seq_offset=" << seq_offset << ")"
                  << " expected_seq_id=" << expected_seq_id
                  << " offset_from_start=" << offset_in_seq << " offset_from_end=" << expected_offset_from_end
                  << " seq_len=" << seq_len
                  << " match_start?=" << (seq_offset == offset_in_seq ? "yes" : "no")
                  << " match_end?=" << (seq_offset == expected_offset_from_end ? "yes" : "no") << std::endl;

        if (seq_id != expected_seq_id) {
            matches_offset_from_end = false;
            matches_offset_from_start = false;
            break;
        }
        if (seq_offset != expected_offset_from_end) matches_offset_from_end = false;
        if (seq_offset != offset_in_seq) matches_offset_from_start = false;
    }

    std::cout << "[Samples_Offset_From_End] Result: matches_offset_from_end=" << (matches_offset_from_end ? "YES" : "no")
              << " matches_offset_from_start=" << (matches_offset_from_start ? "YES" : "no") << std::endl;

    ASSERT_TRUE(matches_offset_from_end || matches_offset_from_start)
        << "Stored sample offsets match neither offset-from-start nor offset-from-end (SA/seq_id mismatch?)";

    // r-index (as built/loaded) uses offset-from-START; assert that and document it
    ASSERT_TRUE(matches_offset_from_start)
        << "Head/tail samples do not match offset-from-start (SA/seq_id or ordering mismatch?)";
    std::cout << "[Samples_Offset_From_End] Assertion: r-index samples use offset-from-START (beginning of sequence)." << std::endl;

    // Verify every position: stored seq_offset should equal offset_from_start
    for (size_t i = 0; i < r_index_sa.size(); ++i) {
        size_t packed = r_index_sa[i];
        size_t seq_id = r_index.seqId(packed);
        size_t seq_offset = r_index.seqOffset(packed);

        size_t text_pos = expected_sa[i];
        size_t expected_seq_id = static_cast<size_t>(seq_indices[text_pos]);
        size_t offset_in_seq = text_pos - seq_start[expected_seq_id];

        ASSERT_EQ(seq_id, expected_seq_id) << "At BWT position i=" << i << ": seq_id mismatch";
        ASSERT_EQ(seq_offset, offset_in_seq)
            << "At BWT position i=" << i << ": seq_offset should be offset-from-START (got " << seq_offset
            << ", expected offset_from_start=" << offset_in_seq << ")";
    }

    size_t n_runs = r_index.tot_runs();
    std::cout << "[Samples_Offset_From_End] Head samples (tot_runs=" << n_runs << "):" << std::endl;
    for (size_t run_id = 0; run_id < n_runs; ++run_id) {
        size_t head_packed = r_index.getSample(run_id);
        size_t head_seq_id = r_index.seqId(head_packed);
        size_t head_offset = r_index.seqOffset(head_packed);

        std::cout << "  run_id=" << run_id << " head(seq_id=" << head_seq_id << ", seq_offset=" << head_offset << ")" << std::endl;

        ASSERT_LT(head_seq_id, static_cast<size_t>(n_seq)) << "Invalid head sample seq_id at run " << run_id;
        ASSERT_LT(head_offset, r_index.header.max_length) << "Invalid head sample offset at run " << run_id;
    }
    std::cout << "[Samples_Offset_From_End] Done." << std::endl;
}

TEST(FMDINDEX_Test, CompareSampledKmersNoN_Bidirectional_Encoded) {
    std::string rlbwt_file = "../test_data/bidirectional_test/contigs_xy.rl_bwt";
    FastLocate built(rlbwt_file);
    std::string enc_path = "./tmp_bidirectional.ri";
    {
        std::ofstream out(enc_path, std::ios::binary | std::ios::trunc);
        ASSERT_TRUE(out.good());
        built.serialize_encoded(out);
    }
    FastLocate r_index;
    {
        std::ifstream in(enc_path, std::ios::binary);
        ASSERT_TRUE(in.good());
        r_index.load_encoded(in);
    }
    r_index.initialize_complement_table();

    std::ifstream in_txt("../test_data/bidirectional_test/contigs_xy");
    ASSERT_TRUE(in_txt.is_open()) << "Could not open input text file";
    std::string full_text, line;
    while (std::getline(in_txt, line)) { if (!line.empty()) full_text += line; }
    ASSERT_GE(full_text.size(), 20) << "Input text is too short to sample from";

    const size_t k = 12;
    const size_t num_samples = 100;
    std::mt19937 rng(11);
    std::uniform_int_distribution<size_t> dist(0, full_text.size() - k);

    for (size_t i = 0; i < num_samples; ++i) {
        size_t pos = dist(rng);
        std::string kmer = full_text.substr(pos, k);
        if (kmer.find('N') != std::string::npos) { i--; continue; }

        std::string revcomp = kmer;
        std::reverse(revcomp.begin(), revcomp.end());
        for (char& c : revcomp) { c = r_index.complement(c); }

        panindexer::FastLocate::bi_interval int_kmer = {0, 0, r_index.bwt_size()};
        for (int j = kmer.size() - 1; j >= 0; --j) { int_kmer = r_index.backward_extend_encoded(int_kmer, kmer[j]); }

        panindexer::FastLocate::bi_interval int_rc = {0, 0, r_index.bwt_size()};
        for (int j = revcomp.size() - 1; j >= 0; --j) { int_rc = r_index.backward_extend_encoded(int_rc, revcomp[j]); }

        ASSERT_EQ(int_kmer.forward, int_rc.reverse) << "FMD symmetry violated between '" << kmer << "' and '" << revcomp << "'";
        ASSERT_EQ(int_kmer.reverse, int_rc.forward) << "FMD symmetry violated between '" << kmer << "' and '" << revcomp << "'";
        ASSERT_EQ(int_kmer.size, int_rc.size) << "FMD symmetry violated between '" << kmer << "' and '" << revcomp << "'";
    }
}

