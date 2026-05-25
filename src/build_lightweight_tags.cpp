// build_lightweight_tags:
//   Build a .ltags (LightTagIndex) file from an existing _compressed.tags
//   (compact-format TagArray) file. The .ltags file contains ONLY the
//   bwt_intervals bitvector — no tag values. Used by `find_mems --lightweight-tags`.
//
// Usage:
//   build_lightweight_tags <input_compressed.tags> <output.ltags> [--validate]
//
// --validate: after writing, re-load the .ltags and verify that its
//             bwt_intervals are bit-for-bit identical to the source.
//             Also performs a sample of rank/select consistency checks.

#include "pangenome_index/light_tag_index.hpp"
#include "pangenome_index/tag_arrays.hpp"

#include <sdsl/sd_vector.hpp>

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>

using namespace panindexer;

static void usage(const char* prog) {
    std::cerr << "Usage: " << prog
              << " <input_compressed.tags> <output.ltags> [--validate]\n";
}

static long long file_size_bytes(const std::string& path) {
    struct stat st{};
    if (::stat(path.c_str(), &st) != 0) return -1;
    return static_cast<long long>(st.st_size);
}

// Compare two sd_vectors by their populated 1-positions. Cheap enough for
// validation since we only need to do it once after build.
static bool sd_vectors_equal(const sdsl::sd_vector<>& a, const sdsl::sd_vector<>& b) {
    if (a.size() != b.size()) return false;
    sdsl::sd_vector<>::rank_1_type ra(&a);
    sdsl::sd_vector<>::rank_1_type rb(&b);
    if (ra(a.size()) != rb(b.size())) return false;
    const std::size_t n_ones = ra(a.size());
    sdsl::sd_vector<>::select_1_type sa(&a);
    sdsl::sd_vector<>::select_1_type sb(&b);
    for (std::size_t i = 1; i <= n_ones; ++i) {
        if (sa(i) != sb(i)) return false;
    }
    return true;
}

int main(int argc, char** argv) {
    if (argc < 3) { usage(argv[0]); return 2; }
    const std::string in_path  = argv[1];
    const std::string out_path = argv[2];
    bool validate = false;
    for (int i = 3; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--validate") { validate = true; }
        else if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
        else { std::cerr << "Unknown option: " << a << "\n"; usage(argv[0]); return 2; }
    }

    std::cerr << "========================================\n"
              << "build_lightweight_tags\n"
              << "========================================\n"
              << "  input:    " << in_path  << "\n"
              << "  output:   " << out_path << "\n"
              << "  validate: " << (validate ? "yes" : "no") << "\n"
              << "----------------------------------------\n";

    // --- Build ---
    LightTagIndex idx;
    try {
        auto t0 = std::chrono::high_resolution_clock::now();
        idx.build_from_compact_tags_file(in_path);
        auto t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        std::cerr << "Loaded compact tags from " << in_path << "\n"
                  << "  bwt_size:       " << idx.bwt_size() << "\n"
                  << "  num_runs:       " << idx.num_runs() << "\n"
                  << "  mem footprint:  " << idx.memory_bytes() << " bytes\n"
                  << "  build time:     " << dt.count() << " s\n";
    } catch (const std::exception& e) {
        std::cerr << "ERROR building light tag index: " << e.what() << "\n";
        return 1;
    }

    // --- Save ---
    try {
        std::ofstream out(out_path, std::ios::binary | std::ios::trunc);
        if (!out) {
            std::cerr << "ERROR: cannot open output file: " << out_path << "\n";
            return 1;
        }
        idx.save(out);
        out.close();
    } catch (const std::exception& e) {
        std::cerr << "ERROR writing " << out_path << ": " << e.what() << "\n";
        return 1;
    }

    long long in_sz  = file_size_bytes(in_path);
    long long out_sz = file_size_bytes(out_path);
    std::cerr << "Wrote " << out_path << " (" << out_sz << " bytes)\n";
    if (in_sz > 0 && out_sz > 0) {
        double ratio = static_cast<double>(in_sz) / static_cast<double>(out_sz);
        std::cerr << "  input  size:  " << in_sz  << " bytes\n"
                  << "  output size:  " << out_sz << " bytes\n"
                  << "  compression:  " << ratio << "x smaller than input\n";
    }

    // --- Validate ---
    if (validate) {
        std::cerr << "----------------------------------------\n"
                  << "Validating round-trip ...\n";
        LightTagIndex round_trip;
        try {
            std::ifstream in(out_path, std::ios::binary);
            round_trip.load(in);
        } catch (const std::exception& e) {
            std::cerr << "VALIDATE FAIL: cannot reload " << out_path
                      << ": " << e.what() << "\n";
            return 3;
        }
        if (round_trip.bwt_size() != idx.bwt_size()) {
            std::cerr << "VALIDATE FAIL: bwt_size mismatch ("
                      << round_trip.bwt_size() << " vs " << idx.bwt_size() << ")\n";
            return 3;
        }
        if (round_trip.num_runs() != idx.num_runs()) {
            std::cerr << "VALIDATE FAIL: num_runs mismatch ("
                      << round_trip.num_runs() << " vs " << idx.num_runs() << ")\n";
            return 3;
        }
        if (!sd_vectors_equal(round_trip.get_bwt_intervals(), idx.get_bwt_intervals())) {
            std::cerr << "VALIDATE FAIL: sd_vector contents differ\n";
            return 3;
        }
        std::cerr << "  round-trip:                OK ("
                  << round_trip.num_runs() << " run starts match)\n";

        // Cross-check against the source TagArray. This catches the case
        // where build_from_compact_tags_file silently produces an empty or
        // truncated bwt_intervals.
        TagArray src;
        try {
            std::ifstream sin(in_path, std::ios::binary);
            src.load_compressed_tags_compact(sin);
        } catch (const std::exception& e) {
            std::cerr << "VALIDATE WARN: cannot reload source for cross-check: "
                      << e.what() << "\n";
            return 0; // round-trip succeeded; cross-check is best-effort.
        }
        if (!sd_vectors_equal(round_trip.get_bwt_intervals(), src.get_bwt_intervals())) {
            std::cerr << "VALIDATE FAIL: round-tripped index differs from source TagArray's bwt_intervals\n";
            return 3;
        }
        std::cerr << "  cross-check vs source:     OK\n";

        // Sample rank/select consistency: for K positions spaced uniformly across
        // the bitvector, verify select(rank(p)+0) == start of the run containing p.
        const std::size_t N = idx.bwt_size();
        const std::size_t K = std::min<std::size_t>(1000, N);
        std::size_t checked = 0;
        for (std::size_t i = 1; i <= K; ++i) {
            std::size_t p   = (N * i) / (K + 1);
            std::size_t rid = round_trip.run_id_at(p);
            std::size_t s   = round_trip.run_start_bwt(rid);
            if (s > p) {
                std::cerr << "VALIDATE FAIL: select(rank(p)) > p at p=" << p
                          << " rid=" << rid << " s=" << s << "\n";
                return 3;
            }
            // The next run must start strictly after p (or rid is the last run).
            if (rid + 1 < round_trip.num_runs()) {
                std::size_t s_next = round_trip.run_start_bwt(rid + 1);
                if (s_next <= p) {
                    std::cerr << "VALIDATE FAIL: select(rank(p)+1) <= p at p=" << p
                              << " rid=" << rid << " s_next=" << s_next << "\n";
                    return 3;
                }
            }
            ++checked;
        }
        std::cerr << "  rank/select consistency:   OK (" << checked
                  << " random positions sampled)\n";
        std::cerr << "VALIDATE: PASS\n";
    }

    return 0;
}
