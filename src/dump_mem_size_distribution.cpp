// One-off diagnostic to characterize the MEM-size distribution divergence
// between the legacy 3-phase finder and the flipped SA-carry finder,
// specifically to understand Risk E in DESIGN_FLIPPED_MEM.md sec 5.3.
//
// Emits one TSV row per MEM (from both finders), each row tagged with the
// finder that produced it plus a classification of whether that MEM is
// left-maximal (via count_encoded on the left-extension). Downstream awk
// pipelines produce the actual histograms.
//
// TSV schema (writes to stdout, header on first line):
//   finder   read_id   read_len   mem_start   mem_end   mem_len   mem_size
//     left_maximal   in_other_finder
//
//   finder            "legacy" | "flipped"
//   left_maximal      "yes" | "no" (yes if start==0 or count_encoded of
//                      pattern[start-1..end) is < min_occ)
//   in_other_finder   "yes" | "no" (set-membership check on (start,end,bwt_start,size))
//
// Usage:
//   bin/dump_mem_size_distribution <r_index_file> <reads_file>
//                                   [min_len=50] [min_occ=1] [max_reads=0=all]
//
// See VALIDATION_GUIDE.md for the canonical dataset paths to feed this.
//
// Not a pass/fail test; there is no exit code contract beyond
//   0 = ran to completion; nonzero = argv/io error.

#include <pangenome_index/algorithm.hpp>
#include <pangenome_index/r-index.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>

using panindexer::FastLocate;
using panindexer::MEM;
using panindexer::MEM_with_sa;

namespace {

// Identical semantics to test_flipped_mems.cpp's is_left_maximal: uses
// count_encoded, which is O(|substring|) backward-extends. Only called on
// emitted MEMs so total cost stays proportional to the MEM set, not the
// text length.
bool is_left_maximal(FastLocate& idx, const std::string& pattern,
                     size_t start, size_t end, size_t min_occ) {
    if (start == 0) return true;
    std::string ext = pattern.substr(start - 1, end - start + 1);
    auto range = idx.count_encoded(ext);
    if (range.first > range.second) return true;
    size_t occ = range.second - range.first + 1;
    return occ < min_occ;
}

using MemKey = std::tuple<size_t, size_t, size_t, int64_t>;

struct MemKeyHash {
    size_t operator()(const MemKey& k) const noexcept {
        // Splice the four fields into one 64-bit hash. Collisions here just
        // fall back to std::unordered_set's tuple equality; the goal is
        // distribution, not cryptography.
        auto h = std::hash<size_t>{};
        return h(std::get<0>(k)) ^ (h(std::get<1>(k)) << 1)
             ^ (h(std::get<2>(k)) << 2) ^ (h(static_cast<size_t>(std::get<3>(k))) << 3);
    }
};

MemKey key_of(const MEM& m)         { return {m.start, m.end, m.bwt_start, m.size}; }
MemKey key_of(const MEM_with_sa& m) { return {m.start, m.end, m.bwt_start, m.size}; }

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <r_index_file> <reads_file>"
                     " [min_len=50] [min_occ=1] [max_reads=0=all]" << std::endl;
        return 2;
    }
    const std::string ri_path = argv[1];
    const std::string reads_path = argv[2];
    const size_t min_len = (argc >= 4) ? std::stoul(argv[3]) : 50;
    const size_t min_occ = (argc >= 5) ? std::stoul(argv[4]) : 1;
    const size_t max_reads = (argc >= 6) ? std::stoul(argv[5]) : 0;

    std::cerr << "Loading r-index from: " << ri_path << std::endl;
    FastLocate idx;
    {
        std::ifstream in(ri_path, std::ios::binary);
        if (!in) { std::cerr << "ERROR: cannot open " << ri_path << std::endl; return 2; }
        idx.load_encoded(in);
    }
    std::cerr << "bwt_size=" << idx.bwt_size() << " num_runs=" << idx.size() << std::endl;

    std::ifstream reads(reads_path);
    if (!reads) { std::cerr << "ERROR: cannot open " << reads_path << std::endl; return 2; }
    std::cerr << "min_len=" << min_len << " min_occ=" << min_occ
              << " max_reads=" << (max_reads == 0 ? "all" : std::to_string(max_reads))
              << std::endl;

    // TSV header on stdout so downstream awk / pandas can consume it.
    std::cout << "finder\tread_id\tread_len\tmem_start\tmem_end\tmem_len\tmem_size"
                 "\tleft_maximal\tin_other_finder\n";

    std::string read;
    int read_id = 0;
    size_t processed = 0;
    while (std::getline(reads, read)) {
        if (read.empty()) continue;
        read_id++;
        if (max_reads > 0 && processed >= max_reads) break;
        processed++;

        auto legacy = panindexer::find_all_mems(read, min_len, min_occ, idx);
        auto flipped = panindexer::find_all_mems_flipped(read, min_len, min_occ, idx);

        std::unordered_set<MemKey, MemKeyHash> legacy_keys, flipped_keys;
        legacy_keys.reserve(legacy.size());
        flipped_keys.reserve(flipped.size());
        for (const auto& m : legacy)  legacy_keys.insert(key_of(m));
        for (const auto& m : flipped) flipped_keys.insert(key_of(m));

        auto emit = [&](const char* finder, size_t s, size_t e, size_t bwt_start,
                        int64_t size, bool in_other) {
            bool lm = is_left_maximal(idx, read, s, e, min_occ);
            std::cout << finder << '\t' << read_id << '\t' << read.size()
                      << '\t' << s << '\t' << e << '\t' << (e - s)
                      << '\t' << size
                      << '\t' << (lm ? "yes" : "no")
                      << '\t' << (in_other ? "yes" : "no") << '\n';
        };

        for (const auto& m : legacy) {
            emit("legacy", m.start, m.end, m.bwt_start, m.size,
                 flipped_keys.count(key_of(m)) > 0);
        }
        for (const auto& m : flipped) {
            emit("flipped", m.start, m.end, m.bwt_start, m.size,
                 legacy_keys.count(key_of(m)) > 0);
        }

        if (processed % 1000 == 0) {
            std::cerr << "  processed " << processed << " reads" << std::endl;
        }
    }

    std::cerr << "done: " << processed << " reads processed" << std::endl;
    return 0;
}
