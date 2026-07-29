// Unit test for FastLocate::leftmost_c_in_interval and the per-character
// run-start bitvector infrastructure added in commit 3de0ddf.
//
// Verifies the primitive against a brute-force scan using bwt_char_at_encoded
// on the same interval. Any mismatch is a bug; test exits non-zero.
//
// Usage: bin/test_r_index_sa <r_index_file>
//
// IMPORTANT: use trusted r-index files only (e.g. the yeast chrII index at
// runs/v1-current/yeast235_chrII_100kb_normalized.ri, or HPRCv1 chr6 on
// vesuvio). Do NOT use xy-test/xy.ri -- it is a known-bad fixture.
//
// Rationale: this primitive is used by the SA-carrying rule in the flipped
// MEM finder (see DESIGN_FLIPPED_MEM.md section 4). Correctness is
// non-negotiable -- wrong results here would cascade into wrong MEM
// enumerations, wrong tag emissions, wrong biological output.

#include <pangenome_index/r-index.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using panindexer::FastLocate;

namespace {

// The 5 DNA characters the primitive supports.
constexpr char DNA_CHARS[5] = {'A', 'C', 'G', 'T', 'N'};

// Brute-force ground truth: scan [sp, ep] for first position with BWT == c.
// O(ep - sp + 1) per query. Since bwt_char_at_encoded is itself a predecessor
// query + block scan (not O(1)), the caller must keep intervals small.
size_t leftmost_c_naive(const FastLocate& idx, size_t c, size_t sp, size_t ep) {
    if (sp > ep) return SIZE_MAX;
    if (ep >= idx.sequence_size) return SIZE_MAX;  // out of range
    for (size_t p = sp; p <= ep; ++p) {
        if (idx.bwt_char_at_encoded(p) == c) {
            return p;
        }
    }
    return SIZE_MAX;
}

// Maximum interval size for which we invoke the naive baseline. Larger
// intervals would make the naive scan the bottleneck. Chosen so a single
// naive scan runs in well under 1 ms on typical hardware.
constexpr size_t MAX_NAIVE_INTERVAL = 4096;

struct TestStats {
    size_t total = 0;
    size_t passed = 0;
    size_t failed = 0;
    std::vector<std::string> failure_details;  // limit to first 10

    void record(bool ok, const std::string& detail) {
        total++;
        if (ok) {
            passed++;
        } else {
            failed++;
            if (failure_details.size() < 10) {
                failure_details.push_back(detail);
            }
        }
    }
};

// Run one test case and record the result.
void check_one(FastLocate& idx, size_t c, size_t sp, size_t ep, TestStats& stats,
               const std::string& label) {
    size_t got = idx.leftmost_c_in_interval(c, sp, ep);
    size_t expected = leftmost_c_naive(idx, c, sp, ep);
    bool ok = (got == expected);
    if (!ok) {
        std::string detail = label + " c='" + std::string(1, (char)c) +
                             "' sp=" + std::to_string(sp) +
                             " ep=" + std::to_string(ep) +
                             " got=" + (got == SIZE_MAX ? "SIZE_MAX" : std::to_string(got)) +
                             " expected=" + (expected == SIZE_MAX ? "SIZE_MAX" : std::to_string(expected));
        stats.record(false, detail);
    } else {
        stats.record(true, "");
    }
}

void run_edge_cases(FastLocate& idx, TestStats& stats) {
    const size_t n = idx.sequence_size;
    std::cerr << "  edge cases (bwt_size=" << n << ")..." << std::endl;

    for (char c_char : DNA_CHARS) {
        size_t c = static_cast<size_t>(c_char);

        // Single-position interval at position 0.
        check_one(idx, c, 0, 0, stats, "edge:sp=ep=0");

        // Single-position interval at last valid position.
        if (n > 0) {
            check_one(idx, c, n - 1, n - 1, stats, "edge:sp=ep=last");
        }

        // Whole-BWT interval (only when tractable for the naive baseline).
        if (n > 0 && n <= MAX_NAIVE_INTERVAL) {
            check_one(idx, c, 0, n - 1, stats, "edge:whole");
        }

        // Invalid interval: sp > ep. Primitive must return SIZE_MAX.
        check_one(idx, c, 5, 3, stats, "edge:sp>ep");

        // Sample point-intervals across the BWT. Adaptive stride so this
        // stays fast on large indexes -- aim for ~1000 samples per character.
        size_t stride = std::max<size_t>(1, n / 1000);
        for (size_t p = 0; p < n; p += stride) {
            check_one(idx, c, p, p, stats, "edge:point");
        }
    }
}

void run_random_tests(FastLocate& idx, size_t n_tests, uint32_t seed, TestStats& stats) {
    const size_t n = idx.sequence_size;
    if (n == 0) {
        std::cerr << "  skipping random tests (empty BWT)" << std::endl;
        return;
    }
    std::cerr << "  " << n_tests << " random (c, sp, ep) tests (seed=" << seed << ")..." << std::endl;

    std::mt19937 rng(seed);
    std::uniform_int_distribution<size_t> char_dist(0, 4);
    std::uniform_int_distribution<size_t> pos_dist(0, n - 1);
    // Cap interval size so the naive baseline stays cheap.
    // Distribution: mix of small (typical) and larger (up to cap) intervals.
    std::uniform_int_distribution<size_t> len_dist(1, MAX_NAIVE_INTERVAL);

    for (size_t i = 0; i < n_tests; ++i) {
        size_t c = static_cast<size_t>(DNA_CHARS[char_dist(rng)]);
        size_t sp = pos_dist(rng);
        size_t len = len_dist(rng);
        size_t ep = std::min(sp + len - 1, n - 1);
        check_one(idx, c, sp, ep, stats, "random#" + std::to_string(i));
    }
}

// Verify structural invariant: every position marked in run_starts_by_char[c]
// has BWT[pos] == c. This tests the build correctness of the bitvectors
// themselves, independent of the query logic.
void verify_marker_positions(FastLocate& idx, TestStats& stats) {
    std::cerr << "  verifying marker positions match BWT..." << std::endl;
    // Force the bitvectors to be built by calling the query once.
    (void)idx.leftmost_c_in_interval(static_cast<size_t>('A'), 0, 0);

    // Iterate marked positions of each character and verify BWT[pos] == c.
    // We access the bitvectors through the public API only via
    // leftmost_c_in_interval, so we instead sample: iterate every position in
    // the BWT and count runs by character, comparing against the primitive's
    // claim of "leftmost c in [pos, pos]".
    // This is O(n) but only runs once and gives a strong structural check.
    const size_t n = idx.sequence_size;
    // Skip if huge; for xy.ri this should be small.
    if (n > 5'000'000) {
        std::cerr << "    (skipping full-BWT structural check: n=" << n
                  << " too large)" << std::endl;
        return;
    }
    size_t char_mismatches = 0;
    for (size_t p = 0; p < n; ++p) {
        size_t actual_c = idx.bwt_char_at_encoded(p);
        // For each of the 5 DNA chars, verify:
        //   leftmost_c_in_interval(c, p, p) == p iff actual_c == c
        for (char c_char : DNA_CHARS) {
            size_t c = static_cast<size_t>(c_char);
            size_t got = idx.leftmost_c_in_interval(c, p, p);
            bool expects_hit = (actual_c == c);
            bool got_hit = (got == p);
            if (expects_hit != got_hit) {
                if (char_mismatches < 10) {
                    std::cerr << "    MISMATCH: pos=" << p << " actual_c='"
                              << (char)actual_c << "' testing_c='" << c_char
                              << "' got=" << (got == SIZE_MAX ? "SIZE_MAX" : std::to_string(got))
                              << " expected_hit=" << expects_hit << std::endl;
                }
                char_mismatches++;
                stats.record(false, "struct:pos=" + std::to_string(p) + " c='" + c_char + "'");
            } else {
                stats.record(true, "");
            }
        }
    }
    if (char_mismatches == 0) {
        std::cerr << "    OK: " << (n * 5) << " point queries all consistent" << std::endl;
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <r_index_file> [num_random_tests=10000] [seed=42]"
                  << std::endl;
        return 2;
    }
    const std::string ri_path = argv[1];
    const size_t n_random = (argc >= 3) ? std::stoul(argv[2]) : 10000;
    const uint32_t seed = (argc >= 4) ? static_cast<uint32_t>(std::stoul(argv[3])) : 42;

    std::cerr << "Loading r-index from: " << ri_path << std::endl;
    FastLocate idx;
    std::ifstream in(ri_path, std::ios::binary);
    if (!in) {
        std::cerr << "ERROR: cannot open " << ri_path << std::endl;
        return 2;
    }
    idx.load_encoded(in);
    std::cerr << "Loaded. bwt_size=" << idx.sequence_size
              << " num_runs=" << idx.size() << std::endl;

    TestStats stats;

    std::cerr << "\n=== Structural verification ===" << std::endl;
    verify_marker_positions(idx, stats);

    std::cerr << "\n=== Edge cases ===" << std::endl;
    run_edge_cases(idx, stats);

    std::cerr << "\n=== Random tests ===" << std::endl;
    run_random_tests(idx, n_random, seed, stats);

    std::cerr << "\n=== Summary ===" << std::endl;
    std::cerr << "Total:  " << stats.total << std::endl;
    std::cerr << "Passed: " << stats.passed << std::endl;
    std::cerr << "Failed: " << stats.failed << std::endl;

    if (stats.failed > 0) {
        std::cerr << "\nFirst failures:" << std::endl;
        for (const auto& d : stats.failure_details) {
            std::cerr << "  " << d << std::endl;
        }
        std::cerr << "\nFAIL" << std::endl;
        return 1;
    }
    std::cerr << "\nPASS" << std::endl;
    return 0;
}
