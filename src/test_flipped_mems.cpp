// Unit test for find_all_mems_flipped: the right-to-left, pure-backward-extend
// MEM finder that carries SA[sp] through Step 1' (see DESIGN_FLIPPED_MEM.md
// section 5.2 Stage 2).
//
// Verification bar (per the design doc):
//   1. Bitwise-identical (start, end, bwt_start, size) MEM sets vs the legacy
//      find_all_mems, per read. MEMs are compared as sorted sets because
//      emission order differs (legacy: left-to-right; flipped: right-to-left).
//   2. For every emitted MEM, the flipped path's sa_sp equals the ground-truth
//      locate_sa_value(bwt_start) computed via the legacy locate walk. This is
//      the algorithm-level correctness check that ties Stage 2's output back
//      to the Stage 1 primitive-level guarantee.
//
// Any divergence is a hard fail. The MEM set is a downstream input to
// find_mems -> gafpack -> validate_gaf; a wrong MEM here would corrupt
// biology exactly the way VALIDATION_GUIDE.md warns against.
//
// Usage:
//   bin/test_flipped_mems <r_index_file> <reads_file>
//                         [min_len=30] [min_occ=1] [max_reads=0=all]
//
// IMPORTANT: use the canonical validation datasets listed in
// VALIDATION_GUIDE.md (do NOT use xy-test/xy.ri -- known-bad fixture, and do
// NOT use shipped subsets like xy-test/yeast_500reads.txt -- iterate on the
// canonical file with max_reads instead). Recommended:
//   .ri:    ../mem-projection/pangenome-pipeline/runs/v2-yeast235/
//           yeast235_chrII_100kb_normalized.ri
//   reads:  ../mem-projection/yeast-235/yeast-235-chrI/
//           S288C_chrII_N100K_R1_200_reads.txt
//
// For quick iteration, pass max_reads=500 (or any N) instead of splitting the
// reads file into a separate on-disk sample; that keeps the sample
// deterministic (the file's first N lines) and rooted in the canonical
// dataset the guide references.

#include <pangenome_index/algorithm.hpp>
#include <pangenome_index/r-index.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using panindexer::FastLocate;
using panindexer::MEM;
using panindexer::MEM_with_sa;
using size_type = FastLocate::size_type;

namespace {

// Ground truth locate: same walk as find_mems.cpp:locate_sa_value. Kept
// private to this test so it stays a self-contained oracle.
size_type locate_sa_value_naive(const FastLocate& idx, size_t bwt_pos) {
    size_t run_id = 0;
    size_t offset_of_first = 0;
    idx.run_id_and_offset_at(bwt_pos, run_id, offset_of_first);
    size_type sa = idx.getSample(run_id);
    while (offset_of_first < bwt_pos) {
        sa = idx.locateNext(sa);
        offset_of_first++;
    }
    return sa;
}

// Sortable key over the four fields that must match byte-for-byte between the
// two finders (per the Stage 2 bar). sa_sp is deliberately excluded from
// equality -- it's compared separately against locate_sa_value ground truth.
using MemKey = std::tuple<size_t, size_t, size_t, int64_t>;

MemKey key_of(const MEM& m) {
    return std::make_tuple(m.start, m.end, m.bwt_start, m.size);
}
MemKey key_of(const MEM_with_sa& m) {
    return std::make_tuple(m.start, m.end, m.bwt_start, m.size);
}

std::string format_mem(const MemKey& k) {
    std::ostringstream os;
    os << "(start=" << std::get<0>(k)
       << ", end=" << std::get<1>(k)
       << ", bwt_start=" << std::get<2>(k)
       << ", size=" << std::get<3>(k) << ")";
    return os.str();
}

struct Stats {
    size_t reads_processed = 0;
    size_t legacy_mems_total = 0;
    size_t flipped_mems_total = 0;
    // Distinguish TWO kinds of set-difference cases (see check_read):
    //   * flipped_missed:  MEM in legacy \ flipped that is actually left-maximal
    //                      -> flipped bug, hard fail
    //   * legacy_spurious: MEM in legacy \ flipped that is NOT left-maximal
    //                      -> legacy bug (Risk E in DESIGN_FLIPPED_MEM.md), tolerated
    //   * flipped_extra:   MEM in flipped \ legacy -> flipped bug, hard fail
    size_t reads_with_flipped_missed = 0;
    size_t reads_with_flipped_extra = 0;
    size_t reads_with_legacy_spurious = 0;
    size_t total_flipped_missed = 0;
    size_t total_flipped_extra = 0;
    size_t total_legacy_spurious = 0;
    size_t reads_with_sa_mismatch = 0;
    size_t sa_mismatches = 0;
    size_t sa_checks = 0;
    // Cap detailed reporting so log output stays bounded on large inputs.
    std::vector<std::string> failure_details;
    std::vector<std::string> spurious_examples;  // for informational context

    void record_fail(const std::string& d) {
        if (failure_details.size() < 20) failure_details.push_back(d);
    }
    void record_spurious(const std::string& d) {
        if (spurious_examples.size() < 10) spurious_examples.push_back(d);
    }
};

// Returns true iff the substring pattern[start..end) is left-maximal, meaning
// either it starts at position 0 or its left-extension pattern[start-1..end)
// has strictly fewer than min_occ occurrences in the index.
//
// Uses the r-index count() primitive, which is O(|substring|) backward-extends
// -- affordable here because we only call this on MEMs that already diverged
// between the two finders (a small subset in practice).
bool is_left_maximal(FastLocate& idx, const std::string& pattern,
                     size_t start, size_t end, size_t min_occ) {
    if (start == 0) return true;
    // count_encoded takes a non-const string& (a signature quirk of the
    // legacy count()/count_encoded() pair). Materialize the left-extended
    // substring locally so we can pass an lvalue reference.
    std::string ext = pattern.substr(start - 1, end - start + 1);
    auto range = idx.count_encoded(ext);
    // count_encoded returns {first, second}; empty range means first > second.
    if (range.first > range.second) return true;  // 0 occurrences
    size_t occ = range.second - range.first + 1;
    return occ < min_occ;
}

// Compare one read. Returns true iff both the MEM sets match and every
// sa_sp checks out.
bool check_read(FastLocate& idx, const std::string& read, int read_id,
                size_t min_len, size_t min_occ, Stats& stats) {
    auto legacy = panindexer::find_all_mems(read, min_len, min_occ, idx);
    auto flipped = panindexer::find_all_mems_flipped(read, min_len, min_occ, idx);
    stats.legacy_mems_total += legacy.size();
    stats.flipped_mems_total += flipped.size();

    // --- Set equality check on the four legacy MEM fields ---------------
    std::vector<MemKey> lkeys, fkeys;
    lkeys.reserve(legacy.size());
    fkeys.reserve(flipped.size());
    for (const auto& m : legacy) lkeys.push_back(key_of(m));
    for (const auto& m : flipped) fkeys.push_back(key_of(m));
    std::sort(lkeys.begin(), lkeys.end());
    std::sort(fkeys.begin(), fkeys.end());

    // --- Diff the two sets ----------------------------------------------
    // Note on classification below: legacy has a known Step 3 bug (Risk E in
    // DESIGN_FLIPPED_MEM.md sec 5.3) where extending to pattern[len] reads a
    // null byte that the encoded index resolves to symbol code 0 -- which
    // aliases to the first DNA character. That makes Step 3 continue instead
    // of terminating, and the outer loop then advances by 1 anchor per
    // iteration, emitting one spurious non-left-maximal MEM per position for
    // every anchor whose true MEM extends to end-of-pattern. Those spurious
    // emissions are what this test tolerates below -- confirmed by checking
    // left-maximality via count_encoded on the left-extension.
    std::vector<MemKey> only_legacy, only_flipped;
    std::set_difference(lkeys.begin(), lkeys.end(), fkeys.begin(), fkeys.end(),
                        std::back_inserter(only_legacy));
    std::set_difference(fkeys.begin(), fkeys.end(), lkeys.begin(), lkeys.end(),
                        std::back_inserter(only_flipped));

    size_t local_flipped_missed = 0;
    size_t local_legacy_spurious = 0;
    for (const auto& k : only_legacy) {
        size_t s = std::get<0>(k), e = std::get<1>(k);
        if (is_left_maximal(idx, read, s, e, min_occ)) {
            // Legacy emitted a MEM that IS left-maximal but flipped missed
            // it. This is a genuine flipped bug.
            local_flipped_missed++;
            std::ostringstream os;
            os << "read#" << read_id << " FLIPPED MISSED left-maximal MEM: "
               << format_mem(k);
            stats.record_fail(os.str());
        } else {
            // Legacy emitted a non-left-maximal MEM -> spurious per definition.
            local_legacy_spurious++;
            if (stats.total_legacy_spurious < 5) {
                std::ostringstream os;
                os << "read#" << read_id << " legacy spurious (not left-maximal): "
                   << format_mem(k);
                stats.record_spurious(os.str());
            }
        }
    }

    if (local_flipped_missed) stats.reads_with_flipped_missed++;
    if (local_legacy_spurious) stats.reads_with_legacy_spurious++;
    stats.total_flipped_missed += local_flipped_missed;
    stats.total_legacy_spurious += local_legacy_spurious;

    if (!only_flipped.empty()) {
        // Any MEM flipped emits that legacy doesn't is a hard failure --
        // flipped inventing MEMs from thin air would corrupt downstream.
        stats.reads_with_flipped_extra++;
        stats.total_flipped_extra += only_flipped.size();
        std::ostringstream os;
        os << "read#" << read_id << " FLIPPED EXTRA (not in legacy): "
           << format_mem(only_flipped.front())
           << " (+" << (only_flipped.size() - 1) << " more)";
        stats.record_fail(os.str());
    }

    bool set_ok = (local_flipped_missed == 0) && only_flipped.empty();

    // --- sa_sp check on the flipped emissions ---------------------------
    // Even when the MEM set matches, sa_sp is still checked independently
    // -- getting the interval right but the SA wrong would silently drop the
    // primary win of this design.
    bool sa_ok = true;
    for (const auto& m : flipped) {
        stats.sa_checks++;
        size_type gt = locate_sa_value_naive(idx, m.bwt_start);
        if (m.sa_sp != gt) {
            sa_ok = false;
            stats.sa_mismatches++;
            std::ostringstream os;
            os << "read#" << read_id << " SA MISMATCH mem=" << format_mem(key_of(m))
               << " sa_sp_flipped=" << m.sa_sp << " sa_ground_truth=" << gt;
            stats.record_fail(os.str());
        }
    }
    if (!sa_ok) stats.reads_with_sa_mismatch++;

    return set_ok && sa_ok;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <r_index_file> <reads_file>"
                     " [min_len=30] [min_occ=1] [max_reads=0=all]" << std::endl;
        return 2;
    }
    const std::string ri_path = argv[1];
    const std::string reads_path = argv[2];
    const size_t min_len = (argc >= 4) ? std::stoul(argv[3]) : 30;
    const size_t min_occ = (argc >= 5) ? std::stoul(argv[4]) : 1;
    const size_t max_reads = (argc >= 6) ? std::stoul(argv[5]) : 0;  // 0 = all

    std::cerr << "Loading r-index from: " << ri_path << std::endl;
    FastLocate idx;
    {
        std::ifstream in(ri_path, std::ios::binary);
        if (!in) {
            std::cerr << "ERROR: cannot open " << ri_path << std::endl;
            return 2;
        }
        idx.load_encoded(in);
    }
    std::cerr << "Loaded. bwt_size=" << idx.bwt_size()
              << " num_runs=" << idx.size() << std::endl;

    std::ifstream reads(reads_path);
    if (!reads) {
        std::cerr << "ERROR: cannot open " << reads_path << std::endl;
        return 2;
    }
    std::cerr << "min_len=" << min_len << " min_occ=" << min_occ
              << " max_reads=" << (max_reads == 0 ? "all" : std::to_string(max_reads))
              << std::endl;

    Stats stats;
    std::string read;
    int read_id = 0;
    while (std::getline(reads, read)) {
        if (read.empty()) continue;
        read_id++;
        if (max_reads > 0 && stats.reads_processed >= max_reads) break;
        stats.reads_processed++;
        check_read(idx, read, read_id, min_len, min_occ, stats);
        // Cheap progress ping every 100 reads so long runs are observable.
        if (stats.reads_processed % 100 == 0) {
            std::cerr << "  processed " << stats.reads_processed
                      << " reads (legacy=" << stats.legacy_mems_total
                      << " flipped=" << stats.flipped_mems_total
                      << " flipped_missed=" << stats.total_flipped_missed
                      << " flipped_extra=" << stats.total_flipped_extra
                      << " legacy_spurious=" << stats.total_legacy_spurious
                      << " sa_mismatches=" << stats.sa_mismatches << ")" << std::endl;
        }
    }

    std::cerr << "\n=== Summary ===" << std::endl;
    std::cerr << "Reads processed:                    " << stats.reads_processed << std::endl;
    std::cerr << "Legacy MEMs emitted total:          " << stats.legacy_mems_total << std::endl;
    std::cerr << "Flipped MEMs emitted total:         " << stats.flipped_mems_total << std::endl;
    std::cerr << "Reads with FLIPPED_MISSED (bug):    " << stats.reads_with_flipped_missed << std::endl;
    std::cerr << "Reads with FLIPPED_EXTRA  (bug):    " << stats.reads_with_flipped_extra << std::endl;
    std::cerr << "Reads with LEGACY_SPURIOUS (tolerated): " << stats.reads_with_legacy_spurious << std::endl;
    std::cerr << "Total MEMs flipped-missed:          " << stats.total_flipped_missed << std::endl;
    std::cerr << "Total MEMs flipped-extra:           " << stats.total_flipped_extra << std::endl;
    std::cerr << "Total MEMs legacy-spurious:         " << stats.total_legacy_spurious << std::endl;
    std::cerr << "Reads with SA mismatch:             " << stats.reads_with_sa_mismatch << std::endl;
    std::cerr << "Total sa_sp checks:                 " << stats.sa_checks << std::endl;
    std::cerr << "Total sa_sp mismatches:             " << stats.sa_mismatches << std::endl;

    if (!stats.failure_details.empty()) {
        std::cerr << "\nFailure details:" << std::endl;
        for (const auto& d : stats.failure_details) std::cerr << "  " << d << std::endl;
    }
    if (!stats.spurious_examples.empty()) {
        std::cerr << "\nLegacy-spurious examples (informational, not failures):" << std::endl;
        for (const auto& d : stats.spurious_examples) std::cerr << "  " << d << std::endl;
    }

    // Fail only on the two conditions that indicate a bug in the NEW code:
    //   * flipped missed a left-maximal MEM
    //   * flipped emitted a MEM legacy didn't
    //   * any sa_sp disagreed with locate_sa_value ground truth
    // Legacy-spurious MEMs are noted but tolerated (documented legacy bug).
    if (stats.total_flipped_missed > 0 || stats.total_flipped_extra > 0 ||
        stats.sa_mismatches > 0) {
        std::cerr << "\nFAIL" << std::endl;
        return 1;
    }
    if (stats.reads_processed == 0) {
        std::cerr << "\nERROR: no reads processed" << std::endl;
        return 1;
    }
    if (stats.flipped_mems_total == 0) {
        std::cerr << "\nERROR: zero MEMs emitted -- test is not verifying anything" << std::endl;
        return 1;
    }
    std::cerr << "\nPASS" << std::endl;
    return 0;
}
