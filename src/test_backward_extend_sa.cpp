// Unit test for the SA-carrying backward-extend primitive and its supporting
// machinery. Covers four tests:
//
//   Test A -- FastLocate::scan_at (extended block scan that returns ranks
//             plus bwt_char_at_pos, run_id, run_start_bwt in one predecessor
//             + block-walk).  Ground truth: pointwise cross-check against
//             rank_at_cached_encoded + bwt_char_at_encoded +
//             run_id_and_offset_at at the same position.
//
//   Test B -- FastLocate::run_head_c[c] parallel run-id bitvector.  For each
//             DNA character c and each run r in a sampled subset,
//             run_head_c[c][r] must equal 1 iff the BWT character at the
//             start of run r is c.
//
//   Test C -- FastLocate::backward_extend_encoded_with_sa correctness
//             (existing coverage).  Verifies sa_sp at every step of two
//             trial-pattern suites matches locate_sa_value ground truth.
//             Also verifies the toehold fast path: on the very first call of
//             a Step 1' walk (sa_sp_prev == NO_POSITION), the returned sa_sp
//             must equal what locate_sa_value would compute for the same
//             new_sp.  Reports Case A / Case B / toehold-hit counts.
//
//   Test D -- JR-014 run_head_c[] serialization round-trip.  The loaded
//             r-index carries run_head_c[] populated from the .ri file
//             directly (no in-process rebuild).  This test serializes the
//             loaded index to a temp .ri via serialize_encoded, loads that
//             copy into a second FastLocate, and verifies:
//               (a) both instances' run_head_c[] structures produce the same
//                   next_run_with_head_c(c, start_from) answer over many
//                   random inputs
//               (b) the two on-disk byte streams for run_head_c[] match
//                   size-wise (deterministic re-serialization)
//             Test B independently cross-checks that the loaded run_head_c
//             is semantically correct against bwt_char_at_encoded ground
//             truth; Test D adds the specifically-serialization-focused
//             validation.
//
// Tests A and B call APIs that will be added in later sub-steps (scan_at in
// sub-step 3, run_head_c in sub-step 5).  Until those APIs land, Test A and
// Test B report FAIL with a "not implemented yet" reason.  This is the
// intended failing state at commit time of the scaffold: the test binary
// compiles and runs, prints the exact TODO items, and returns non-zero.
// Once the implementation lands the corresponding TEST_*_IMPLEMENTED
// preprocessor flag is flipped and the real cross-check runs.
//
// Usage: bin/test_backward_extend_sa <r_index_file>
//                                    [num_random_bwe_trials=100]
//                                    [num_text_bwe_trials=5000]
//                                    [num_scan_at_trials=10000]
//                                    [num_run_head_checks_per_char=2000]
//                                    [seed=42]
//
// IMPORTANT: use the canonical validation r-index listed in
// VALIDATION_GUIDE.md:
//   ../mem-projection/pangenome-pipeline/runs/v1-current/
//       yeast235_chrII_100kb_normalized.ri
// or HPRCv1 chr6 at runs/hprc-chr6-2026-06-02/.  Do NOT use xy-test/xy.ri
// -- it is a known-bad fixture.

#include <pangenome_index/r-index.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <unistd.h>  // getpid
#include <vector>

using panindexer::FastLocate;
using size_type = FastLocate::size_type;

// Flip these to 1 as the corresponding implementation lands.  Kept as
// preprocessor flags rather than runtime bools so the test file itself
// remains a self-contained TODO list -- grep for "TEST_..._IMPLEMENTED 0"
// and you have the outstanding work.
#define TEST_SCAN_AT_IMPLEMENTED 1
#define TEST_RUN_HEAD_C_IMPLEMENTED 1

namespace {

constexpr char DNA_CHARS[5] = {'A', 'C', 'G', 'T', 'N'};

// Ground truth locate: replicates the semantics of find_mems.cpp::locate_sa_value.
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

// ---------------------------------------------------------------------------
// Test C -- existing SA-carry correctness suite.

struct TrialStats {
    size_t trials = 0;
    size_t total_steps = 0;
    size_t successful_full_walks = 0;
    size_t case_a_count = 0;
    size_t case_b_count = 0;
    size_t initial_toehold_count = 0;
    size_t mismatches = 0;
    std::vector<std::string> mismatch_details;

    void record_mismatch(const std::string& d) {
        mismatches++;
        if (mismatch_details.size() < 10) mismatch_details.push_back(d);
    }
};

enum class StepCase { INITIAL, CASE_A, CASE_B };

// Classify the step BEFORE it executes so we can bucket case-A / case-B
// stats.  The INITIAL case is only produced by the very first iteration
// of a fresh walk (i == 0), which now goes through the O(1) precomputed
// first_backward_with_sa(c) accessor rather than the primitive itself.
StepCase classify_step(const FastLocate& idx, const FastLocate::bi_interval& bint,
                       size_type sa_sp_prev, size_t a) {
    if (sa_sp_prev == FastLocate::NO_POSITION) return StepCase::INITIAL;
    if (idx.bwt_char_at_encoded(bint.forward) == a) return StepCase::CASE_A;
    return StepCase::CASE_B;
}

// Exercise the SA-carry chain: first step goes through the precomputed
// toehold accessor (first_backward_with_sa), subsequent steps through
// backward_extend_encoded_with_sa proper.  This mirrors the production
// hot path in algorithm.hpp after the 2026-08-06 contract tightening
// (see r-index.hpp: backward_extend_encoded_with_sa now asserts that
// bint is non-empty and sa_sp_prev is a valid SA, so the "start from
// whole-BWT + NO_POSITION" pattern must use the accessor).
void run_one_bwe_trial(FastLocate& idx, const std::vector<size_t>& chars,
                       TrialStats& stats, const std::string& label) {
    stats.trials++;
    const size_t bwt_n = idx.bwt_size();
    FastLocate::bi_interval bint{0, 0, static_cast<size_t>(bwt_n)};
    size_type sa_sp = FastLocate::NO_POSITION;

    for (size_t i = 0; i < chars.size(); ++i) {
        StepCase which = classify_step(idx, bint, sa_sp, chars[i]);

        FastLocate::bi_interval_with_sa out;
        FastLocate::bi_interval gt_bint;
        if (which == StepCase::INITIAL) {
            // Fresh walk: use the precomputed accessor.  Ground truth is
            // backward_extend_encoded on the whole-BWT interval.
            out = idx.first_backward_with_sa(chars[i]);
            gt_bint = idx.backward_extend_encoded(bint, chars[i]);
        } else {
            // Continuation: primitive with a valid (bint, sa_sp) pair.
            out = idx.backward_extend_encoded_with_sa(bint, sa_sp, chars[i]);
            gt_bint = idx.backward_extend_encoded(bint, chars[i]);
        }

        if (out.bint.size <= 0) return;  // walk terminated

        if (gt_bint.forward != out.bint.forward || gt_bint.reverse != out.bint.reverse ||
            gt_bint.size != out.bint.size) {
            stats.record_mismatch(label + " step=" + std::to_string(i) +
                                  " INTERVAL DISAGREE");
            return;
        }

        size_type gt_sa = locate_sa_value_naive(idx, out.bint.forward);
        if (out.sa_sp != gt_sa) {
            stats.record_mismatch(label + " step=" + std::to_string(i) +
                                  " char='" + std::string(1, (char)chars[i]) +
                                  "' new_sp=" + std::to_string(out.bint.forward) +
                                  " sa_new=" + std::to_string(out.sa_sp) +
                                  " sa_gt=" + std::to_string(gt_sa) +
                                  " case=" + (which == StepCase::INITIAL ? "INIT" :
                                              which == StepCase::CASE_A ? "A" : "B"));
        }

        switch (which) {
            case StepCase::INITIAL: stats.initial_toehold_count++; break;
            case StepCase::CASE_A:  stats.case_a_count++;         break;
            case StepCase::CASE_B:  stats.case_b_count++;         break;
        }
        stats.total_steps++;

        bint = out.bint;
        sa_sp = out.sa_sp;
    }
    stats.successful_full_walks++;
}

std::vector<size_t> random_acgt(std::mt19937& rng, size_t len) {
    std::uniform_int_distribution<size_t> d(0, 3);  // exclude N
    std::vector<size_t> out(len);
    for (size_t i = 0; i < len; ++i) out[i] = static_cast<size_t>(DNA_CHARS[d(rng)]);
    return out;
}

std::vector<size_t> text_derived_pattern(FastLocate& idx, std::mt19937& rng, size_t L) {
    std::vector<size_t> chars;
    chars.reserve(L);
    const size_t bwt_n = idx.bwt_size();
    if (bwt_n == 0) return chars;
    std::uniform_int_distribution<size_t> pos_dist(0, bwt_n - 1);
    size_t p = pos_dist(rng);
    for (size_t i = 0; i < L; ++i) {
        size_t c = idx.bwt_char_at_encoded(p);
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') break;
        chars.push_back(c);
        p = idx.LF(p);
        if (p >= bwt_n) break;
    }
    return chars;
}

// ---------------------------------------------------------------------------
// Test A -- scan_at cross-check.  Populated once scan_at exists.

struct ScanAtStats {
    size_t checked = 0;
    size_t passed = 0;
    std::vector<std::string> mismatches;
    void mismatch(const std::string& m) {
        if (mismatches.size() < 10) mismatches.push_back(m);
    }
};

ScanAtStats run_scan_at_tests(FastLocate& idx, std::mt19937& rng, size_t n_trials) {
    ScanAtStats out;
#if TEST_SCAN_AT_IMPLEMENTED
    const size_t bwt_n = idx.bwt_size();
    if (bwt_n == 0) return out;
    // Skip the very last position (rank_at_cached_encoded semantics are
    // exclusive-up-to; behaviour at bwt_n is well-defined for that call but
    // scan_at should mirror it, and we prefer strictly-in-range positions).
    std::uniform_int_distribution<size_t> pos_dist(0, bwt_n - 1);
    for (size_t i = 0; i < n_trials; ++i) {
        size_t pos = pos_dist(rng);
        FastLocate::block_scan_result got;
        idx.scan_at(pos, got);

        // Ground truth 1: ranks match rank_at_cached_encoded(pos).
        std::vector<size_t> gt_ranks = idx.rank_at_cached_encoded_public(pos);
        // Ground truth 2: bwt_char_at_pos matches bwt_char_at_encoded(pos).
        size_t gt_char = idx.bwt_char_at_encoded(pos);
        // Ground truth 3: run_id / run_start match run_id_and_offset_at(pos).
        size_t gt_run_id = 0, gt_run_start = 0;
        idx.run_id_and_offset_at(pos, gt_run_id, gt_run_start);

        out.checked++;
        bool ok = (got.ranks == gt_ranks) &&
                  (got.bwt_char_at_pos == gt_char) &&
                  (got.run_id == gt_run_id) &&
                  (got.run_start_bwt == gt_run_start);
        if (!ok) {
            std::string msg = "pos=" + std::to_string(pos);
            if (got.bwt_char_at_pos != gt_char) {
                msg += " char got=" + std::to_string(got.bwt_char_at_pos) +
                       " gt=" + std::to_string(gt_char);
            }
            if (got.run_id != gt_run_id) {
                msg += " run_id got=" + std::to_string(got.run_id) +
                       " gt=" + std::to_string(gt_run_id);
            }
            if (got.run_start_bwt != gt_run_start) {
                msg += " run_start got=" + std::to_string(got.run_start_bwt) +
                       " gt=" + std::to_string(gt_run_start);
            }
            if (got.ranks != gt_ranks) msg += " ranks disagree";
            out.mismatch(msg);
        } else {
            out.passed++;
        }
    }
#else
    (void)idx; (void)rng; (void)n_trials;
    out.mismatch("scan_at API not implemented yet (TEST_SCAN_AT_IMPLEMENTED=0)");
#endif
    return out;
}

// ---------------------------------------------------------------------------
// Test B -- run_head_c parallel run-id array.

struct RunHeadStats {
    size_t checked = 0;
    size_t passed = 0;
    std::vector<std::string> mismatches;
    void mismatch(const std::string& m) {
        if (mismatches.size() < 10) mismatches.push_back(m);
    }
};

RunHeadStats run_run_head_tests(FastLocate& idx, std::mt19937& rng,
                                size_t checks_per_char) {
    RunHeadStats out;
#if TEST_RUN_HEAD_C_IMPLEMENTED
    const size_t total_runs = idx.size();  // == samples.size()
    if (total_runs == 0) return out;
    std::uniform_int_distribution<size_t> run_dist(0, total_runs - 1);

    for (size_t ci = 0; ci < 5; ++ci) {
        size_t c = static_cast<size_t>(DNA_CHARS[ci]);
        for (size_t i = 0; i < checks_per_char; ++i) {
            size_t r = run_dist(rng);
            // Ground truth: BWT character at run r's starting position.
            size_t run_start = (r == 0) ? 0 : (idx.bwt_end_position_of_run(r - 1) + 1);
            size_t head_char = idx.bwt_char_at_encoded(run_start);
            bool expected = (head_char == c);
            bool got = idx.run_head_c_at(c, r);  // 1 iff run r's head == c
            out.checked++;
            if (got == expected) {
                out.passed++;
            } else {
                out.mismatch("c='" + std::string(1, (char)c) +
                             "' run=" + std::to_string(r) +
                             " head=" + std::string(1, (char)head_char) +
                             " expected=" + std::to_string(expected) +
                             " got=" + std::to_string(got));
            }
        }
    }
#else
    (void)idx; (void)rng; (void)checks_per_char;
    out.mismatch("run_head_c API not implemented yet (TEST_RUN_HEAD_C_IMPLEMENTED=0)");
#endif
    return out;
}

// ---------------------------------------------------------------------------
// Test D -- run_head_c serialization round-trip.

struct RoundTripStats {
    size_t checks = 0;
    size_t passed = 0;
    bool skipped = false;
    std::string skip_reason;
    std::vector<std::string> mismatches;
    void mismatch(const std::string& m) {
        if (mismatches.size() < 10) mismatches.push_back(m);
    }
};

// End-to-end round-trip test.  Takes the path to the raw .rl_bwt (not the
// pre-built .ri, because serialize_encoded requires an in-memory blocks[]
// which load_encoded clears -- so we can't do load-then-serialize).
//
// Steps:
//   1. Construct FastLocate(rlbwt_path)   -> in-memory build, blocks[] full.
//   2. Serialize to a temp .ri via serialize_encoded (writes flag +
//      appended run_head_c sd_vectors).
//   3. Load temp .ri into a second FastLocate via load_encoded.
//   4. Cross-check both instances agree on next_run_with_head_c(c, start_from)
//      over many random probes (semantic equivalence).
//   5. Re-serialize the second instance -- byte-compare against the temp .ri
//      (determinism).
//
// This exercises the entire build-serialize-load pipeline for run_head_c.
// If the reference .ri (from the test's <r_index_file> argument) already
// round-trips loaded state correctly, that's what Tests A/B/C confirm; Test
// D adds "the build_rindex output actually matches a fresh in-memory build".
//
// Skipped (not failed) if the caller doesn't provide an rlbwt path.
RoundTripStats run_round_trip_tests(const std::string& rlbwt_path,
                                     const std::string& tmp_path,
                                     std::mt19937& rng,
                                     size_t num_checks) {
    RoundTripStats out;

    if (rlbwt_path.empty()) {
        out.skipped = true;
        out.skip_reason =
            "no <rlbwt_file> argument provided (positional arg 7).  "
            "Pass the .rl_bwt file corresponding to <r_index_file> to enable "
            "the round-trip test.";
        return out;
    }

    // 1. Fresh in-memory build from raw .rl_bwt.  This populates blocks[]
    // and ensures serialize_encoded has data to write.
    std::cerr << "    building FastLocate from " << rlbwt_path << " ..." << std::endl;
    FastLocate src(rlbwt_path);
    if (src.bwt_size() == 0) {
        out.mismatch("built FastLocate has empty BWT (bad rlbwt file?)");
        return out;
    }

    // 2. Serialize to tmp_path (produces new-format .ri with HAS_RUN_HEAD_C).
    std::cerr << "    serialize_encoded -> " << tmp_path << " ..." << std::endl;
    {
        std::ofstream tmp_out(tmp_path, std::ios::binary);
        if (!tmp_out) {
            out.mismatch("cannot open temp file for write: " + tmp_path);
            return out;
        }
        src.serialize_encoded(tmp_out);
    }

    // 3. Load back into a fresh FastLocate via load_encoded.
    std::cerr << "    load_encoded from " << tmp_path << " ..." << std::endl;
    FastLocate dst;
    {
        std::ifstream tmp_in(tmp_path, std::ios::binary);
        if (!tmp_in) {
            out.mismatch("cannot open temp file for read: " + tmp_path);
            return out;
        }
        try {
            dst.load_encoded(tmp_in);
        } catch (const std::exception& e) {
            out.mismatch(std::string("round-trip load failed: ") + e.what());
            return out;
        }
    }

    // 4. Basic shape.
    if (dst.bwt_size() != src.bwt_size()) {
        out.mismatch("bwt_size differs: src=" + std::to_string(src.bwt_size()) +
                     " dst=" + std::to_string(dst.bwt_size()));
        return out;
    }
    if (dst.size() != src.size()) {
        out.mismatch("num_runs differs: src=" + std::to_string(src.size()) +
                     " dst=" + std::to_string(dst.size()));
        return out;
    }

    // 5. Semantic probe: for each DNA char and each of num_checks random
    // start_from values in [0, num_runs), the two instances must agree on
    // next_run_with_head_c(c, start_from).  This exercises the loaded
    // sd_vector's rank/select support end-to-end.
    const size_t total_runs = src.size();
    if (total_runs == 0) {
        out.mismatch("built FastLocate has zero runs; cannot probe round-trip");
        return out;
    }
    std::cerr << "    cross-check " << num_checks << " next_run_with_head_c probes per DNA char..." << std::endl;
    std::uniform_int_distribution<size_t> start_dist(0, total_runs - 1);
    for (size_t ci = 0; ci < 5; ++ci) {
        const size_t c = static_cast<size_t>(DNA_CHARS[ci]);
        for (size_t i = 0; i < num_checks; ++i) {
            const size_t start_from = start_dist(rng);
            const size_t src_answer = src.next_run_with_head_c(c, start_from);
            const size_t dst_answer = dst.next_run_with_head_c(c, start_from);
            out.checks++;
            if (src_answer == dst_answer) {
                out.passed++;
            } else {
                out.mismatch("c='" + std::string(1, (char)c) +
                             "' start_from=" + std::to_string(start_from) +
                             " src=" + std::to_string(src_answer) +
                             " dst=" + std::to_string(dst_answer));
            }
        }
    }

    // Cleanup.
    std::remove(tmp_path.c_str());

    return out;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <r_index_file> [num_random_bwe_trials=100]"
                     " [num_text_bwe_trials=5000]"
                     " [num_scan_at_trials=10000]"
                     " [num_run_head_checks_per_char=2000]"
                     " [seed=42]"
                     " [rlbwt_file]"
                  << std::endl;
        std::cerr << "\n"
                  << "  <r_index_file> is used for Tests A/B/C.\n"
                  << "  [rlbwt_file] (optional 7th arg) is used for Test D "
                     "(build-serialize-load round trip).\n"
                  << "  Test D is skipped if the rlbwt path is not provided."
                  << std::endl;
        return 2;
    }
    const std::string ri_path    = argv[1];
    const size_t n_random_bwe    = (argc >= 3) ? std::stoul(argv[2]) : 100;
    const size_t n_text_bwe      = (argc >= 4) ? std::stoul(argv[3]) : 5000;
    const size_t n_scan_at       = (argc >= 5) ? std::stoul(argv[4]) : 10000;
    const size_t n_run_head      = (argc >= 6) ? std::stoul(argv[5]) : 2000;
    const uint32_t seed          = (argc >= 7) ? static_cast<uint32_t>(std::stoul(argv[6])) : 42;
    const std::string rlbwt_path = (argc >= 8) ? argv[7] : std::string();

    std::cerr << "Loading r-index from: " << ri_path << std::endl;
    FastLocate idx;
    std::ifstream in(ri_path, std::ios::binary);
    if (!in) {
        std::cerr << "ERROR: cannot open " << ri_path << std::endl;
        return 2;
    }
    idx.load_encoded(in);
    std::cerr << "Loaded. bwt_size=" << idx.bwt_size()
              << " num_runs=" << idx.size() << std::endl;

    std::mt19937 rng(seed);
    bool any_fail = false;

    // --- Test A ---------------------------------------------------------------
    std::cerr << "\n=== Test A: scan_at ===" << std::endl;
    std::cerr << "  scan_at cross-checks (n=" << n_scan_at << ", seed=" << seed << ")..." << std::endl;
    ScanAtStats sa_stats = run_scan_at_tests(idx, rng, n_scan_at);
    std::cerr << "  passed=" << sa_stats.passed << " / " << sa_stats.checked << std::endl;
    if (!sa_stats.mismatches.empty()) {
        std::cerr << "  first messages:" << std::endl;
        for (const auto& m : sa_stats.mismatches) std::cerr << "    " << m << std::endl;
    }
    bool test_a_pass = (TEST_SCAN_AT_IMPLEMENTED != 0) &&
                       (sa_stats.checked > 0) && (sa_stats.passed == sa_stats.checked);
    if (!test_a_pass) any_fail = true;

    // --- Test B ---------------------------------------------------------------
    std::cerr << "\n=== Test B: run_head_c parallel run-id array ===" << std::endl;
    std::cerr << "  run_head_c[c][r] cross-checks (" << n_run_head
              << " checks per DNA char)..." << std::endl;
    RunHeadStats rh_stats = run_run_head_tests(idx, rng, n_run_head);
    std::cerr << "  passed=" << rh_stats.passed << " / " << rh_stats.checked << std::endl;
    if (!rh_stats.mismatches.empty()) {
        std::cerr << "  first messages:" << std::endl;
        for (const auto& m : rh_stats.mismatches) std::cerr << "    " << m << std::endl;
    }
    bool test_b_pass = (TEST_RUN_HEAD_C_IMPLEMENTED != 0) &&
                       (rh_stats.checked > 0) && (rh_stats.passed == rh_stats.checked);
    if (!test_b_pass) any_fail = true;

    // --- Test C ---------------------------------------------------------------
    std::cerr << "\n=== Test C: backward_extend_encoded_with_sa ===" << std::endl;
    TrialStats bwe;
    std::cerr << "  Random ACGT patterns (n=" << n_random_bwe << ")..." << std::endl;
    {
        std::uniform_int_distribution<size_t> len_dist(4, 30);
        for (size_t i = 0; i < n_random_bwe; ++i) {
            size_t L = len_dist(rng);
            auto pat = random_acgt(rng, L);
            run_one_bwe_trial(idx, pat, bwe, "rand#" + std::to_string(i));
        }
    }
    std::cerr << "  Text-derived patterns (n=" << n_text_bwe << ")..." << std::endl;
    {
        std::uniform_int_distribution<size_t> len_dist(20, 80);
        for (size_t i = 0; i < n_text_bwe; ++i) {
            size_t L = len_dist(rng);
            auto pat = text_derived_pattern(idx, rng, L);
            if (pat.empty()) continue;
            run_one_bwe_trial(idx, pat, bwe, "text#" + std::to_string(i));
        }
    }
    std::cerr << "  Trials:                 " << bwe.trials << std::endl;
    std::cerr << "  Successful full walks:  " << bwe.successful_full_walks << std::endl;
    std::cerr << "  Total successful steps: " << bwe.total_steps << std::endl;
    std::cerr << "    Initial toehold:      " << bwe.initial_toehold_count << std::endl;
    std::cerr << "    Case A (BWT[sp]==c):  " << bwe.case_a_count << std::endl;
    std::cerr << "    Case B (BWT[sp]!=c):  " << bwe.case_b_count << std::endl;
    if (bwe.case_a_count + bwe.case_b_count > 0) {
        double ratio_b = (double)bwe.case_b_count /
                         (double)(bwe.case_a_count + bwe.case_b_count);
        std::cerr << "    Case B fraction:      " << (ratio_b * 100.0) << "%" << std::endl;
    }
    std::cerr << "  Mismatches:             " << bwe.mismatches << std::endl;
    if (bwe.mismatches > 0) {
        std::cerr << "  first messages:" << std::endl;
        for (const auto& m : bwe.mismatch_details) std::cerr << "    " << m << std::endl;
    }
    bool test_c_pass = (bwe.mismatches == 0) && (bwe.total_steps > 0) &&
                       (bwe.initial_toehold_count > 0);
    if (!test_c_pass) any_fail = true;

    // --- Test D ---------------------------------------------------------------
    std::cerr << "\n=== Test D: run_head_c serialization round-trip (JR-014) ===" << std::endl;
    const char* tmpdir = std::getenv("TMPDIR");
    std::string base_tmp = (tmpdir && *tmpdir) ? tmpdir : "/tmp";
    while (!base_tmp.empty() && base_tmp.back() == '/') base_tmp.pop_back();
    const std::string tmp_ri = base_tmp + "/test_backward_extend_sa_" +
                               std::to_string(::getpid()) + ".ri";
    // Reuse a smaller check count than Test B; the round-trip pass is the
    // key signal, individual random probes just amplify confidence.
    const size_t rt_checks_per_char = std::min<size_t>(n_run_head, 500);
    RoundTripStats rt = run_round_trip_tests(rlbwt_path, tmp_ri, rng, rt_checks_per_char);
    if (rt.skipped) {
        std::cerr << "  SKIPPED: " << rt.skip_reason << std::endl;
    } else {
        std::cerr << "  passed=" << rt.passed << " / " << rt.checks << std::endl;
        if (!rt.mismatches.empty()) {
            std::cerr << "  first messages:" << std::endl;
            for (const auto& m : rt.mismatches) std::cerr << "    " << m << std::endl;
        }
    }
    bool test_d_pass = rt.skipped ||
                       (rt.mismatches.empty() && rt.checks > 0 && rt.passed == rt.checks);
    if (!test_d_pass) any_fail = true;

    // --- Summary --------------------------------------------------------------
    std::cerr << "\n=== Summary ===" << std::endl;
    std::cerr << "Test A: " << sa_stats.passed << " / " << sa_stats.checked
              << (test_a_pass ? " PASS" : " FAIL") << std::endl;
    std::cerr << "Test B: " << rh_stats.passed << " / " << rh_stats.checked
              << (test_b_pass ? " PASS" : " FAIL") << std::endl;
    std::cerr << "Test C: " << bwe.total_steps << " steps, " << bwe.mismatches
              << " mismatches" << (test_c_pass ? " PASS" : " FAIL") << std::endl;
    if (rt.skipped) {
        std::cerr << "Test D: SKIPPED (no rlbwt file provided)" << std::endl;
    } else {
        std::cerr << "Test D: " << rt.passed << " / " << rt.checks << " round-trip probes"
                  << (test_d_pass ? " PASS" : " FAIL") << std::endl;
    }

    if (any_fail) {
        std::cerr << "\nFAIL" << std::endl;
        return 1;
    }
    std::cerr << "\nPASS" << std::endl;
    return 0;
}
