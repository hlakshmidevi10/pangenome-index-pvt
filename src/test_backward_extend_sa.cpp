// Unit test for FastLocate::backward_extend_encoded_with_sa, the SA-carrying
// backward-extend primitive added for the flipped MEM finder (see
// DESIGN_FLIPPED_MEM.md sections 4 and 5.2 Stage 1).
//
// Verification strategy (mirrors the pattern established in test_r_index_sa.cpp):
//
//   At every step of every trial walk we run backward-extend two ways:
//     (a) the new SA-carrying primitive, which returns (new_interval, sa_sp);
//     (b) the existing backward_extend_encoded followed by a fresh
//         locate_sa_value(new_sp), which is the ground truth.
//   Assert sa_sp == ground_truth at every non-empty step. Any mismatch is a
//   bug -- Step 1' of the flipped MEM finder consumes this value directly to
//   seed downstream tag emission, and a wrong SA here would corrupt biology.
//
// Trial mix:
//   1) Text-derived patterns. Pick a random BWT position, walk LF backward L
//      steps to collect a genuine substring of the text; then feed it back
//      through the SA-carrying backward-extend. These walks are guaranteed
//      to succeed at every step because the substring is a real substring
//      of the text (occurrence count >= 1). Exercises Case A most heavily
//      -- long real substrings tend to sit in wide c-heavy intervals.
//   2) Random ACGT patterns. Short (typically 4-20 char) random strings.
//      Most fail before reaching full length; still useful because failures
//      exercise the empty-result path and shorter walks exercise Case B when
//      the interval is narrow and mostly non-target-character.
//
// Also counts Case A vs Case B frequencies on successful steps. If Case B
// dominates in practice (Risk B in DESIGN_FLIPPED_MEM.md section 5.3), the
// per-step cost model of the flipped finder shifts. The ratio itself is a
// design input, not a pass/fail criterion for this stage.
//
// Usage: bin/test_backward_extend_sa <r_index_file> [num_random_trials=5000]
//                                    [num_text_trials=5000] [seed=42]
//
// IMPORTANT: use the canonical validation r-index listed in
// VALIDATION_GUIDE.md:
//   ../mem-projection/pangenome-pipeline/runs/v2-yeast235/
//       yeast235_chrII_100kb_normalized.ri
// or HPRCv1 chr6 on vesuvio (runs/hprc-chr6-2026-06-02/). Do NOT use
// xy-test/xy.ri -- it is a known-bad fixture.

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
using size_type = FastLocate::size_type;

namespace {

constexpr char DNA_CHARS[4] = {'A', 'C', 'G', 'T'};  // exclude N for pattern generation
                                                     // (matches how real reads are typically drawn)

// Ground truth locate: replicates the semantics of find_mems.cpp:locate_sa_value.
// Kept private to this test file so the test remains a self-contained oracle
// even if the find_mems.cpp helper changes shape.
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

struct TrialStats {
    size_t trials = 0;
    size_t total_steps = 0;              // total successful backward-extend steps executed
    size_t successful_full_walks = 0;    // trials that never hit empty interval
    size_t case_a_count = 0;             // BWT[old_sp] == a on the step
    size_t case_b_count = 0;             // BWT[old_sp] != a on the step
    size_t initial_toehold_count = 0;    // first-step "seed via locate_sa_value" invocations
    size_t mismatches = 0;               // sa_sp mismatches vs ground truth
    std::vector<std::string> mismatch_details;  // cap at first 10

    void record_mismatch(const std::string& d) {
        mismatches++;
        if (mismatch_details.size() < 10) mismatch_details.push_back(d);
    }
};

// Classify which case (A / B / initial) a step would take, without touching
// the new primitive. Used only for logging/statistics, so its implementation
// mirrors the primitive's own case discrimination for reporting purposes.
enum class StepCase { INITIAL, CASE_A, CASE_B };

StepCase classify_step(const FastLocate& idx, const FastLocate::bi_interval& bint,
                       size_type sa_sp_prev, size_t a) {
    if (sa_sp_prev == FastLocate::NO_POSITION) return StepCase::INITIAL;
    if (idx.bwt_char_at_encoded(bint.forward) == a) return StepCase::CASE_A;
    return StepCase::CASE_B;
}

// Run one trial: repeatedly backward-extend by `chars` (in order chars[0], chars[1], ...
// -- caller is responsible for reversing if the input pattern is left-to-right in text).
// At each step compare sa_sp against ground truth.
void run_one_trial(FastLocate& idx, const std::vector<size_t>& chars,
                   TrialStats& stats, const std::string& label) {
    stats.trials++;
    const size_t bwt_n = idx.bwt_size();
    FastLocate::bi_interval bint{0, 0, static_cast<size_t>(bwt_n)};
    size_type sa_sp = FastLocate::NO_POSITION;

    for (size_t i = 0; i < chars.size(); ++i) {
        StepCase which = classify_step(idx, bint, sa_sp, chars[i]);

        auto out = idx.backward_extend_encoded_with_sa(bint, sa_sp, chars[i]);

        if (out.bint.size <= 0) {
            // Walk terminated (interval collapsed). Not a bug; just end trial.
            return;
        }

        // Ground truth: run the plain extend + locate on the freshly-computed
        // sp. Any discrepancy is a real bug in the SA-carry logic.
        FastLocate::bi_interval gt_bint = idx.backward_extend_encoded(bint, chars[i]);
        // Structural sanity: the interval math must be identical.
        if (gt_bint.forward != out.bint.forward || gt_bint.reverse != out.bint.reverse ||
            gt_bint.size != out.bint.size) {
            stats.record_mismatch(label + " step=" + std::to_string(i) +
                                  " INTERVAL DISAGREE: new=(k=" + std::to_string(out.bint.forward) +
                                  ",l=" + std::to_string(out.bint.reverse) +
                                  ",s=" + std::to_string(out.bint.size) + ") gt=(k=" +
                                  std::to_string(gt_bint.forward) + ",l=" +
                                  std::to_string(gt_bint.reverse) + ",s=" +
                                  std::to_string(gt_bint.size) + ")");
            return;  // no point continuing this trial
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

        // Bookkeeping: which case did this step actually exercise?
        switch (which) {
            case StepCase::INITIAL: stats.initial_toehold_count++; break;
            case StepCase::CASE_A:  stats.case_a_count++; break;
            case StepCase::CASE_B:  stats.case_b_count++; break;
        }
        stats.total_steps++;

        // Advance state for next step.
        bint = out.bint;
        sa_sp = out.sa_sp;
    }
    stats.successful_full_walks++;
}

// Generate a random ACGT pattern of length len (uniform).
std::vector<size_t> random_acgt(std::mt19937& rng, size_t len) {
    std::uniform_int_distribution<size_t> d(0, 3);
    std::vector<size_t> out(len);
    for (size_t i = 0; i < len; ++i) out[i] = static_cast<size_t>(DNA_CHARS[d(rng)]);
    return out;
}

// Walk L steps of LF from a random BWT position, collecting the BWT characters
// visited. The collected sequence, when fed as the input to a chain of
// backward extends, corresponds to a real substring of the underlying text
// and is therefore guaranteed to be findable at every step (occurrence count
// >= 1 throughout the walk).
//
// This is the "text-derived" trial mode -- it exercises the deep-walk case
// where Case A dominates, which is the main win path for the flipped finder.
std::vector<size_t> text_derived_pattern(FastLocate& idx, std::mt19937& rng, size_t L) {
    std::vector<size_t> chars;
    chars.reserve(L);
    const size_t bwt_n = idx.bwt_size();
    if (bwt_n == 0) return chars;

    // Skip endmarker positions when picking start: endmarker rows have SA
    // values at seq-boundaries and their BWT character is NENDMARKER, which
    // won't be findable via a DNA extend anyway.
    std::uniform_int_distribution<size_t> pos_dist(0, bwt_n - 1);
    size_t p = pos_dist(rng);
    for (size_t i = 0; i < L; ++i) {
        size_t c = idx.bwt_char_at_encoded(p);
        // Bail out early if we hit a non-DNA character (e.g., NENDMARKER at a
        // sequence boundary). The already-collected prefix is still a valid
        // real substring.
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') break;
        chars.push_back(c);
        p = idx.LF(p);
        if (p >= bwt_n) break;
    }
    return chars;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <r_index_file> [num_random_trials=5000] [num_text_trials=5000] [seed=42]"
                  << std::endl;
        return 2;
    }
    const std::string ri_path = argv[1];
    const size_t n_random = (argc >= 3) ? std::stoul(argv[2]) : 5000;
    const size_t n_text = (argc >= 4) ? std::stoul(argv[3]) : 5000;
    const uint32_t seed = (argc >= 5) ? static_cast<uint32_t>(std::stoul(argv[4])) : 42;

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
    TrialStats stats;

    // --- Suite 1: random ACGT patterns of length 4-30 ---------------------
    std::cerr << "\n=== Random ACGT patterns (n=" << n_random << ") ===" << std::endl;
    {
        std::uniform_int_distribution<size_t> len_dist(4, 30);
        for (size_t i = 0; i < n_random; ++i) {
            size_t L = len_dist(rng);
            auto pat = random_acgt(rng, L);
            run_one_trial(idx, pat, stats, "rand#" + std::to_string(i));
        }
    }

    // --- Suite 2: text-derived patterns (guaranteed long walks) -----------
    std::cerr << "\n=== Text-derived patterns (n=" << n_text << ") ===" << std::endl;
    {
        // Length 20-80 covers typical read-length backward walks used by the
        // MEM finder inner loop; adjust upward to stress-test if desired.
        std::uniform_int_distribution<size_t> len_dist(20, 80);
        for (size_t i = 0; i < n_text; ++i) {
            size_t L = len_dist(rng);
            auto pat = text_derived_pattern(idx, rng, L);
            if (pat.empty()) continue;  // hit a boundary immediately
            run_one_trial(idx, pat, stats, "text#" + std::to_string(i));
        }
    }

    // --- Report -----------------------------------------------------------
    std::cerr << "\n=== Summary ===" << std::endl;
    std::cerr << "Trials:                  " << stats.trials << std::endl;
    std::cerr << "Successful full walks:   " << stats.successful_full_walks << std::endl;
    std::cerr << "Total successful steps:  " << stats.total_steps << std::endl;
    std::cerr << "  Initial toehold (seed): " << stats.initial_toehold_count << std::endl;
    std::cerr << "  Case A (BWT[sp]==c):    " << stats.case_a_count << std::endl;
    std::cerr << "  Case B (BWT[sp]!=c):    " << stats.case_b_count << std::endl;
    if (stats.case_a_count + stats.case_b_count > 0) {
        double ratio_b = (double)stats.case_b_count /
                         (double)(stats.case_a_count + stats.case_b_count);
        std::cerr << "  Case B fraction:        " << (ratio_b * 100.0) << "%" << std::endl;
    }
    std::cerr << "Mismatches:              " << stats.mismatches << std::endl;

    if (stats.mismatches > 0) {
        std::cerr << "\nFirst mismatches:" << std::endl;
        for (const auto& d : stats.mismatch_details) {
            std::cerr << "  " << d << std::endl;
        }
        std::cerr << "\nFAIL" << std::endl;
        return 1;
    }
    if (stats.total_steps == 0) {
        std::cerr << "\nERROR: no successful backward-extend steps executed -- "
                     "test is not actually verifying anything." << std::endl;
        return 1;
    }
    std::cerr << "\nPASS" << std::endl;
    return 0;
}
