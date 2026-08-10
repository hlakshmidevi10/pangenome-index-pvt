#ifndef PANINDEXER_TAG_HEAD_SAMPLES_QUERY_HPP
#define PANINDEXER_TAG_HEAD_SAMPLES_QUERY_HPP

// Query helpers for TagHeadSamples that also touch the r-index and
// LightTagIndex. Kept in a separate header from tag_head_samples.hpp so the
// pure data class stays free of the r-index dependency and can be included
// standalone (e.g., for tests, tools that only inspect the samples file).
//
// The one function here, seed_sa_via_samples(), resolves the packed SA value
// at the head of tag-run `rid`:
//
//   1. Fast path: if the head was kept by the deletion rule, return
//      TagHeadSamples::try_sa_at_tag_head(rid) directly.
//
//   2. Slow path: LF-walk from run_start_bwt(rid) up to sample_period() steps.
//      At each step check for
//        (a) BWT-run head       -> SA available from r_index.samples[]
//        (b) kept tag-run head  -> SA from TagHeadSamples::try_sa_at_tag_head
//      Return the landed sample + step count.
//
// This is the runtime analog of the verification spot-check in
// build_tag_head_samples.cpp:550-571 (which used the naive LF()); this
// version uses LF_scan so the BWT-run-head check ("is scan.run_start_bwt
// == cur") is free as a byproduct of the fused char+rank computation.
//
// Precondition: the samples file was built against the same r-index +
// LightTagIndex that are passed here. The build guarantees that the walk
// terminates within sample_period() steps for any deleted tag-run head.
// If the invariant is violated (e.g., corrupted file, mismatched inputs)
// the function aborts.

#include "pangenome_index/r-index.hpp"
#include "pangenome_index/light_tag_index.hpp"
#include "pangenome_index/tag_head_samples.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>

namespace panindexer {

// Result of a samples-based SA resolution. Callers use `sa_value`; the
// other fields are lightweight instrumentation for aggregated profiling
// (see MEMProfilingStats.samples_* counters in find_mems.cpp).
struct SamplesResolution {
    FastLocate::size_type sa_value;
    // How many LF steps were taken (0 on fast-path hit; 1..sample_period()
    // on slow-path). Reported so callers can accumulate a step histogram.
    std::size_t lf_steps;
    // true iff the fast path (kept tag-head lookup) resolved the query.
    // false iff we fell through to the LF-walk slow path.
    bool hit_fast_path;
};

// Resolve packed SA[run_start_bwt(rid)] via samples (fast) or bounded
// LF-walk (slow, <= samples.sample_period() steps). Never returns
// TagHeadSamples::UNKNOWN.
//
// r_index passed non-const because scan_at (called by LF_scan) is a
// non-const member of FastLocate.
inline SamplesResolution seed_sa_via_samples(
    FastLocate& r_index,
    const LightTagIndex& ltag,
    const TagHeadSamples& samples,
    std::size_t rid)
{
    using size_type = FastLocate::size_type;

    // Fast path.
    uint64_t fast = samples.try_sa_at_tag_head(rid);
    if (fast != TagHeadSamples::UNKNOWN) {
        return SamplesResolution{static_cast<size_type>(fast), 0, true};
    }

    // Slow path: LF-walk from the head of this tag run.
    std::size_t cur = ltag.run_start_bwt(rid);
    const std::size_t max_steps = samples.sample_period();

    FastLocate::block_scan_result scan;
    for (std::size_t step = 0; step < max_steps; ++step) {
        std::size_t lf_pos = r_index.LF_scan(cur, scan);

        // (a) BWT-run head check -- free from scan_at (populates run_start_bwt).
        //     If cur is a BWT-run head, r_index.samples[scan.run_id] gives the
        //     packed SA at cur, and SA at the original head is that + step.
        if (scan.run_start_bwt == cur) {
            size_type sa = r_index.getSample(scan.run_id);
            return SamplesResolution{
                sa + static_cast<size_type>(step),
                step,
                false
            };
        }

        // (b) Kept tag-run head check. run_id_at is an sd_vector rank
        //     (O(log(n/m))); we accept the cost since the fast path already
        //     failed and we're in the slow path budget.
        std::size_t cur_rid = ltag.run_id_at(cur);
        if (ltag.run_start_bwt(cur_rid) == cur) {
            uint64_t sa_at_cur = samples.try_sa_at_tag_head(cur_rid);
            if (sa_at_cur != TagHeadSamples::UNKNOWN) {
                return SamplesResolution{
                    static_cast<size_type>(sa_at_cur)
                        + static_cast<size_type>(step),
                    step,
                    false
                };
            }
        }

        cur = lf_pos;
    }

    // Invariant violation: the build phase guarantees a kept sample or
    // BWT-run head lies within sample_period() LF steps of every deleted
    // tag-run head. Reaching here means the file is mismatched or corrupted.
    std::cerr << "FATAL: seed_sa_via_samples: no anchor within "
              << max_steps << " LF steps from rid=" << rid
              << " (bwt_pos=" << ltag.run_start_bwt(rid) << ")\n";
    std::abort();
}

} // namespace panindexer

#endif
