#ifndef PANINDEXER_TAG_HEAD_SAMPLES_HPP
#define PANINDEXER_TAG_HEAD_SAMPLES_HPP

// TagHeadSamples: sr-index-inspired sampled SA table over tag-run heads.
//
// Loads a .tag_samples.s<N> file produced by build_tag_head_samples (v2,
// text-order deletion rule). Provides O(rank) lookup of SA values at kept
// tag-run head positions.
//
// File format (v2, produced by build_tag_head_samples.cpp):
//   [8B]  magic       = "THSAMP\0\0"
//   [8B]  version     = 2
//   [8B]  s           (stride threshold, in SA/text units)
//   [8B]  num_tag_runs
//   [8B]  kept_count
//   [8B]  bwt_size
//   [...] sdsl::sd_vector<>::serialize  -- kept_marker over tag-run rids
//   [...] sdsl::int_vector<0>::serialize -- kept-tag-head SA values
//                                            (ordered by rid ascending)
//
// Query semantics:
//   try_sa_at_tag_head(rid) returns UNKNOWN if the tag-run head at rid was
//   deleted (implicit-deleted for coincidence with a BWT-run head, or
//   explicit-deleted by the deletion rule). Otherwise returns the packed SA
//   value stored at that tag head.
//
// Callers handling the UNKNOWN case must perform an LF-walk of up to
// sample_period() steps from run_start_bwt(rid), checking at each step for
// (a) a BWT-run head (SA available from r_index.samples[]) or (b) a kept
// tag-run head. Reference implementation lives in build_tag_head_samples.cpp
// (verification spot-check loop, lines 550-571).
//
// See RESEARCH_JOURNAL.md JR-016 (storage validation) and JR-017 (LF cost).

#include <sdsl/sd_vector.hpp>
#include <sdsl/int_vector.hpp>

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <limits>
#include <string>

namespace panindexer {

class TagHeadSamples {
public:
    // 8-byte magic: ASCII "THSAMP\0\0".
    static constexpr char     MAGIC[8] = {'T','H','S','A','M','P','\0','\0'};
    static constexpr uint64_t VERSION  = 2;

    // Sentinel returned by try_sa_at_tag_head() when the queried rid was
    // deleted (either implicitly by coincidence with a BWT-run head, or
    // explicitly by the deletion rule).
    static constexpr uint64_t UNKNOWN = std::numeric_limits<uint64_t>::max();

    TagHeadSamples();
    ~TagHeadSamples() = default;

    // Disable copy (rank support holds a pointer into kept_marker; copying
    // would leave the copy's rank support dangling).
    TagHeadSamples(const TagHeadSamples&) = delete;
    TagHeadSamples& operator=(const TagHeadSamples&) = delete;
    TagHeadSamples(TagHeadSamples&&) = default;
    TagHeadSamples& operator=(TagHeadSamples&&) = default;

    // Load from a .tag_samples.s<N> file. Validates magic, version, and
    // structural consistency (kept_count matches sd_vector 1-count,
    // sa_values length matches kept_count).
    void load(std::istream& in);

    // Fast-path lookup. Returns UNKNOWN if rid is deleted (kept_marker[rid]
    // is 0); otherwise the stored packed SA value.
    inline uint64_t try_sa_at_tag_head(std::size_t rid) const {
        if (rid >= num_tag_runs_) return UNKNOWN;
        if (!kept_marker_[rid])   return UNKNOWN;
        // rank_1(rid) counts 1-bits strictly before rid; since rid itself is
        // a 1-bit, that count is the slot index in sa_values_.
        std::size_t slot = kept_rank_(rid);
        return static_cast<uint64_t>(sa_values_[slot]);
    }

    // Metadata accessors.
    inline std::size_t sample_period() const { return s_; }
    inline std::size_t num_tag_runs() const  { return num_tag_runs_; }
    inline std::size_t kept_count() const    { return kept_count_; }
    inline std::size_t bwt_size() const      { return bwt_size_; }

    // Approximate in-memory footprint, bytes.
    std::size_t memory_bytes() const;

private:
    std::size_t s_ = 0;
    std::size_t num_tag_runs_ = 0;
    std::size_t kept_count_ = 0;
    std::size_t bwt_size_ = 0;
    sdsl::sd_vector<>              kept_marker_;
    sdsl::sd_vector<>::rank_1_type kept_rank_;
    sdsl::int_vector<0>            sa_values_;

    // (Re)bind rank support after kept_marker_ is mutated.
    void rebind_supports();
};

} // namespace panindexer

#endif
