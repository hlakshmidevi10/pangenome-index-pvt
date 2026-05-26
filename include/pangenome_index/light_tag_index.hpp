#ifndef PANINDEXER_LIGHT_TAG_INDEX_HPP
#define PANINDEXER_LIGHT_TAG_INDEX_HPP

// LightTagIndex: a minimal tag index containing ONLY the BWT-position run-start
// bitvector (bwt_intervals from the full TagArray). No tag values (node_id,
// offset, strand) are stored, so it cannot answer "what graph position does
// BWT position p map to". It only answers:
//   - run_id_at(p):   which tag run contains BWT position p
//   - run_start_bwt(r): BWT start position of run r
//   - num_runs()
//   - bwt_size()
//
// This is sufficient for find_mems --lightweight-tags, which only needs to know
// where tag-run boundaries fall in order to emit one (seq_id, path_bp) entry
// per intersecting tag run (no graph-position dedup).
//
// File format (.ltags):
//   [8 bytes]  magic   = "LTAGS\0\0\0"      (LightTagIndex::MAGIC)
//   [8 bytes]  version = 1                  (LightTagIndex::VERSION)
//   [8 bytes]  bwt_size                     (sd_vector size in bits)
//   [...]      sdsl::sd_vector<>::serialize output (bwt_intervals)

#include <sdsl/sd_vector.hpp>
#include <iosfwd>
#include <string>
#include <cstdint>
#include <cstddef>

namespace panindexer {

class LightTagIndex {
public:
    // 8-byte magic: ASCII "LTAGS\0\0\0". On little-endian disks this is the
    // byte sequence 'L','T','A','G','S',0,0,0 which equals 0x0000005347415414
    // when interpreted as a host uint64_t on little-endian platforms. We write
    // and read it as raw bytes (memcpy) so endianness does not matter.
    static constexpr char     MAGIC[8] = {'L','T','A','G','S','\0','\0','\0'};
    static constexpr uint64_t VERSION  = 1;

    LightTagIndex();
    ~LightTagIndex() = default;

    // Disable copy (sd_vector supports it, but the rank/select supports hold
    // pointers into bwt_intervals — copying would silently leave them dangling).
    LightTagIndex(const LightTagIndex&) = delete;
    LightTagIndex& operator=(const LightTagIndex&) = delete;
    LightTagIndex(LightTagIndex&&) = default;
    LightTagIndex& operator=(LightTagIndex&&) = default;

    // Serialize to / load from .ltags file (binary, includes magic + version header).
    void save(std::ostream& out) const;
    void load(std::istream& in);

    // Build by copying bwt_intervals out of an existing compact TagArray file.
    // Convenience wrapper: opens the file, loads a full TagArray, copies bwt_intervals.
    void build_from_compact_tags_file(const std::string& compact_tags_path);

    // Direct construction from an sd_vector (used by tests / alternate builders).
    void set_bwt_intervals(sdsl::sd_vector<>&& v);

    // BWT pos -> run id (0-based). Caller must ensure 0 <= bwt_pos < bwt_size().
    inline std::size_t run_id_at(std::size_t bwt_pos) const {
        return bwt_intervals_rank(bwt_pos + 1) - 1;
    }
    // run id -> BWT start position of that run.
    // Caller must ensure 0 <= run_id < num_runs().
    inline std::size_t run_start_bwt(std::size_t run_id) const {
        return bwt_intervals_select(run_id + 1);
    }

    inline std::size_t num_runs() const  { return num_runs_cached; }
    inline std::size_t bwt_size() const  { return bwt_intervals.size(); }

    // Approximate in-memory footprint, bytes. Sum of sdsl::size_in_bytes over
    // bwt_intervals + its rank/select supports.
    std::size_t memory_bytes() const;

    // Read-only access for validators / tests.
    const sdsl::sd_vector<>& get_bwt_intervals() const { return bwt_intervals; }

private:
    sdsl::sd_vector<>                bwt_intervals;
    sdsl::sd_vector<>::rank_1_type   bwt_intervals_rank;
    sdsl::sd_vector<>::select_1_type bwt_intervals_select;
    std::size_t                      num_runs_cached = 0;

    // (Re)bind the rank/select supports after bwt_intervals is mutated.
    void rebind_supports();
};

} // namespace panindexer
#endif
