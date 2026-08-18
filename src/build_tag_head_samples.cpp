// build_tag_head_samples:
//   Build sampled SA tables over tag-run heads using an sr-index-inspired
//   deletion rule in TEXT (SA) ORDER, so that from any deleted tag-run head
//   t, walking LF backward at most s steps lands on a stored sample (either
//   a kept tag-run head OR a BWT-run head, the latter being always sampled
//   in the r-index).
//
// Note on SA representation: throughout this binary, "SA value" means the
// PACKED value (seq_id * max_length + seq_offset) that the r-index uses
// internally. locateNext() and getSample() both operate on packed values.
// Adjacent BWT positions within a single sequence differ by 1 in packed
// space, so the "LF decrements SA by 1" invariant is preserved for the
// deletion-rule bound (< s in packed distance == < s in LF steps within a
// sequence). BWT runs never cross sequence boundaries (endmarker chars
// interrupt runs), so intra-run walks always stay in one sequence.
//
// Deletion rule (single left-to-right pass over candidates sorted by SA
// value, in ascending text order):
//   Enumerate candidates = tag-run heads (deletable) + BWT-run heads
//   (unremovable anchors). Sort by SA value. For each tag-run head s_i,
//   delete iff (SA[s_{i+1}] - SA[s_prev_kept] < s), where s_prev_kept is
//   the most recent kept sample in the scan. BWT-run heads are always kept.
//
// Guarantee: LF(pos) decrements SA by exactly 1. Walking LF backward from
// any deleted tag-run head therefore visits BWT positions whose SA values
// are contiguous and decreasing. The deletion rule ensures a kept sample
// exists at some SA in [SA[t] - s + 1, SA[t]], so the walk lands on that
// sample's exact BWT position within < s LF steps.
//
// Output file format (per s value): <base>.tag_samples.s<s>
//   [8B]  magic   = "THSAMP\0\0"
//   [8B]  version = 2  (v2 = text-order deletion; v1 was BWT-order)
//   [8B]  s (stride threshold, in SA / text units)
//   [8B]  num_tag_runs   (from LightTagIndex, for consistency check)
//   [8B]  kept_count     (number of 1-bits in kept_marker; also length of sa_values)
//   [8B]  bwt_size       (n)
//   [...] sd_vector<>::serialize -- kept_marker over tag-run rids
//   [...] int_vector<0>::serialize -- kept-tag-head SA values, ordered by rid ascending
//         (so query is: slot = kept_marker.rank_1(rid); return sa_values[slot])
//
// Coincident case: if a tag-run head sits at the same BWT position as a
// BWT-run head, its SA is already available from r_index.samples[]. We treat
// it as anchor-only (tag-head implicitly deleted from kept_marker). Query
// slow-path handles this via the "is pos a BWT-run head?" check at step 0.

#include "pangenome_index/r-index.hpp"
#include "pangenome_index/light_tag_index.hpp"

#include <sdsl/sd_vector.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace panindexer;

namespace {

constexpr char     THSAMP_MAGIC[8] = {'T','H','S','A','M','P','\0','\0'};
constexpr uint64_t THSAMP_VERSION  = 2;  // v2 = text-order deletion

void usage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " <ri_file> <ltags_file> <output_base> "
        "[--s LIST] [--skip-verify]\n"
        "\n"
        "  <ri_file>       Input r-index (.ri)\n"
        "  <ltags_file>    Input light tag index (.ltags)\n"
        "  <output_base>   Output filename prefix; writes <base>.tag_samples.s<N> per s\n"
        "\n"
        "Options:\n"
        "  --s LIST        Comma-separated s values (default: 32,64,128,256)\n"
        "  --skip-verify   Skip post-build spot-check (default: verify 1000 kept + 100 deleted)\n"
        "  -h, --help      Show this help\n";
}

long long file_size_bytes(const std::string& path) {
    struct stat st{};
    if (::stat(path.c_str(), &st) != 0) return -1;
    return static_cast<long long>(st.st_size);
}

std::vector<size_t> parse_s_list(const std::string& s) {
    std::vector<size_t> out;
    std::string tok;
    std::stringstream ss(s);
    while (std::getline(ss, tok, ',')) {
        if (tok.empty()) continue;
        long long v = std::stoll(tok);
        if (v <= 0) throw std::runtime_error("s values must be positive, got: " + tok);
        out.push_back(static_cast<size_t>(v));
    }
    if (out.empty()) throw std::runtime_error("--s list is empty");
    std::sort(out.begin(), out.end());
    return out;
}

std::string human_bytes(std::size_t bytes) {
    const char* units[] = {"B", "KB", "MB", "GB", "TB"};
    double v = static_cast<double>(bytes);
    int u = 0;
    while (v >= 1024.0 && u < 4) { v /= 1024.0; ++u; }
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << v << " " << units[u];
    return oss.str();
}

// Compute SA at BWT position `pos` via the today-standard path.
// Reference implementation, used ONLY for spot-check verification.
uint64_t locate_sa_value_ref(FastLocate& r_index, size_t pos) {
    size_t run_id = 0;
    size_t offset_of_first = 0;
    r_index.run_id_and_offset_at(pos, run_id, offset_of_first);
    FastLocate::size_type sa = r_index.getSample(run_id);
    while (offset_of_first < pos) {
        sa = r_index.locateNext(sa);
        offset_of_first++;
    }
    return static_cast<uint64_t>(sa);
}

// Return true iff `pos` is the head (leftmost BWT position) of some BWT run.
// Cheap: one run_id_and_offset_at lookup; head iff offset_of_first == pos.
bool is_bwt_run_head(FastLocate& r_index, size_t pos) {
    size_t rid = 0, off = 0;
    r_index.run_id_and_offset_at(pos, rid, off);
    return off == pos;
}

// Return true iff `pos` is a tag-run head. bwt_rid_out set when true.
bool is_bwt_run_head_with_id(FastLocate& r_index, size_t pos, size_t& rid_out) {
    size_t off = 0;
    r_index.run_id_and_offset_at(pos, rid_out, off);
    return off == pos;
}

bool is_tag_run_head(const LightTagIndex& ltag, size_t pos, size_t& rid_out) {
    // Guard against pos >= bwt_size (sentinel region).
    if (pos >= ltag.bwt_size()) return false;
    rid_out = ltag.run_id_at(pos);
    return ltag.run_start_bwt(rid_out) == pos;
}

// A merged candidate for the deletion scan. Sorted by sa_value ascending.
// `kind`: 0 = BWT anchor (unremovable), 1 = tag-run head (deletable).
// For kind=0, `rid` = bwt_run_id (unused post-scan, kept for debug).
// For kind=1, `rid` = tag_run_id, referenced by kept_marker output.
struct Candidate {
    uint64_t sa_value;
    uint64_t rid;
    uint8_t  kind;  // 0=anchor, 1=tag
    uint8_t  _pad[7];  // pad to 24 bytes for alignment cleanliness

    bool operator<(const Candidate& o) const { return sa_value < o.sa_value; }
};
static_assert(sizeof(Candidate) == 24, "Candidate size unexpected");

// Iterator over BWT-run heads. Emits (bwt_pos, bwt_run_id) in BWT-position
// ascending order by decoding blocks sequentially.
struct BwtRunHeadIter {
    FastLocate* r_index;
    size_t num_blocks;
    size_t block_idx = 0;
    size_t local_run_idx_in_block = 0;
    size_t global_run_id = 0;
    size_t cur_block_start = 0;
    size_t cur_offset_in_block = 0;
    std::vector<std::pair<size_t, size_t>> block_runs;
    bool done_flag = false;

    explicit BwtRunHeadIter(FastLocate& ri) : r_index(&ri) {
        num_blocks = ri.num_blocks();
        load_current_block();
    }

    void load_current_block() {
        while (block_idx < num_blocks) {
            block_runs.clear();
            r_index->get_block_runs(block_idx, block_runs);
            if (!block_runs.empty()) {
                cur_block_start = r_index->blocks_start_select_1(block_idx + 1);
                local_run_idx_in_block = 0;
                cur_offset_in_block = 0;
                return;
            }
            block_idx++;
        }
        done_flag = true;
    }

    bool done() const { return done_flag; }
    size_t pos() const { return cur_block_start + cur_offset_in_block; }
    size_t run_id() const { return global_run_id; }
    size_t run_length() const { return block_runs[local_run_idx_in_block].second; }

    void advance() {
        if (done_flag) return;
        cur_offset_in_block += block_runs[local_run_idx_in_block].second;
        local_run_idx_in_block++;
        global_run_id++;
        if (local_run_idx_in_block >= block_runs.size()) {
            block_idx++;
            load_current_block();
        }
    }
};

}  // namespace

int main(int argc, char** argv) {
    if (argc < 4) { usage(argv[0]); return 2; }

    std::string ri_path = argv[1];
    std::string ltags_path = argv[2];
    std::string out_base = argv[3];
    std::vector<size_t> s_values = {32, 64, 128, 256};
    bool skip_verify = false;

    for (int i = 4; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
        else if (a == "--s" && i + 1 < argc) { s_values = parse_s_list(argv[++i]); }
        else if (a == "--skip-verify") { skip_verify = true; }
        else { std::cerr << "Unknown option: " << a << "\n"; usage(argv[0]); return 2; }
    }

    std::cerr << "========================================\n"
              << "build_tag_head_samples (v2, text-order)\n"
              << "========================================\n"
              << "  ri:          " << ri_path << "\n"
              << "  ltags:       " << ltags_path << "\n"
              << "  output base: " << out_base << "\n"
              << "  s values:   ";
    for (auto s : s_values) std::cerr << " " << s;
    std::cerr << "\n"
              << "  verify:      " << (skip_verify ? "skipped" : "1000 kept + 100 deleted per s") << "\n"
              << "----------------------------------------\n";

    // --- Load r-index ---
    FastLocate r_index;
    {
        auto t0 = std::chrono::high_resolution_clock::now();
        std::cerr << "Loading r-index..." << std::endl;
        std::ifstream in(ri_path, std::ios::binary);
        if (!in) { std::cerr << "ERROR: cannot open " << ri_path << "\n"; return 1; }
        r_index.load_encoded(in);
        auto t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        std::cerr << "  BWT size:      " << r_index.bwt_size() << "\n"
                  << "  BWT runs:      " << r_index.tot_runs() << "\n"
                  << "  num_blocks:    " << r_index.num_blocks() << "\n"
                  << "  samples.size():" << r_index.samples.size() << "\n"
                  << "  samples.width:" << r_index.samples.width() << " bits\n"
                  << "  load time:     " << dt.count() << " s\n";
    }

    LightTagIndex ltag;
    {
        auto t0 = std::chrono::high_resolution_clock::now();
        std::cerr << "Loading light tag index..." << std::endl;
        std::ifstream in(ltags_path, std::ios::binary);
        if (!in) { std::cerr << "ERROR: cannot open " << ltags_path << "\n"; return 1; }
        ltag.load(in);
        auto t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> dt = t1 - t0;
        std::cerr << "  bwt_size:      " << ltag.bwt_size() << "\n"
                  << "  num_tag_runs:  " << ltag.num_runs() << "\n"
                  << "  load time:     " << dt.count() << " s\n";
    }

    if (ltag.bwt_size() != r_index.bwt_size() && ltag.bwt_size() != r_index.bwt_size() + 1) {
        std::cerr << "ERROR: BWT size mismatch: ri=" << r_index.bwt_size()
                  << " vs ltag=" << ltag.bwt_size() << " (expected equal or +1)\n";
        return 1;
    }
    const size_t n = r_index.bwt_size();
    const size_t num_tag_runs = ltag.num_runs();
    const size_t num_bwt_runs = r_index.tot_runs();

    std::cerr << "----------------------------------------\n";

    // =========================================================
    // Phase 1: enumerate all candidates with their SA values.
    // =========================================================
    //
    // Walk the BWT sequentially in BWT-position order. Maintain a running
    // SA cursor initialized at each BWT-run head via samples[bwt_rid]. At
    // each BWT run start, emit a candidate (SA=sample, kind=anchor). While
    // walking locateNext through the run, at every BWT position that is
    // ALSO a tag-run head, emit a candidate (SA=cursor, kind=tag).
    //
    // Coincident case: if a BWT-run head coincides with a tag-run head
    // (i.e., pos is start of both), we emit ONLY the anchor candidate --
    // the tag head is implicitly deleted (its SA is available for free).

    auto t_phase1_start = std::chrono::high_resolution_clock::now();
    std::cerr << "Phase 1: enumerate candidates with SA values..." << std::endl;

    std::vector<Candidate> candidates;
    // Upper bound: num_tag_runs + num_bwt_runs. Reserve a bit less to save
    // RAM under coincident collapse; we'll grow as needed.
    candidates.reserve(num_tag_runs + num_bwt_runs);

    // Iterator over tag-run heads in BWT-position order.
    size_t next_tag_rid = 0;
    auto next_tag_pos = [&]() -> size_t {
        while (next_tag_rid < num_tag_runs) {
            size_t p = ltag.run_start_bwt(next_tag_rid);
            if (p < n) return p;
            next_tag_rid++;  // skip out-of-range sentinel(s)
        }
        return SIZE_MAX;
    };

    BwtRunHeadIter bwt_iter(r_index);
    size_t total_locate_next_phase1 = 0;
    const size_t samples_size = r_index.samples.size();
    const size_t bwt_n        = r_index.bwt_size();
    const size_t ltag_num_runs = ltag.num_runs();
    size_t bwt_runs_walked = 0;
    // Progress cadence: log ~40 lines end-to-end so a multi-hour walk still
    // gives a heartbeat. Round to a power-of-two-ish for cheap modulo.
    const size_t progress_every = std::max<size_t>(1'000'000, samples_size / 40);

    while (!bwt_iter.done()) {
        size_t run_start = bwt_iter.pos();
        size_t run_len   = bwt_iter.run_length();
        size_t run_end   = run_start + run_len;   // exclusive
        size_t bwt_rid   = bwt_iter.run_id();

        // Defensive: BwtRunHeadIter numbers global_run_id densely 0..N-1 by
        // walking blocks sequentially. samples[] is sized to r_index.tot_runs().
        // If chr1 has any partially-filled encoded blocks (or any drift
        // between block-walk and samples-array sizing) this fires with a
        // diagnostic instead of a bare sdsl assert.
        if (bwt_rid >= samples_size) {
            std::cerr << "\nERROR: BwtRunHeadIter overshot samples[]:\n"
                      << "  bwt_rid (from iterator) = " << bwt_rid << "\n"
                      << "  samples.size()          = " << samples_size << "\n"
                      << "  bwt_runs_walked so far  = " << bwt_runs_walked << "\n"
                      << "  block_idx               = " << bwt_iter.block_idx << "\n"
                      << "  num_blocks              = " << r_index.num_blocks() << "\n"
                      << "  local_run_idx_in_block  = " << bwt_iter.local_run_idx_in_block << "\n"
                      << "  block_runs.size()       = " << bwt_iter.block_runs.size() << "\n"
                      << "  run_start (BWT pos)     = " << run_start << "\n"
                      << "  run_len                 = " << run_len << "\n";
            return 5;
        }
        if (run_end > bwt_n + 1) {
            // ltag.bwt_size() can be r_index.bwt_size() OR r_index.bwt_size()+1
            // (see line-283 check); anything beyond is a sign of block-walk drift.
            std::cerr << "\nERROR: run_end (" << run_end << ") exceeds bwt_size+1 ("
                      << (bwt_n + 1) << ") at bwt_rid=" << bwt_rid
                      << " run_start=" << run_start << " run_len=" << run_len << "\n";
            return 5;
        }

        // Emit anchor at run_start.
        FastLocate::size_type sa_cursor = r_index.getSample(bwt_rid);
        candidates.push_back({static_cast<uint64_t>(sa_cursor),
                              static_cast<uint64_t>(bwt_rid), 0, {}});

        // Check tag-head at run_start: if coincident with anchor, skip
        // (implicit-delete) by advancing tag iterator past it.
        {
            size_t tp = next_tag_pos();
            if (tp == run_start) {
                if (next_tag_rid >= ltag_num_runs) {
                    std::cerr << "\nERROR: tag rid overflow at coincident-collapse: "
                              << "next_tag_rid=" << next_tag_rid
                              << " ltag.num_runs()=" << ltag_num_runs << "\n";
                    return 5;
                }
                next_tag_rid++;
            }
        }

        // Walk locateNext through the interior of this BWT run, emitting
        // tag-head candidates at their exact BWT positions.
        size_t cur_pos = run_start;
        while (true) {
            size_t tp = next_tag_pos();
            if (tp == SIZE_MAX || tp >= run_end) break;
            // Defensive: tp must lie inside the current BWT run (strict lower
            // bound already enforced by run_start emit above).
            if (tp < run_start) {
                std::cerr << "\nERROR: tag pos " << tp << " precedes current run_start "
                          << run_start << " (rid=" << next_tag_rid << ")\n";
                return 5;
            }
            // Advance the cursor to tp.
            while (cur_pos < tp) {
                sa_cursor = r_index.locateNext(sa_cursor);
                total_locate_next_phase1++;
                cur_pos++;
            }
            // Emit tag candidate at cur_pos = tp.
            if (next_tag_rid >= ltag_num_runs) {
                std::cerr << "\nERROR: tag rid overflow in interior emit: "
                          << "next_tag_rid=" << next_tag_rid
                          << " ltag.num_runs()=" << ltag_num_runs
                          << " tp=" << tp << " run=[" << run_start << "," << run_end << ")\n";
                return 5;
            }
            candidates.push_back({static_cast<uint64_t>(sa_cursor),
                                  static_cast<uint64_t>(next_tag_rid), 1, {}});
            next_tag_rid++;
        }

        bwt_runs_walked++;
        if (bwt_runs_walked % progress_every == 0) {
            auto t_now = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> dt = t_now - t_phase1_start;
            double pct = 100.0 * static_cast<double>(bwt_runs_walked) / static_cast<double>(samples_size);
            std::cerr << "  [phase1] walked " << bwt_runs_walked << " / " << samples_size
                      << " BWT runs (" << std::fixed << std::setprecision(2) << pct
                      << "%), pos=" << run_start
                      << ", candidates=" << candidates.size()
                      << ", locateNext=" << total_locate_next_phase1
                      << ", t=" << std::fixed << std::setprecision(1) << dt.count() << "s\n";
        }

        bwt_iter.advance();
    }

    auto t_phase1_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> phase1_dt = t_phase1_end - t_phase1_start;
    std::cerr << "  candidates emitted: " << candidates.size()
              << " (anchors + tag-heads, before coincident collapse counting)\n"
              << "  locateNext calls:   " << total_locate_next_phase1 << "\n"
              << "  phase 1 time:       " << phase1_dt.count() << " s ("
              << (total_locate_next_phase1 ? phase1_dt.count() * 1e6 / total_locate_next_phase1 : 0.0)
              << " us/locateNext)\n";

    // Sort by SA ascending. Ties should be impossible (SA is a bijection).
    auto t_sort_start = std::chrono::high_resolution_clock::now();
    std::cerr << "  sorting " << candidates.size() << " candidates by SA..." << std::endl;
    std::sort(candidates.begin(), candidates.end());
    auto t_sort_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> sort_dt = t_sort_end - t_sort_start;
    std::cerr << "  sort time:          " << sort_dt.count() << " s\n";

    // Sanity: no ties.
    for (size_t i = 1; i < candidates.size(); ++i) {
        if (candidates[i].sa_value == candidates[i-1].sa_value) {
            std::cerr << "ERROR: SA collision at index " << i
                      << " (SA=" << candidates[i].sa_value << ")\n";
            return 3;
        }
    }
    std::cerr << "  SA uniqueness:      OK\n"
              << "----------------------------------------\n";

    // =========================================================
    // Phase 2 + 3 + 4: for each s value, apply deletion, build, verify, write.
    // =========================================================
    for (size_t s : s_values) {
        std::cerr << "\n=========================================\n"
                  << "s = " << s << " (text-order distance)\n"
                  << "-----------------------------------------\n";
        auto t_pass_start = std::chrono::high_resolution_clock::now();

        // Phase 2: deletion scan in text (SA) order.
        // For each tag-head candidate s_i, delete iff
        //   SA[s_{i+1}] - SA[prev_kept] < s.
        // BWT anchors are always kept.
        //
        // Output: two parallel vectors indexed by "kept slot":
        //   kept_rid[k]  = tag_run_id of the k-th surviving tag-run head
        //                  (in text order, NOT rid order)
        //   kept_sa[k]   = SA value of that tag-run head
        // We'll re-sort by rid before serialize (int_vector indexed by rank).
        std::vector<std::pair<uint64_t, uint64_t>> kept_tag;  // (rid, sa)
        kept_tag.reserve(candidates.size() / 8);

        size_t prev_kept_sa = 0;
        bool have_prev = false;
        size_t deleted_count = 0;

        for (size_t i = 0; i < candidates.size(); ++i) {
            const Candidate& c = candidates[i];
            if (c.kind == 0) {
                // BWT anchor: always keep. NOT written to sa_values (r-index
                // already has samples[]). Just updates prev_kept_sa.
                prev_kept_sa = c.sa_value;
                have_prev = true;
                continue;
            }
            // Tag-run head candidate. Peek next candidate.
            uint64_t next_sa = (i + 1 < candidates.size()) ? candidates[i + 1].sa_value : UINT64_MAX;
            bool keep = true;
            if (have_prev && next_sa != UINT64_MAX) {
                if (next_sa - prev_kept_sa < s) keep = false;
            }
            if (keep) {
                kept_tag.push_back({c.rid, c.sa_value});
                prev_kept_sa = c.sa_value;
                have_prev = true;
            } else {
                deleted_count++;
            }
        }

        const size_t total_tag_candidates = kept_tag.size() + deleted_count;
        std::cerr << "  tag-head candidates: " << total_tag_candidates
                  << " (excludes " << (num_tag_runs - total_tag_candidates)
                  << " coincident with BWT anchors)\n"
                  << "  kept: " << kept_tag.size()
                  << " (" << std::fixed << std::setprecision(2)
                  << 100.0 * kept_tag.size() / (total_tag_candidates ? total_tag_candidates : 1)
                  << "% of candidates, "
                  << 100.0 * kept_tag.size() / (num_tag_runs ? num_tag_runs : 1)
                  << "% of all tag runs)\n"
                  << "  deleted: " << deleted_count << "\n";

        // Sort kept_tag by rid ascending for output ordering.
        auto t_sort2_start = std::chrono::high_resolution_clock::now();
        std::sort(kept_tag.begin(), kept_tag.end(),
                  [](const std::pair<uint64_t, uint64_t>& a,
                     const std::pair<uint64_t, uint64_t>& b) {
                      return a.first < b.first;
                  });
        auto t_sort2_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> sort2_dt = t_sort2_end - t_sort2_start;
        std::cerr << "  sort by rid:        " << sort2_dt.count() << " s\n";

        // Phase 3: build sd_vector kept_marker + int_vector sa_values.
        const size_t kept_count = kept_tag.size();
        sdsl::sd_vector<> kept_marker;
        {
            sdsl::sd_vector_builder builder(num_tag_runs, kept_count);
            for (const auto& p : kept_tag) builder.set(p.first);
            kept_marker = sdsl::sd_vector<>(builder);
        }

        // SA values are PACKED (seq_id * max_length + seq_offset), not raw
        // text positions. Their range is [0, n_seq * max_length), which can
        // be considerably larger than n. Use the same bit width as the
        // r-index's own samples[] table so we never truncate.
        const size_t bits_per_sa = r_index.samples.width();
        sdsl::int_vector<0> sa_iv(kept_count, 0, bits_per_sa);
        for (size_t k = 0; k < kept_count; ++k) sa_iv[k] = kept_tag[k].second;

        const size_t sd_bytes = sdsl::size_in_bytes(kept_marker);
        const size_t iv_bytes = sdsl::size_in_bytes(sa_iv);

        std::cerr << "  bits/sa:            " << bits_per_sa << "\n"
                  << "  kept_marker (sdv):  " << human_bytes(sd_bytes)
                  << " (" << (kept_count ? sd_bytes * 8.0 / kept_count : 0.0) << " bits/kept)\n"
                  << "  sa_values (int_v):  " << human_bytes(iv_bytes) << "\n";

        // Verification passes.
        if (!skip_verify && kept_count > 0) {
            std::cerr << "  verifying..." << std::endl;
            std::mt19937_64 rng(0xC0FFEEuLL + s);

            // Support for kept_marker: need rank + select
            sdsl::sd_vector<>::rank_1_type km_rank(&kept_marker);
            sdsl::sd_vector<>::select_1_type km_select(&kept_marker);

            // (a) Spot-check kept-tag-head SAs vs locate_sa_value_ref.
            size_t n_kept_checks = std::min<size_t>(1000, kept_count);
            std::uniform_int_distribution<size_t> uni_kept(0, kept_count - 1);
            for (size_t t = 0; t < n_kept_checks; ++t) {
                size_t k = uni_kept(rng);
                uint64_t rid = kept_tag[k].first;
                uint64_t stored_sa = static_cast<uint64_t>(sa_iv[k]);
                size_t pos = ltag.run_start_bwt(rid);
                uint64_t ref_sa = locate_sa_value_ref(r_index, pos);
                if (stored_sa != ref_sa) {
                    std::cerr << "  VERIFY FAIL (kept SA): rid=" << rid
                              << " pos=" << pos
                              << " stored=" << stored_sa
                              << " ref=" << ref_sa << "\n";
                    return 4;
                }
            }
            std::cerr << "  (a) kept-SA spot check: " << n_kept_checks << "/"
                      << n_kept_checks << " match\n";

            // (b) Simulate query for random deleted tag heads: LF walk from
            // t up to s steps, verify we land on a kept tag-head OR BWT-run
            // head whose SA + steps == SA[t]. Ground-truth SA[t] computed via
            // locate_sa_value_ref.
            //
            // We need a source of deleted tag heads. Iterate a random subset
            // of tag-run rids and skip any that are in kept_marker.
            std::uniform_int_distribution<size_t> uni_rid(0, num_tag_runs - 1);
            size_t n_del_checks_target = 100;
            size_t n_del_checked = 0;
            size_t n_del_attempts = 0;
            const size_t max_attempts = n_del_checks_target * 20;
            while (n_del_checked < n_del_checks_target && n_del_attempts < max_attempts) {
                size_t rid = uni_rid(rng);
                n_del_attempts++;
                if (kept_marker[rid]) continue;  // this rid is kept; skip
                size_t pos = ltag.run_start_bwt(rid);
                if (pos >= n) continue;  // sentinel
                uint64_t sa_true = locate_sa_value_ref(r_index, pos);

                // Simulate query: LF-walk from pos up to s steps.
                // Note: LF at step 0 is `pos` itself. Check whether pos is
                // a BWT-run head (which would mean rid is a coincident tag
                // head, expected). If not, walk LF.
                bool found = false;
                size_t cur = pos;
                uint64_t reconstructed = 0;
                for (size_t step = 0; step <= s; ++step) {
                    size_t bwt_rid = 0;
                    if (is_bwt_run_head_with_id(r_index, cur, bwt_rid)) {
                        uint64_t landed_sa = r_index.getSample(bwt_rid);
                        reconstructed = landed_sa + step;
                        found = true;
                        break;
                    }
                    size_t tag_rid = 0;
                    if (is_tag_run_head(ltag, cur, tag_rid) && kept_marker[tag_rid]) {
                        size_t slot = km_rank(tag_rid);
                        uint64_t landed_sa = static_cast<uint64_t>(sa_iv[slot]);
                        reconstructed = landed_sa + step;
                        found = true;
                        break;
                    }
                    if (step == s) break;  // don't LF past the cap
                    cur = r_index.LF(cur);
                }
                if (!found) {
                    std::cerr << "  VERIFY FAIL (unreachable): rid=" << rid
                              << " pos=" << pos << " SA_true=" << sa_true
                              << " no sample within " << s << " LF steps\n";
                    return 4;
                }
                if (reconstructed != sa_true) {
                    std::cerr << "  VERIFY FAIL (SA mismatch): rid=" << rid
                              << " pos=" << pos
                              << " SA_true=" << sa_true
                              << " reconstructed=" << reconstructed << "\n";
                    return 4;
                }
                n_del_checked++;
            }
            std::cerr << "  (b) deleted-head LF-walk: " << n_del_checked << "/"
                      << n_del_checks_target << " reachable within s and correct ("
                      << n_del_attempts << " random rids sampled to find "
                      << n_del_checked << " deleted)\n";
        }

        // Serialize.
        std::ostringstream out_name;
        out_name << out_base << ".tag_samples.s" << s;
        std::string out_path = out_name.str();
        {
            std::ofstream out(out_path, std::ios::binary | std::ios::trunc);
            if (!out) { std::cerr << "ERROR: cannot open " << out_path << "\n"; return 1; }
            out.write(THSAMP_MAGIC, sizeof(THSAMP_MAGIC));
            uint64_t ver = THSAMP_VERSION;
            uint64_t s_u64 = s;
            uint64_t ntr = num_tag_runs;
            uint64_t kc = kept_count;
            uint64_t nb = n;
            out.write(reinterpret_cast<const char*>(&ver),   sizeof(ver));
            out.write(reinterpret_cast<const char*>(&s_u64), sizeof(s_u64));
            out.write(reinterpret_cast<const char*>(&ntr),   sizeof(ntr));
            out.write(reinterpret_cast<const char*>(&kc),    sizeof(kc));
            out.write(reinterpret_cast<const char*>(&nb),    sizeof(nb));
            if (!out) { std::cerr << "ERROR: header write failed\n"; return 1; }
            kept_marker.serialize(out);
            sa_iv.serialize(out);
            if (!out) { std::cerr << "ERROR: body write failed\n"; return 1; }
        }
        long long file_sz = file_size_bytes(out_path);
        std::cerr << "  wrote:              " << out_path
                  << " (" << human_bytes(file_sz) << ", " << file_sz << " bytes)\n";

        auto t_pass_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> pass_dt = t_pass_end - t_pass_start;
        std::cerr << "  s=" << s << " pass time: " << pass_dt.count() << " s\n";
    }

    std::cerr << "\n----------------------------------------\n"
              << "Done.\n";
    return 0;
}
