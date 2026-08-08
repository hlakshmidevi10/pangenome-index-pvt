// build_tag_head_samples:
//   Build a sampled SA table over tag-run heads, using an sr-index-inspired
//   deletion rule to keep only a subset of samples while guaranteeing that
//   any tag-run head is at most s LF steps from either a kept tag-run-head
//   sample OR a BWT-run head (which are always sampled in the r-index).
//
// Deletion rule (per user spec, in BWT-position space, single left-to-right pass):
//   Enumerate merged candidates (tag-run heads + BWT-run heads) in BWT-position
//   order. BWT-run heads are unremovable anchors. For each tag-run head t_i,
//   delete it iff (next_candidate.pos - prev_kept.pos < s).
//
// Output file format (per s value): <output_base>.tag_samples.s<s>
//   [8B]  magic   = "THSAMP\0\0"
//   [8B]  version = 1
//   [8B]  s (stride threshold, in BWT positions)
//   [8B]  num_tag_runs   (from LightTagIndex, for consistency check)
//   [8B]  kept_count     (equal to number of 1-bits in kept_marker; also length of sa_values)
//   [8B]  bwt_size       (n)
//   [...] sd_vector<>::serialize -- kept_marker over tag-run rids
//   [...] int_vector<0>::serialize -- packed SA values, length = kept_count
//
// Runtime notes:
//   - Enumerating BWT-run heads by decoding blocks in order: cheap (O(num_bwt_runs)).
//   - Enumerating tag-run heads via LightTagIndex::run_start_bwt: cheap (O(num_tag_runs)).
//   - Deletion scan: single linear pass, O(num_tag_runs + num_bwt_runs).
//   - SA computation for kept samples: walks locateNext through each BWT run
//     up to the last kept tag-run head inside that run. Worst case n = BWT
//     size total locateNext calls; expected much less. On HPRC chr6 (n ~ 3e10)
//     this is the bottleneck -- expect 30-120 min per s on Vesuvio.

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
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace panindexer;

namespace {

constexpr char     THSAMP_MAGIC[8] = {'T','H','S','A','M','P','\0','\0'};
constexpr uint64_t THSAMP_VERSION  = 1;

void usage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " <ri_file> <ltags_file> <output_base> "
        "[--s LIST] [--limit-bwt N]\n"
        "\n"
        "  <ri_file>       Input r-index (.ri)\n"
        "  <ltags_file>    Input light tag index (.ltags)\n"
        "  <output_base>   Output filename prefix; writes <base>.tag_samples.s<N> per s\n"
        "\n"
        "Options:\n"
        "  --s LIST        Comma-separated s values (default: 32,64,128,256)\n"
        "  --limit-bwt N   Stop processing after BWT position N (for smoke tests)\n"
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

// Human-readable byte size.
std::string human_bytes(std::size_t bytes) {
    const char* units[] = {"B", "KB", "MB", "GB", "TB"};
    double v = static_cast<double>(bytes);
    int u = 0;
    while (v >= 1024.0 && u < 4) { v /= 1024.0; ++u; }
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << v << " " << units[u];
    return oss.str();
}

// Enumerator over BWT-run heads (start BWT position, incrementing run id).
// Emits (bwt_pos, bwt_run_id) in BWT-position order by iterating blocks.
struct BwtRunHeadIter {
    const FastLocate* r_index;
    size_t num_blocks;
    size_t block_idx = 0;
    size_t local_run_idx_in_block = 0;
    size_t global_run_id = 0;
    size_t cur_block_start = 0;
    size_t cur_offset_in_block = 0;
    std::vector<std::pair<size_t, size_t>> block_runs; // (symbol, length) pairs
    bool done_flag = false;

    explicit BwtRunHeadIter(const FastLocate& ri) : r_index(&ri) {
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

    // Current head position and run id.
    size_t pos() const { return cur_block_start + cur_offset_in_block; }
    size_t run_id() const { return global_run_id; }

    void advance() {
        if (done_flag) return;
        // Advance past current run.
        cur_offset_in_block += block_runs[local_run_idx_in_block].second;
        local_run_idx_in_block++;
        global_run_id++;
        if (local_run_idx_in_block >= block_runs.size()) {
            block_idx++;
            load_current_block();
        }
    }
};

} // namespace

int main(int argc, char** argv) {
    if (argc < 4) { usage(argv[0]); return 2; }

    std::string ri_path = argv[1];
    std::string ltags_path = argv[2];
    std::string out_base = argv[3];
    std::vector<size_t> s_values = {32, 64, 128, 256};
    size_t limit_bwt = 0; // 0 means no limit

    for (int i = 4; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
        else if (a == "--s" && i + 1 < argc) {
            s_values = parse_s_list(argv[++i]);
        }
        else if (a == "--limit-bwt" && i + 1 < argc) {
            limit_bwt = static_cast<size_t>(std::stoll(argv[++i]));
        }
        else {
            std::cerr << "Unknown option: " << a << "\n";
            usage(argv[0]);
            return 2;
        }
    }

    std::cerr << "========================================\n"
              << "build_tag_head_samples\n"
              << "========================================\n"
              << "  ri:          " << ri_path << "\n"
              << "  ltags:       " << ltags_path << "\n"
              << "  output base: " << out_base << "\n"
              << "  s values:   ";
    for (auto s : s_values) std::cerr << " " << s;
    std::cerr << "\n";
    if (limit_bwt > 0) std::cerr << "  --limit-bwt: " << limit_bwt << "\n";
    std::cerr << "----------------------------------------\n";

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
                  << "  load time:     " << dt.count() << " s\n";
    }

    // --- Load ltags ---
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

    // ltag.bwt_size() is bwt_intervals.size(), which is n+1 (one extra
    // trailing position — see mem-projection/pangenome-pipeline/CLAUDE.md
    // "bwt_intervals size == n+1"). Accept either exact match or off-by-one.
    if (ltag.bwt_size() != r_index.bwt_size() && ltag.bwt_size() != r_index.bwt_size() + 1) {
        std::cerr << "ERROR: BWT size mismatch between r-index and ltags ("
                  << r_index.bwt_size() << " vs " << ltag.bwt_size()
                  << "; expected equal or +1)\n";
        return 1;
    }
    const size_t n = r_index.bwt_size();
    const size_t num_tag_runs = ltag.num_runs();
    const size_t num_bwt_runs = r_index.tot_runs();
    const size_t effective_n = (limit_bwt > 0 && limit_bwt < n) ? limit_bwt : n;

    std::cerr << "----------------------------------------\n";

    // Process each s value independently. Each pass:
    //   1. Enumerate candidates in BWT-position order (merge scan).
    //   2. Apply deletion rule; record kept tag-run rids and their (bwt_pos, offset within
    //      containing BWT run) for the SA-walk phase.
    //   3. Walk BWT runs via locateNext, computing SA for each kept tag-run head.
    //   4. Build sd_vector kept_marker + int_vector sa_values; serialize.
    //
    // Memory: we store `kept_rids` (uint64 each) + `kept_bwt_pos` (uint64 each) as
    // in-memory buffers. At kept ~500M, that's ~8 GB per buffer = 16 GB. On
    // machines where that's a problem, we could stream directly to disk, but
    // for a one-off analysis on Vesuvio's 1 TB RAM this is fine.

    for (size_t s : s_values) {
        std::cerr << "\n=========================================\n"
                  << "s = " << s << "\n"
                  << "-----------------------------------------\n";

        auto t_pass_start = std::chrono::high_resolution_clock::now();

        // Buffers to collect kept tag-run heads.
        std::vector<uint64_t> kept_rids;
        std::vector<uint64_t> kept_bwt_pos;
        // Optimistic reserve: assume up to 100% kept.
        // Users with tight RAM can subset with --limit-bwt to smoke-test first.
        kept_rids.reserve(num_tag_runs / 4);
        kept_bwt_pos.reserve(num_tag_runs / 4);

        // Deletion-rule scan: enumerate merged sorted candidates. To decide
        // whether to keep tag-run head t_i, we need next_candidate.pos (peek).
        // We maintain a 3-slot lookahead over the tag-run head stream and
        // BWT-run head stream, always knowing "curr" and "next".
        auto t_scan_start = std::chrono::high_resolution_clock::now();

        BwtRunHeadIter bwt_iter(r_index);
        size_t next_tag_rid = 0;
        // The ltag bwt_intervals has one extra trailing 1-bit at position n
        // (sentinel). Skip any tag-run head that would land at pos >= n.
        auto next_tag_pos_fn = [&]() -> size_t {
            while (next_tag_rid < num_tag_runs) {
                size_t p = ltag.run_start_bwt(next_tag_rid);
                if (p < n) return p;
                // Skip out-of-range sentinel.
                next_tag_rid++;
            }
            return SIZE_MAX;
        };

        size_t prev_kept_pos = 0;
        // The very first candidate at BWT pos 0 is always a BWT-run start
        // (either the endmarker run or the first run of the smallest character).
        // We initialize prev_kept_pos = 0 to reflect the mandatory BWT anchor
        // at position 0; subsequent scan starts from candidates >= 1.
        //
        // But we need to handle general enumeration. We'll do a proper merge.
        bool have_prev = false;

        auto candidate_pos = [&](bool prefer_bwt_on_tie) -> size_t {
            size_t bp = bwt_iter.done() ? SIZE_MAX : bwt_iter.pos();
            size_t tp = next_tag_pos_fn();
            if (bp == SIZE_MAX && tp == SIZE_MAX) return SIZE_MAX;
            if (bp == SIZE_MAX) return tp;
            if (tp == SIZE_MAX) return bp;
            if (bp < tp) return bp;
            if (tp < bp) return tp;
            // Equal positions: shouldn't happen because BWT-run starts and
            // tag-run starts should coincide (a tag run always aligns with a
            // BWT-run boundary in the build), in which case they're the same
            // structural position. To be safe, prefer bwt to avoid double-emit.
            return bp;
        };

        // Merge scan: at each step, take the smaller of (bwt_iter.pos, next_tag_pos).
        // If they're equal, treat as a single position; the BWT-run head keeps and
        // we skip the coincident tag-run head (delete-neutral: we would emit an SA
        // for the tag-run head elsewhere anyway via the BWT-anchor's sample).
        //
        // Actually, per user's design: if a tag-run head coincides with a BWT-run
        // head, its SA is already available for free via samples[rid], so the
        // "kept tag-run heads" set can skip it -- there's nothing to sample.
        // But conceptually we might want to note this so a query can shortcut.
        // For storage counting, we'll consider such coincident tag-run heads
        // as always-deleted (no sample needed).

        while (true) {
            size_t bp = bwt_iter.done() ? SIZE_MAX : bwt_iter.pos();
            size_t tp = next_tag_pos_fn();
            size_t curr = std::min(bp, tp);
            if (curr == SIZE_MAX) break;
            if (limit_bwt > 0 && curr >= limit_bwt) break;

            if (bp == tp) {
                // Coincident BWT-anchor + tag-run head. Treat as anchor.
                // Tag-run head at this position is implicitly deleted (SA free from anchor).
                prev_kept_pos = curr;
                have_prev = true;
                bwt_iter.advance();
                next_tag_rid++;
                continue;
            }
            if (bp < tp) {
                // BWT anchor. Always keep.
                prev_kept_pos = curr;
                have_prev = true;
                bwt_iter.advance();
                continue;
            }
            // tp < bp: tag-run head candidate.
            // Look at next candidate after curr = tp (advancing tag iter mentally).
            size_t bp_after = bp; // bwt_iter still points to bp, which is >= tp
            size_t tp_after = SIZE_MAX;
            // Look ahead over any trailing sentinel entries.
            {
                size_t look = next_tag_rid + 1;
                while (look < num_tag_runs) {
                    size_t p = ltag.run_start_bwt(look);
                    if (p < n) { tp_after = p; break; }
                    look++;
                }
            }
            size_t next_after = std::min(bp_after, tp_after);

            bool keep = true;
            if (have_prev && next_after != SIZE_MAX) {
                if (next_after - prev_kept_pos < s) keep = false;
            }
            if (keep) {
                kept_rids.push_back(next_tag_rid);
                kept_bwt_pos.push_back(tp);
                prev_kept_pos = tp;
                have_prev = true;
            }
            next_tag_rid++;
        }

        auto t_scan_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> scan_dt = t_scan_end - t_scan_start;

        const size_t kept_count = kept_rids.size();
        std::cerr << "  candidates scanned: tag_runs<=" << num_tag_runs
                  << ", bwt_runs<=" << num_bwt_runs << "\n"
                  << "  kept tag-run heads: " << kept_count
                  << " / " << num_tag_runs
                  << " (" << std::fixed << std::setprecision(2)
                  << 100.0 * kept_count / (num_tag_runs ? num_tag_runs : 1)
                  << "% kept, " << 100.0 * (num_tag_runs - kept_count) / (num_tag_runs ? num_tag_runs : 1)
                  << "% deleted)\n"
                  << "  deletion scan time: " << scan_dt.count() << " s\n";

        // Now compute SA for each kept tag-run head. Walk BWT runs in order
        // via locateNext, maintaining a running SA cursor initialized at each
        // BWT-run's sample. When we cross a kept tag-run head, record SA.

        // Group kept tag heads by containing BWT run: we need to walk each BWT
        // run's positions with an SA cursor. We know kept_bwt_pos is already
        // sorted (we appended during a left-to-right scan). We'll iterate BWT
        // runs in tandem with a pointer into kept_bwt_pos.

        auto t_sa_start = std::chrono::high_resolution_clock::now();

        std::vector<uint64_t> sa_values(kept_count, 0);
        size_t kept_idx = 0;
        size_t total_locate_next = 0;

        BwtRunHeadIter bwt_iter2(r_index);
        // We also need each BWT run's length. get_block_runs gives (symbol, length)
        // per run in a block, so let's iterate blocks and process runs manually.
        // Since we're going through BwtRunHeadIter which already tracks the block
        // and per-run offsets, we can leverage it: at each run start, we know the
        // run's length is block_runs[local_run_idx_in_block].second.

        size_t bwt_run_id_counter = 0;
        (void)bwt_run_id_counter;

        while (!bwt_iter2.done() && kept_idx < kept_count) {
            size_t run_bwt_start = bwt_iter2.pos();
            size_t run_len = bwt_iter2.block_runs[bwt_iter2.local_run_idx_in_block].second;
            size_t run_bwt_end = run_bwt_start + run_len; // exclusive
            size_t bwt_run_id = bwt_iter2.run_id();

            if (limit_bwt > 0 && run_bwt_start >= limit_bwt) break;

            // Any kept tag-run heads inside [run_bwt_start, run_bwt_end)?
            if (kept_bwt_pos[kept_idx] >= run_bwt_end) {
                // Nothing to sample in this BWT run, skip.
                bwt_iter2.advance();
                continue;
            }

            // Walk this BWT run with locateNext.
            FastLocate::size_type sa_cursor = r_index.getSample(bwt_run_id);
            size_t cur_pos = run_bwt_start;
            while (kept_idx < kept_count && kept_bwt_pos[kept_idx] < run_bwt_end) {
                size_t target = kept_bwt_pos[kept_idx];
                while (cur_pos < target) {
                    sa_cursor = r_index.locateNext(sa_cursor);
                    total_locate_next++;
                    cur_pos++;
                }
                sa_values[kept_idx] = static_cast<uint64_t>(sa_cursor);
                kept_idx++;
            }

            bwt_iter2.advance();
        }

        auto t_sa_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> sa_dt = t_sa_end - t_sa_start;
        std::cerr << "  SA computation:     " << sa_dt.count()
                  << " s (" << total_locate_next << " locateNext calls, "
                  << (total_locate_next ? sa_dt.count() * 1e6 / total_locate_next : 0.0)
                  << " us/call)\n";

        if (kept_idx != kept_count) {
            std::cerr << "  WARNING: SA computation covered only " << kept_idx
                      << " / " << kept_count << " kept samples "
                      << "(possibly truncated by --limit-bwt or a BWT/tag index mismatch)\n";
        }

        // --- Build sd_vector kept_marker (over tag-run rids) + int_vector sa_values ---
        auto t_build_start = std::chrono::high_resolution_clock::now();

        // sd_vector from sorted positions:
        // Use sd_vector_builder for O(kept) construction.
        sdsl::sd_vector<> kept_marker;
        {
            sdsl::sd_vector_builder builder(num_tag_runs, kept_count);
            for (uint64_t rid : kept_rids) builder.set(rid);
            kept_marker = sdsl::sd_vector<>(builder);
        }
        // Free rid buffer.
        std::vector<uint64_t>().swap(kept_rids);

        // Pack SA values into int_vector with bit-width = ceil(log2(n+1))
        const size_t bits_per_sa = (n <= 1) ? 1 : static_cast<size_t>(std::ceil(std::log2(static_cast<double>(n + 1))));
        sdsl::int_vector<0> sa_iv(kept_count, 0, bits_per_sa);
        for (size_t i = 0; i < kept_count; ++i) sa_iv[i] = sa_values[i];
        std::vector<uint64_t>().swap(sa_values);
        std::vector<uint64_t>().swap(kept_bwt_pos);

        auto t_build_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> build_dt = t_build_end - t_build_start;

        const size_t sd_bytes = sdsl::size_in_bytes(kept_marker);
        const size_t iv_bytes = sdsl::size_in_bytes(sa_iv);

        std::cerr << "  bits/sa:            " << bits_per_sa << "\n"
                  << "  kept_marker (sdv):  " << human_bytes(sd_bytes)
                  << " (" << sd_bytes << " bytes, "
                  << (kept_count ? sd_bytes * 8.0 / kept_count : 0.0) << " bits/kept)\n"
                  << "  sa_values (int_v):  " << human_bytes(iv_bytes)
                  << " (" << iv_bytes << " bytes)\n"
                  << "  build time:         " << build_dt.count() << " s\n";

        // --- Serialize ---
        std::ostringstream out_name;
        out_name << out_base << ".tag_samples.s" << s;
        std::string out_path = out_name.str();
        {
            std::ofstream out(out_path, std::ios::binary | std::ios::trunc);
            if (!out) { std::cerr << "ERROR: cannot open " << out_path << " for writing\n"; return 1; }
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
        std::cerr << "  s=" << s << " total pass time: " << pass_dt.count() << " s\n";
    }

    std::cerr << "\n----------------------------------------\n"
              << "Done.\n";
    return 0;
}
