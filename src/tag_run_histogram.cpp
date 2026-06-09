// tag_run_histogram: stream every tag run in a *_compressed.tags file and
// emit summary stats + a log2-bucketed run-length histogram.
//
// Output format (stdout, all to make machine-parsing trivial):
//   #SUMMARY  total_runs  total_length  min  max  mean
//   #HIST     bucket_lo   bucket_hi     count  cumulative_count  cumulative_fraction
//
// Percentiles are intentionally NOT computed here; derive them post-hoc from
// the cumulative_count column in the #HIST rows. The histogram is exact at
// bucket edges, so for log2 buckets we already have ~64 percentile-friendly
// breakpoints across the full range.
//
// Memory: O(1) per pass; runtime: O(tag_runs) -- a few minutes on HPRC chr6.
//
// Usage: tag_run_histogram <compressed_tags.tags>

#include <iostream>
#include <fstream>
#include <string>
#include <cstdint>
#include <array>
#include <limits>

#include "pangenome_index/tag_arrays.hpp"

using namespace std;
using namespace panindexer;

static void usage(const char* prog) {
    cerr << "Usage: " << prog << " <compressed_tags.tags>\n"
         << "  Streams every tag run, emits summary + log2-bucketed run-length histogram on stdout.\n";
}

int main(int argc, char** argv) {
    if (argc != 2) {
        usage(argv[0]);
        return 1;
    }
    const string tags_path = argv[1];

    TagArray tag_array;
    {
        ifstream tin(tags_path, ios::binary);
        if (!tin) {
            cerr << "Error: cannot open tags file: " << tags_path << "\n";
            return 1;
        }
        // Use the _compact alias (== _sdsl loader) for consistency with find_mems
        // and the other runtime tools, all of which read the same on-disk format.
        tag_array.load_compressed_tags_compact(tin);
    }

    // 64 log2 buckets: bucket b covers [2^b, 2^(b+1)).
    // bucket 0 = [1,1], bucket 1 = [2,3], bucket 2 = [4,7], ..., bucket 63 = [2^63, 2^64-1].
    // Run lengths are stored as uint64_t in the encoded format; in practice
    // all values will land in buckets <= 40 (≈ 1 trillion) but we size for the type.
    constexpr size_t NUM_BUCKETS = 64;
    array<uint64_t, NUM_BUCKETS> bucket_count{};

    uint64_t total_runs = 0;
    uint64_t total_length = 0;
    uint64_t min_len = numeric_limits<uint64_t>::max();
    uint64_t max_len = 0;

    // Single streaming pass via the public callback iterator.
    // We don't need pos_t here -- only the run length.
    tag_array.for_each_run_compact([&](handlegraph::pos_t /*pos*/, uint64_t run_length) {
        total_runs++;
        total_length += run_length;
        if (run_length < min_len) min_len = run_length;
        if (run_length > max_len) max_len = run_length;

        // log2 bucket: index of the highest set bit. __builtin_clzll is defined
        // for non-zero inputs; run_length should always be >= 1 for a real run.
        if (run_length == 0) {
            // Defensive: a zero-length run would indicate a build bug.
            // Fold it into bucket 0 with a stderr note rather than crash.
            bucket_count[0]++;
            cerr << "warning: encountered tag run with length=0\n";
        } else {
            unsigned b = 63u - static_cast<unsigned>(__builtin_clzll(run_length));
            bucket_count[b]++;
        }
    });

    if (total_runs == 0) {
        cerr << "Error: no tag runs found in " << tags_path << "\n";
        return 2;
    }

    const double mean = static_cast<double>(total_length) / static_cast<double>(total_runs);

    // ---- Output ------------------------------------------------------------
    cout << "#SUMMARY\ttotal_runs\ttotal_length\tmin\tmax\tmean\n";
    cout << "#SUMMARY\t"
         << total_runs   << "\t"
         << total_length << "\t"
         << min_len      << "\t"
         << max_len      << "\t"
         << mean         << "\n";

    // Stop emitting once we pass the bucket containing max_len -- the rest
    // are guaranteed zero and clutter the output.
    const unsigned last_nonzero_bucket =
        63u - static_cast<unsigned>(__builtin_clzll(max_len));

    cout << "#HIST\tbucket_lo\tbucket_hi\tcount\tcumulative_count\tcumulative_fraction\n";
    uint64_t cum = 0;
    for (size_t b = 0; b <= last_nonzero_bucket; b++) {
        cum += bucket_count[b];
        const uint64_t lo = (b == 0) ? 1ULL : (1ULL << b);
        const uint64_t hi = (b == 63) ? numeric_limits<uint64_t>::max()
                                      : ((1ULL << (b + 1)) - 1ULL);
        const double frac = static_cast<double>(cum) / static_cast<double>(total_runs);
        cout << "#HIST\t" << lo << "\t" << hi << "\t" << bucket_count[b]
             << "\t" << cum << "\t" << frac << "\n";
    }

    return 0;
}
