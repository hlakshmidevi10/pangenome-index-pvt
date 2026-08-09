// bench_lf:
//   Microbenchmark for LF-operation latency on a FastLocate r-index.
//   Compares two implementations:
//     - LF()      : today's naive path (two independent block walks -- one for
//                   bwt_char_at_encoded, one for rankAt_encoded).
//     - LF_scan() : fused single scan_at() call producing both char and
//                   ranks in one predecessor + block-walk.
//
//   For each variant, measures two access patterns:
//     - Random  : N random BWT positions, one LF call each.  Models the
//                 per-fresh-walk cost that dominates when tag-head LF walks
//                 start at unrelated positions (sr-index-style SA recovery,
//                 first step per walk).
//     - Chain   : M random starts x K sequential LF steps (cur = LF(cur))
//                 per start.  Models the cost of steps 2..s in an ongoing
//                 walk from a deleted tag-head to a nearby sample.
//
//   Reports median-of-T trials plus stddev.  Cross-checks LF()==LF_scan()
//   for a verification subset before benching.
//
// Usage:
//   bench_lf <ri_file> [--n-random N] [--m-starts M] [--k-steps K]
//            [--trials T] [--seed S] [--verify V] [--skip-verify]
//
//   Defaults: N=1000000, M=10000, K=100, T=5, S=42, V=10000.

#include "pangenome_index/r-index.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace panindexer;

namespace {

void usage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " <ri_file> [options]\n"
        "\n"
        "  <ri_file>          Input r-index (.ri, must be JR-014 format).\n"
        "\n"
        "Options:\n"
        "  --n-random N       Random-access iterations (default 1,000,000)\n"
        "  --m-starts M       Chain-mode start positions   (default 10,000)\n"
        "  --k-steps K        Chain-mode steps per start   (default 100)\n"
        "  --trials T         Timed trials per variant     (default 5)\n"
        "  --seed S           RNG seed                     (default 42)\n"
        "  --verify V         Correctness cross-check pairs (default 10,000)\n"
        "  --skip-verify      Skip correctness check\n"
        "  -h, --help         Show this help\n";
}

double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    size_t n = v.size();
    if (n == 0) return 0.0;
    return (n % 2 == 1) ? v[n/2] : 0.5 * (v[n/2 - 1] + v[n/2]);
}

double stddev(const std::vector<double>& v) {
    if (v.size() < 2) return 0.0;
    double mean = 0.0;
    for (double x : v) mean += x;
    mean /= v.size();
    double ss = 0.0;
    for (double x : v) { double d = x - mean; ss += d * d; }
    return std::sqrt(ss / (v.size() - 1));
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) { usage(argv[0]); return 2; }

    std::string ri_path;
    size_t n_random = 1000000;
    size_t m_starts = 10000;
    size_t k_steps  = 100;
    size_t trials   = 5;
    uint64_t seed   = 42;
    size_t verify_n = 10000;
    bool skip_verify = false;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
        else if (a[0] != '-' && ri_path.empty()) { ri_path = a; }
        else if (a == "--n-random" && i + 1 < argc) { n_random = std::stoull(argv[++i]); }
        else if (a == "--m-starts" && i + 1 < argc) { m_starts = std::stoull(argv[++i]); }
        else if (a == "--k-steps"  && i + 1 < argc) { k_steps  = std::stoull(argv[++i]); }
        else if (a == "--trials"   && i + 1 < argc) { trials   = std::stoull(argv[++i]); }
        else if (a == "--seed"     && i + 1 < argc) { seed     = std::stoull(argv[++i]); }
        else if (a == "--verify"   && i + 1 < argc) { verify_n = std::stoull(argv[++i]); }
        else if (a == "--skip-verify") { skip_verify = true; }
        else { std::cerr << "Unknown option: " << a << "\n"; usage(argv[0]); return 2; }
    }
    if (ri_path.empty()) { usage(argv[0]); return 2; }

    std::cerr << "========================================\n"
              << "bench_lf\n"
              << "========================================\n"
              << "  ri:          " << ri_path << "\n"
              << "  n-random:    " << n_random << "\n"
              << "  m-starts:    " << m_starts << "\n"
              << "  k-steps:     " << k_steps << " (chain total = "
              << m_starts * k_steps << ")\n"
              << "  trials:      " << trials << "\n"
              << "  seed:        " << seed << "\n"
              << "  verify:      " << (skip_verify ? "skipped" : std::to_string(verify_n) + " pairs") << "\n"
              << "----------------------------------------\n";

    // Load r-index.
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
    const size_t n = r_index.bwt_size();
    if (n == 0) { std::cerr << "ERROR: empty BWT\n"; return 1; }
    std::cerr << "----------------------------------------\n";

    // -------------------------
    // Correctness verification
    // -------------------------
    if (!skip_verify) {
        std::cerr << "Verifying LF() == LF_scan() on " << verify_n << " random positions...\n";
        std::mt19937_64 rng(seed);
        std::uniform_int_distribution<size_t> uni(0, n - 1);
        size_t mismatches = 0;
        FastLocate::block_scan_result out;
        for (size_t i = 0; i < verify_n; ++i) {
            size_t pos = uni(rng);
            size_t a = r_index.LF(pos);
            size_t b = r_index.LF_scan(pos, out);
            if (a != b) {
                if (mismatches < 5) {
                    std::cerr << "  MISMATCH at pos=" << pos
                              << " LF=" << a << " LF_scan=" << b << "\n";
                }
                mismatches++;
            }
        }
        if (mismatches > 0) {
            std::cerr << "VERIFY FAIL: " << mismatches << " / " << verify_n
                      << " mismatches\n";
            return 3;
        }
        std::cerr << "  " << verify_n << " / " << verify_n << " match\n"
                  << "----------------------------------------\n";
    }

    // Generate the random position pool ONCE and use for all variants/trials,
    // so measurements are on identical inputs.
    std::cerr << "Generating random position pools...\n";
    std::vector<size_t> random_positions;
    random_positions.reserve(n_random);
    {
        std::mt19937_64 rng(seed ^ 0xB055uLL);
        std::uniform_int_distribution<size_t> uni(0, n - 1);
        for (size_t i = 0; i < n_random; ++i) random_positions.push_back(uni(rng));
    }

    std::vector<size_t> chain_starts;
    chain_starts.reserve(m_starts);
    {
        std::mt19937_64 rng(seed ^ 0xC0DEuLL);
        std::uniform_int_distribution<size_t> uni(0, n - 1);
        for (size_t i = 0; i < m_starts; ++i) chain_starts.push_back(uni(rng));
    }
    std::cerr << "  " << random_positions.size() << " random positions, "
              << chain_starts.size() << " chain starts (" << k_steps << " steps each)\n"
              << "----------------------------------------\n";

    // Warmup: one untimed pass of both variants over both patterns.
    std::cerr << "Warmup pass..." << std::endl;
    {
        size_t sink = 0;
        FastLocate::block_scan_result out;
        for (size_t p : random_positions) sink ^= r_index.LF(p);
        for (size_t p : random_positions) sink ^= r_index.LF_scan(p, out);
        for (size_t s : chain_starts) {
            size_t cur = s;
            for (size_t k = 0; k < k_steps; ++k) cur = r_index.LF(cur);
            sink ^= cur;
        }
        for (size_t s : chain_starts) {
            size_t cur = s;
            for (size_t k = 0; k < k_steps; ++k) cur = r_index.LF_scan(cur, out);
            sink ^= cur;
        }
        // Volatile sink to prevent DCE.
        volatile size_t sink_out = sink;
        (void)sink_out;
    }
    std::cerr << "  done\n----------------------------------------\n";

    // Run T trials of each (variant, pattern) combination and record wall.
    auto bench_random = [&](auto&& lf_fn, const char* label) {
        std::vector<double> ns_per_op(trials);
        for (size_t t = 0; t < trials; ++t) {
            size_t sink = 0;
            auto t0 = std::chrono::high_resolution_clock::now();
            for (size_t p : random_positions) sink ^= lf_fn(p);
            auto t1 = std::chrono::high_resolution_clock::now();
            volatile size_t sink_out = sink;
            (void)sink_out;
            double dt = std::chrono::duration<double, std::nano>(t1 - t0).count();
            ns_per_op[t] = dt / random_positions.size();
        }
        double med = median(ns_per_op);
        double sd  = stddev(ns_per_op);
        std::cout << std::left << std::setw(24) << label
                  << " random  n=" << random_positions.size()
                  << "  median=" << std::fixed << std::setprecision(2) << med << " ns/op"
                  << "  stddev=" << std::setprecision(2) << sd
                  << "  trials: ";
        for (double x : ns_per_op) std::cout << std::fixed << std::setprecision(1) << x << " ";
        std::cout << "\n";
    };

    auto bench_chain = [&](auto&& lf_fn, const char* label) {
        const size_t total_ops = chain_starts.size() * k_steps;
        std::vector<double> ns_per_op(trials);
        for (size_t t = 0; t < trials; ++t) {
            size_t sink = 0;
            auto t0 = std::chrono::high_resolution_clock::now();
            for (size_t s : chain_starts) {
                size_t cur = s;
                for (size_t k = 0; k < k_steps; ++k) cur = lf_fn(cur);
                sink ^= cur;
            }
            auto t1 = std::chrono::high_resolution_clock::now();
            volatile size_t sink_out = sink;
            (void)sink_out;
            double dt = std::chrono::duration<double, std::nano>(t1 - t0).count();
            ns_per_op[t] = dt / total_ops;
        }
        double med = median(ns_per_op);
        double sd  = stddev(ns_per_op);
        std::cout << std::left << std::setw(24) << label
                  << " chain   n=" << total_ops
                  << "  median=" << std::fixed << std::setprecision(2) << med << " ns/op"
                  << "  stddev=" << std::setprecision(2) << sd
                  << "  trials: ";
        for (double x : ns_per_op) std::cout << std::fixed << std::setprecision(1) << x << " ";
        std::cout << "\n";
    };

    std::cout << "\n=== Benchmark results ===\n";
    FastLocate::block_scan_result reusable_out;
    bench_random([&](size_t p){ return r_index.LF(p); },                             "LF (baseline)");
    bench_random([&](size_t p){ return r_index.LF_scan(p, reusable_out); },          "LF_scan (fused)");
    bench_chain ([&](size_t p){ return r_index.LF(p); },                             "LF (baseline)");
    bench_chain ([&](size_t p){ return r_index.LF_scan(p, reusable_out); },          "LF_scan (fused)");

    std::cout << "\nDone.\n";
    return 0;
}
