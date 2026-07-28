# Findings: Seed cost dominates locate work in `find_mems` on long BWT runs

**Branch:** `classify-mem-runs` (based on `find-mems-all-positions` @ `22d9991`)
**Date:** 2026-07-23
**Question posed:** Where is the optimization headroom in `find_mems`'s locate hot path? Specifically, does the sr-index-style `samples[rid]` anchor help skip work between consecutive tag emissions inside a MEM?

## TL;DR

- **`samples[rid]` optimization is dead across all datasets tested.** Savings pin at 0.26-0.27% of walk cost regardless of dataset (yeast n/r=2, HPRCv1 chr6 short-MEM sample, HPRCv1 chr6 full 100k reads). Do not implement.
- **The real target is seed cost.** 32.2% of total walk cost on HPRCv1 chr6 100k reads (34% on yeast). Untouched by `samples[rid]` by construction.
- **Seed cost concentrates in the long-BWT-run tail.** 8.8% of MEMs (those with `avg_bwt_run_length >= 1000` bp) contribute **75.4% of all seed cost** and 24.3% of total walk work.
- **Proposed next optimization (not yet built):** SA subsample inside long BWT runs. Projected ~7% end-to-end `find_mems` wall reduction, ~3-5 MB additional `.ri` storage per chromosome.

## Background: what walk cost is

`dump_mem_info_lightweight` (`src/find_mems.cpp:604`) processes each MEM as follows:

1. **Seed:** call `locate_sa_value(sp)` which finds the BWT run containing `sp`, reads the SA sample stored at the run's start (`samples[run_id]`), and walks `locateNext` forward `sp - run_start_bwt` times. This is the **seed cost**.
2. **Emit** an SA position at the first tag run whose start is <= `sp`.
3. **Walk forward** by `locateNext`, emitting at each subsequent tag-run start, until past the last tag-run start inside `[sp, ep]`. Total forward walk = `last_tag_start - sp` calls. This is the **inter-tag cost**.

Each `locateNext` call costs ~0.24 microseconds (empirically, from prior `FINDINGS_PERF` measurements). Walk cost is therefore a direct wall-time proxy.

Two algorithms are modeled by the classifier:

- **Current:** cost = `seed_cost + (last_tag_start - sp)`.
- **Optimized (samples[rid]):** if any BWT-run-start `b` falls between two consecutive tag emissions `t_{k-1}` and `t_k`, teleport to `samples[rid_of_b]` (O(1) array read) and walk `locateNext` from `b` to `t_k` instead of from `t_{k-1}` to `t_k`. Savings per MEM = `sum_k max(0, b_rightmost_in_gap - t_{k-1})` where `b_rightmost_in_gap` is the largest BWT-run-start in `(t_{k-1}, t_k]`. Seed cost is identical in both algorithms.

## Methodology

Added `--debug-classify=<path.tsv>` to `find_mems` (see `src/find_mems.cpp`, commit `40c7d5c`). For each MEM it emits a 20-column row (schema in `xy-test/verify_partition.awk` comments):

- Identity: `read_id, mem_start, mem_length, bwt_start, bwt_size`
- BWT/tag geometry: `num_bwt_runs, num_tag_runs, avg_bwt_run_length, avg_tags_per_bwt_run`
- Tag partition (priority-ordered, mutually exclusive): `tags_contain_bwt_start, tags_at_bwt_end, tags_strictly_interior, tags_other`
- Anchoring: `tags_span_multi_bwt`
- Walk cost: `seed_cost, inter_current, inter_optimized, current_walk_cost, optimized_walk_cost, savings`

Datasets used:

| Dataset                     | Reads    | MEMs     | Notes                                              |
|-----------------------------|---------:|---------:|----------------------------------------------------|
| Yeast chrII 100 kb, n/r=2   | 500      | 581      | Local smoke; only bucket where n/r is tiny         |
| HPRCv1 chr6 (sample)        | 1,000    | 1,534    | Sanity check at intermediate scale                 |
| HPRCv1 chr6 (full)          | 100,000  | 100,899  | Production-representative                          |

Diagnostic-only; flag off by default; no measurable wall impact when disabled.

## Result 1: `samples[rid]` savings pin at ~0.27% across three regimes

| Dataset                     | seed / walk | savings / walk |
|-----------------------------|------------:|---------------:|
| Yeast chrII 100 kb, n/r=2   |       33.8% |         0.269% |
| HPRCv1 chr6 sample (1k)     |       80.1% |         0.271% |
| HPRCv1 chr6 full (100k)     |       32.2% |         0.257% |

Three datasets, three different `seed/walk` ratios, but the `samples[rid]` savings ratio is essentially invariant at ~0.27%. This is the geometry of tag-vs-BWT-run interleaving speaking: it does not matter how you scale the workload, the fraction of gaps between consecutive tag emissions that happen to straddle a BWT-run-start is a stable ~0.3%.

End-to-end projection: `0.27% * 93% (locate share of MEM-proc) * 34% (MEM-proc share of find_mems) = 0.08% wall reduction`. This is deep in the noise floor. **Not worth building.**

## Result 2: seed cost concentrates in long BWT runs

HPRCv1 chr6 full run, seed cost bucketed by MEM's `avg_bwt_run_length`:

```
avg_bwt_len   MEMs     %MEMs    seed_cost      %seed   avg/MEM   %curwalk
<10             148    0.15%          166      0.00%       1.1     0.00%
10-100       44,800   44.40%      349,457      1.54%       7.8     0.50%
100-1000     47,060   46.64%    5,235,603     23.08%     111.3     7.44%
1000-10000    8,319    8.24%   12,912,951     56.92%    1552.2    18.34%
10000+          572    0.57%    4,186,751     18.46%    7319.5     5.95%
                                                                --------
                                                                  32.22%
```

**Key readings:**

- **8.8% of MEMs (the two rightmost buckets) contribute 75.4% of all seed cost** and 24.3% of total walk work.
- The `10000+` bucket has only 572 MEMs (0.57% of MEMs) but contributes 18.5% of seed cost. Average seed = 7,319 `locateNext` calls per MEM (~1.75 ms of pure seed work each at 0.24 microseconds/call).
- Seed cost per MEM scales roughly linearly with BWT run length, consistent with `E[seed] = L/2` when `sp` lands uniformly in a run of length `L`. In the `1000-10000` bucket, avg seed = 1552 implies `E[L] ~ 3100`, which sits at the geometric center of the window. Model consistent.

The `100-1000` bucket (46.6% of MEMs) contributes only 23% of seed cost, a proportionality gap that confirms the effect is real: seed cost does not scale with MEM count, it scales with the BWT volume that MEMs fall into.

## Density and partition (context)

**HPRCv1 chr6 full:**

- 100,899 MEMs, 21.77M tag runs, 182,171 BWT runs
- Volume-weighted avg BWT run length: **552 bp**
- Aggregate tags per BWT run: **119.5**
- `num_bwt_runs == 1` for 98.5% of MEMs; the remaining 1.5% span 2+ runs (0.2% span 65+ runs)

**Tag partition** (each tag run assigned to exactly one bucket, priority-ordered):

| Bucket                                       | Count      | %       |
|----------------------------------------------|-----------:|--------:|
| b1: `contains >=1 bwt_run_start` (anchor)    | 130,124    | 0.598%  |
| b2: `at_bwt_end`                             | 54,992     | 0.253%  |
| b3: `strictly_interior`                      | 21,588,494 | 99.150% |
| b4: `other`                                  | 53         | 0.000%  |

Bucket b1 is the population that a `samples[rid]` optimization could ever help; it is 0.6% of tag runs. Explains the 0.27% savings ceiling directly.

## Structural obstacle for a related idea (SA-during-backward-search)

An earlier line of investigation asked whether SA values could be carried through the FMD backward search itself (like MONI's approach), avoiding the seed cost entirely. This does not work for our indexer:

- `forward_extend_encoded` (`src/r-index.cpp:787-793`) is implemented as swap -> backward-extend-on-RC -> swap. The returned `bint.forward` is computed by summing RC-side rank differences (`src/r-index.cpp:770-777`), not via an LF walk.
- The identity `SA[LF(k)] = SA[k] - 1` holds only for `backward_extend_encoded`, not for `forward_extend_encoded`. And the three-phase MEM finder (`include/pangenome_index/algorithm.hpp:655-738`) mixes both.
- The permutation between the forward and RC BWT indices is not tracked for `size > 1` intervals, which is our workload.

So MONI-style SA-carrying is a non-starter without redesigning the FMD extension primitives. Seed cost has to be attacked at `locate_sa_value` itself, not upstream.

## Proposed next optimization (not built)

**"SA subsample inside long BWT runs."** Store additional SA samples at fixed strides inside BWT runs above a length threshold.

- Sampling: for BWT runs with length >= `min_run_len` (proposed: 1024), store one extra SA sample every `stride` BWT positions (proposed: 128).
- Seed path: `locate_sa_value(sp)` first checks whether an extra sample exists at or before `sp` inside `sp`'s BWT run; if so, jump to it (O(1) predecessor query on a small sparse bitvector per long run, or a flat sorted array), then `locateNext` walk at most `stride` positions.

**Projected seed cost reduction (back-of-envelope):**

- `1000-10000` bucket: current avg seed 1552 -> capped at ~64 (stride/2) -> saves ~96% of that bucket's seed cost = ~12.4M `locateNext` calls
- `10000+` bucket: current avg 7320 -> capped at ~64 -> saves ~99% = ~4.1M
- Total tail savings: ~16.5M of 22.7M seed_cost = **~73% seed cost reduction**
- Seed was 32.2% of walk -> walk drops by `0.322 * 0.73 = 23.5%`
- End-to-end: `23.5% * 93% * 34% ~ 7.4%` wall reduction for `find_mems`

Two orders of magnitude better than `samples[rid]`.

**Projected storage overhead:** rough upper bound. HPRCv1 chr6 has 182K BWT runs and ~100 Mbp of BWT volume (from `runs * volume-weighted-avg = 182171 * 552`). If ~50% of that volume sits in runs >= 1024 (needs measurement to confirm), that is 50 Mbp / 128 stride = ~390K extra samples = ~3 MB at 8 bytes/sample. Current samples table is ~1.4 MB. So maybe 3-5 MB additional `.ri` size per chromosome. Small.

**Unknowns to resolve before building:**

1. Exact BWT volume in runs >= `min_run_len` (add one more awk aggregation column to nail the storage projection).
2. Whether the `locateNext` calls in the seed walk are cache-cold or partially warm. If cold (predecessor query hits a different `blocks_start_pos` block each MEM), the wall-time saving may be larger than the walk-count model predicts.
3. The 0.24 microseconds/call figure comes from an amortized measurement mixing seed and inter-tag walks. Seed walks are likely more cache-cold. A calibration run comparing predicted vs measured wall reduction for a simple experiment would tighten the projection.

## What was validated but is not directly on the critical path

- Papers reviewed: sr-index (arXiv 2103.15329, 2409.14654) and MONI (Boucher et al., CMB 2021, 10.1089/cmb.2021.0290-2). MONI does SA-during-search but is unidirectional; sr-index is subsampled r-index but does not fold SA into the search. Neither directly transfers.
- Partition invariant `b1 + b2 + b3 + b4 == num_tag_runs` verified on yeast (0 violations across 34,332 tags) and HPRCv1 100k (b4 = 53 tags is a legitimate residual, not a bug).
- Bucket 1 semantics were broadened from "tag start equals BWT run start" to "tag contains any BWT run start" during development (commit `9221d4f`). The broader definition matches the applicability condition of the `samples[rid]` optimization and correctly folds in tags that equal or span a BWT run.

## Relevant code

- `src/find_mems.cpp:604` -- `dump_mem_info_lightweight`, the lightweight locate path being modeled
- `src/find_mems.cpp` (classify additions) -- `MEMRunClassification`, `classify_mem_runs`, `enumerate_bwt_runs`, `enumerate_tag_runs`, `write_classification_row`, `CLASSIFY_TSV_HEADER`
- `include/pangenome_index/algorithm.hpp:655-738` -- `find_mems_function`, the three-phase MEM finder (backward-extend seed, forward-extend right-maximal, backward-extend next-x)
- `include/pangenome_index/light_tag_index.hpp` -- LightTagIndex API (`run_id_at`, `run_start_bwt`, `num_runs`, `bwt_size`)
- `src/r-index.cpp:751` -- `backward_extend_encoded`
- `src/r-index.cpp:787` -- `forward_extend_encoded` (structural obstacle detail above)
- `src/r-index.cpp:1212` -- `run_id_and_offset_at`
- `src/r-index.cpp:1483` -- `locateNext`
- `xy-test/verify_partition.awk` -- aggregate script that produced the tables above

## Commits on this branch

- `b3c2a57` -- initial classify hook
- `605c33c` -- 4-bucket priority-ordered tag partition
- `9221d4f` -- reformulate bucket 1 as "contains any bwt_start"
- `40c7d5c` -- add `avg_bwt_run_length` and `avg_tags_per_bwt_run` (v4 schema)
