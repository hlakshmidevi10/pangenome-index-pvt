# Performance Report Template

A standardized, dataset-centric format for reporting `find_mems` + `gafpack`
performance and MEM-coverage characteristics for a given (graph, reads,
config) combination — optionally comparing multiple `find_mems` retrieval
strategies (flipped / legacy / legacy+`--tag-head-samples=sN`) on the same
dataset.

**This is not `RESEARCH_JOURNAL.md`.** The journal records hypotheses,
experiments, and decisions with reasoning (`## JR-NNN` entries, append-only).
This template records a data snapshot — a structured, comparable report with
no narrative required. A journal entry can *embed* a report built from this
template (e.g. as a "dataset characterization" `## JR-NNN` entry) when the
snapshot is itself worth a permanent record; use the template standalone for
one-off or exploratory analysis that doesn't need journal permanence.

Established 2026-08-21, seeded from the HPRC v1 chr6 (PGGB), MEM_LEN=25
flipped-vs-legacy+samples(s=8) analysis. See `RESEARCH_JOURNAL.md` JR-021 for
that worked example in full and the discussion that produced this format.

---

## Sections & field definitions

### Reference characteristics

| Field | Source | Notes |
|---|---|---|
| Graph | — | e.g. "HPRC v1 chr6 (PGGB)" |
| Nodes | `gbz_stats -g <gbz>` → `Nodes` | |
| Edges | `gbz_stats -g <gbz>` → `Edges` | |
| Graph sequence (deduplicated) | `gbz_stats -g <gbz>` → `Sequence` | bp of unique graph sequence |
| Total haplotype path length | `gbz_stats -p <gbz>` → `Path length` | bp, single orientation, all haplotypes summed |
| Haplotypes / samples / paths | `gbz_stats -i <gbz>` → `Metadata:` line | |

**Gotcha:** `gbz_stats -i` also prints a line labeled `Total length`, which is
the sum of GBWT **node-visits** (BWT length in alphabet symbols), **not** a
bp count — do not use it for sequence-content reporting. Use `-p Path
length` for "how much haplotype sequence does this graph carry" and `-g` for
node/edge/sequence counts. (Documented independently in
`mem-projection/pangenome-pipeline/CLAUDE.md` "Known footguns" — this
template just restates it because it's easy to grab the wrong flag.)

### Read characteristics

| Field | Notes |
|---|---|
| Reads | count × length, single/paired |
| Source | haplotype/sample the reads were pulled from; note in-graph vs out-of-graph |
| Error profile | clean vs noisy; simulator + rate (e.g. wgsim `-e`); only call it "Illumina-like" if you mean a uniform per-base substitution model — it is not a true position-dependent Illumina quality-decay profile |
| Read file | path + size |

### MEM Coverage Statistics

Identical across `find_mems` configs querying the same reads/index at the
same `MEM_LEN`/`MIN_OCC` — the MEM set doesn't depend on flipped vs legacy vs
samples, only the *retrieval method* does. **Verify this** (don't assume it)
via the correctness gate below before reporting the section once for
multiple configs.

| Field | Source | Formula / notes |
|---|---|---|
| Total Reads | find_mems log, "Total number of reads" | |
| Total MEMs | find_mems log, "Total number of MEMs outputted" | |
| Average MEM length | find_mems log, "Average MEM length" | |
| Total Occurrence (Σ mem.size) | find_mems log, "Total occurrence/bwt positions scanned" | |
| Avg Occurrence / MEM | find_mems log, "Average MEM occurrence (mem.size)" | |
| Avg Tag Runs / MEM (SA interval) | find_mems log, "Average tag runs per MEM" | one MEM = one SA interval |
| Average Tag Run Length | find_mems log, "Global n/r ratio (Σmem.size / Σtag_runs)" | **Use Global, not "Mean per-MEM n/r ratio."** The per-MEM mean is outlier-skewed by a long tail of high-occurrence MEMs — e.g. the JR-021 dataset: mean=78.0 vs global=21.6, a 3.6x distortion. This exact footgun is called out in the pipeline `CLAUDE.md`. |
| Total Entries Written (find_mems, pre-dedup) | `SUMMARY.tsv` `bin_bytes / 16` (16-byte v2 record) | **Do not trust find_mems.log's abbreviated "K/M" text** — it appears to truncate rather than round (an exact value of 1,482,232 printed as "1.4M"). Always compute the exact integer from `bin_bytes / 16` in the harness `SUMMARY.tsv`, or the raw `.bin` file size. |
| Total Occurrence on graph after MEM projection (post gafpack dedup) | `SUMMARY.tsv` `gafpack_total_entries`, or the gafpack log | GAF-equivalent entry count after `--dedup-read-node` |

### Performance

Report as a comparison table (config as columns, `Δ` as the last column) when
multiple `find_mems` configs are being compared on the same dataset; a single
column is fine for a one-config report.

| Field | Source |
|---|---|
| find_mems time | `SUMMARY.tsv` `wall_s` for the `find_mems` step, mean ± stdev over N≥3 warm-cache trials |
| find_mems memory (peak RSS) | `SUMMARY.tsv` `maxrss_mb`, or `summarize.py`'s "peak (reported)" |
| gafpack time / memory | same, `gafpack` step |
| Total wall | find_mems wall + gafpack wall (`summarize.py`'s "STEPS 09+10 COMBINED WALL") |
| MEM throughput (MEMs/s) | Total MEMs ÷ Total wall (combined find_mems+gafpack wall — **not** find_mems-only wall) |

**Always use N≥3 warm-cache trials via `perf_harness.sh`, not single-trial
`query.sh` runs**, when reporting timing you intend to stand behind —
single-trial numbers can be off by double-digit percentages on an outlier
trial (see JR-008's legacy trial-1: 46.8s vs a 37.9-39.3s cluster on the
other trials). 2 untimed warmups + N timed trials, contiguous per-config
ordering, is the harness default; don't change it without a reason and don't
skip the warmups.

### Disk storage

| Field | Source |
|---|---|
| r-index (`.ri`) | file size (`ls -la` or `PROVENANCE.txt`), report in MiB/GiB (binary) to match existing journal conventions |
| tag-index, lightweight (`.ltags`) | file size |
| tag-samples (`.tag_samples.s<N>`) | file size, only relevant when comparing a `--tag-head-samples` config |

---

## Correctness gate — run this before trusting any Performance comparison

Confirm the configs being compared are computing the *same thing* before
presenting their numbers side by side:

- Coverage MD5 (`alignment_coverage.csv`) identical across all trials, all configs
- `.bin` record count / `bin_bytes` identical
- `gafpack_total_entries` identical
- `Total MEMs` identical

If any of these differ, the configs are not apples-to-apples comparable and
the Performance table must not be presented as a clean A/B comparison — say
so explicitly instead. Coverage-CSV byte-identity is the strongest available
signal (see `RESEARCH_JOURNAL.md` JR-007) — stronger than `validate_gaf`'s
sampled check, which the perf harness doesn't even run by default (it's
`--gaf`-gated and coverage-only by default per `BUILD_QUERY_LAYOUT.md`).

---

## Worked example

See `RESEARCH_JOURNAL.md` JR-021 for the full HPRC v1 chr6 (PGGB), MEM_LEN=25,
alt-noisy-reads report (flipped vs legacy+samples s=8) that this template was
extracted from.
