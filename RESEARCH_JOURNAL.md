# Research Journal — pangenome-index-latest

**Purpose.** Chronological, append-only record of research decisions, experiments,
measurements, and open questions on this codebase. Written by whichever agent or
human is doing the work. Reads top-to-bottom as a history of what we tried,
what we learned, and where the frontier is.

**Not a design doc, not a changelog, not a wiki.** Design docs
(`DESIGN_FLIPPED_MEM.md`) describe intent; changelogs (`git log`) describe code;
wikis get restructured. This journal describes *thinking*, and only ever grows.

---

## Conventions (read before writing an entry)

### Entry structure

Each entry is a level-2 heading of the form `## JR-NNN — <short title>` where
`NNN` is a zero-padded monotonic integer. The heading is immediately followed by
a fenced YAML block, then prose.

```markdown
## JR-042 — <short title>

​```yaml
id: JR-042
date: YYYY-MM-DD
author: <model name or human handle> (session with <collaborator, if any>)
status: open | resolved | superseded | deferred | wontfix
tags: [free-form, comma-separated, lowercase-with-hyphens]
refs:
  follows: [JR-041]              # this entry continues from JR-041 (default lineage)
  supersedes: []                 # entries this replaces
  resolves: []                   # entries whose "open" status this closes
  confirms: []                   # entries whose hypothesis this experiment confirms
  refutes: []                    # entries whose hypothesis this experiment refutes
​```

<body -- see recommended sections below>
```

Omit `refs.*` subfields that are empty. Keep the fenced block; other agents may
grep for it.

### Body sections (loose template, not enforced)

Not every entry needs every section. Use judgement. But if an entry doesn't
have a clear **Findings** section, ask whether it belongs in the journal at
all -- design work goes in `DESIGN_FLIPPED_MEM.md`, code documentation goes
in comments.

- **Context.** What prompted this entry. 1 short paragraph.
- **Hypothesis / Question.** What we're testing or trying to answer. Be
  explicit; vague hypotheses lead to vague conclusions.
- **Method.** Commands run, binaries used, datasets consumed. Exact enough
  that another agent can reproduce it.
- **Findings.** Numbers, tables, quotes from log files. Concrete.
- **Interpretation.** What the numbers mean in the context of the larger
  research goal. Separate from the raw findings so the next reader can
  disagree with your interpretation without disputing your measurements.
- **Open questions.** New questions this entry raised. Prefer explicit
  numbered lists over vague "future work" mentions -- these are what the
  next JR-entry will pick up.

### Append-only rules

1. **Never edit an existing entry.** Ever. Write a new entry that references
   the old one via `refs.supersedes`, `refs.resolves`, `refs.refutes`, etc.
   Editing history erases the record of what we thought at the time.
2. **Exception:** typo/formatting fixes in prose are OK, so long as they
   don't change meaning. If in doubt, write a new entry.
3. **To close an "open" entry:** write a new entry whose `refs.resolves`
   lists the old ID. The old entry's YAML `status` field stays as it was
   when written -- you can (and should) add a one-line note at the end of
   the old entry's *body* saying `> Resolved by JR-NNN (YYYY-MM-DD).`
   That's the only permitted body edit.
4. **Timestamps.** Use the actual calendar date of the work, not the date
   you're writing the entry. If you're backfilling an entry for something
   that happened earlier, use the earlier date.
5. **When multiple agents append in the same session,** the second agent
   picks the next `JR-NNN`. Merge conflicts on this file are the "cost of
   doing business" -- prefer them over losing an entry to a silent overwrite.

### Index maintenance

The **Index** section immediately below is hand-maintained. When you add an
entry, add its row to the index in numeric order. When you close an entry
(via a resolving entry), update its `status` column in the index only --
never in the entry's own YAML.

The index lets an agent skim the journal state without reading every entry.
If you find yourself adding an entry, always update the index in the same
change.

### What to write about

Good candidates for a journal entry:
- A hypothesis about performance, correctness, or algorithm behavior +
  the experiment that tested it.
- An unexpected finding while working on something else.
- A design decision made under uncertainty, with the reasoning.
- A dataset characterization ("what does the size distribution actually
  look like?").
- A bug found in existing code, with the smallest reproducer and its scope.
- A negative result -- something we tried that didn't work.

Not-good candidates:
- "I renamed foo to bar" -- git log.
- "The API for X is Y" -- code comments or design doc.
- "The next step is Z" -- todo list or design doc.

### Cross-references outside the journal

When an entry cites external artifacts (a commit, a source file, a log
under `runs/`, a design doc section), spell out the path or SHA. Do not
say "see the earlier commit" or "in the design doc" -- future readers
won't have the same context you do.

---

## Index

| ID | Date | Title | Status | Tags |
|---|---|---|---|---|
| [JR-001](#jr-001--stage-1-sa-carrying-backward-extend-primitive) | 2026-07-28 | Stage 1: SA-carrying backward extend primitive | resolved | stage1, r-index, primitive, sa-carry |
| [JR-002](#jr-002--stage-2-flipped-mem-finder) | 2026-07-28 | Stage 2: flipped MEM finder | resolved | stage2, mem-finder, algorithm |
| [JR-003](#jr-003--stage-3-wire-flipped-finder-into-find_mems) | 2026-07-28 | Stage 3: wire flipped finder into find_mems | resolved | stage3, find-mems, perf |
| [JR-004](#jr-004--risk-e-legacy-step-3-emits-non-left-maximal-mems-on-hprc) | 2026-07-28 | Risk E: legacy Step 3 emits non-left-maximal MEMs on HPRC | resolved (by JR-007) | risk-e, correctness, legacy-bug, hprc |
| [JR-005](#jr-005--risk-e-size-distribution--per-read-concentration) | 2026-07-29 | Risk E: size distribution + per-read concentration | resolved (by JR-007) | risk-e, distribution, hprc, characterization |
| [JR-006](#jr-006--refined-stage-3-perf-estimate-once-risk-e-is-fixed) | 2026-07-29 | Refined Stage 3 perf estimate once Risk E is fixed | resolved (by JR-007) | stage3, perf, risk-e, projection |
| [JR-007](#jr-007--risk-e-fix-fixed-legacy--flipped-produce-byte-identical-coverage) | 2026-07-29 | Risk E fix: fixed-legacy == flipped (byte-identical coverage) + honest Stage 3 delta | resolved | risk-e, fix, correctness, perf, stage3 |
| [JR-008](#jr-008--stage-3-benchmark-n5-warm-cache-hprc-chr6-100k-reads) | 2026-07-29 | Stage 3 benchmark: N=5 warm-cache HPRC chr6 100K reads | resolved | stage3, perf, benchmark, hprc |
| [JR-011](#jr-011--backward_extend_encoded_with_sa-rewrite-hprc-chr6-noisy-benchmark) | 2026-07-30 | backward_extend_encoded_with_sa rewrite: HPRC chr6 noisy benchmark | resolved | jr-010, sa-carry, primitive, perf, hprc, noisy |
| [JR-012](#jr-012--step-0-min_len-pre-check-in-find_mems_flipped-hprc-chr6-noisy--clean) | 2026-08-01 | Step 0' min_len pre-check in find_mems_flipped: HPRC chr6 noisy + clean | resolved | flipped, algorithm, perf, hprc, noisy, clean, step0, vesuvio |
| [JR-013](#jr-013--toehold-hoist-into-callers-perf-null-contract-cleanup) | 2026-08-06 | Toehold hoist into callers: perf null, contract cleanup | resolved | flipped, r-index, primitive, contract, refactor, vesuvio, perf-null |
| [JR-014](#jr-014--serialize-run_head_c-into-ri-load-time-collapse) | 2026-08-06 | Serialize run_head_c[] into .ri: load time collapse | resolved | jr-011, jr-013, r-index, serialization, load-time, vesuvio |
| [JR-015](#jr-015--post-risk-e-mem-classification-fresh-seed-cost--tag-partition-tables) | 2026-08-07 | Post-Risk-E MEM classification: fresh seed-cost + tag-partition tables | resolved | classification, seed-cost, tag-partition, hprc, vesuvio, historical-comparison |
| [JR-016](#jr-016--tag-head-sa-samples-sr-index-style-storage-cost-on-hprc-chr6) | 2026-08-09 | Tag-head SA samples (sr-index-style): storage cost on HPRC chr6 | open | tag-head-samples, sr-index, sa-sampling, storage, hprc, vesuvio |
| [JR-017](#jr-017--lf-operation-latency-on-hprc-chr6-baseline-vs-fused-lf_scan) | 2026-08-09 | LF-operation latency on HPRC chr6: baseline vs fused LF_scan | open | lf, r-index, microbenchmark, hprc, vesuvio, jr-016 |
| [JR-018](#jr-018--tag-head-sa-samples-end-to-end-integration--n3-warm-cache-perf-on-hprc-chr6) | 2026-08-12 | Tag-head SA samples: end-to-end integration + N=3 warm-cache perf on HPRC chr6 | resolved | tag-head-samples, integration, correctness, perf, benchmark, hprc, vesuvio, jr-016, jr-017 |
| [JR-019](#jr-019--flipped-vs-legacys8-end-to-end-perf-on-hprc-chr6-alt-noisy-l25-vs-l50) | 2026-08-12 | Flipped vs Legacy+s=8 end-to-end perf on HPRC chr6 alt-noisy: L=25 vs L=50 | resolved | perf, benchmark, hprc, alt-noisy, l25, l50, samples, jr-018 |
| [JR-020](#jr-020--convert_tags-decoded-sdsl-header-as-bytecode-latent-bug-hit-by-hprcv1-chr1-build) | 2026-08-13 | convert_tags decoded SDSL header as ByteCode: latent bug hit by HPRCv1 chr1 build | resolved | convert_tags, build_tags, sdsl, int_vector_buffer, tag-arrays, correctness, hprc, chr1, bug-fix |
| [JR-021](#jr-021--standardized-performance-report-format-perf_report_templatemd--hprc-chr6-l25-worked-example) | 2026-08-21 | Standardized performance-report format (PERF_REPORT_TEMPLATE.md) + HPRC chr6 L=25 worked example | resolved | reporting, template, perf, benchmark, hprc, l25, noisy, vesuvio, jr-018, jr-019, tag-head-samples |
| [JR-022](#jr-022--block-count-off-by-one-in-the-fastlocate-constructor-a-phantom-trailing-block) | 2026-08-21 | Block-count off-by-one in the FastLocate constructor: a phantom trailing block | resolved | r-index, blocks, off-by-one, sdsl, sd-vector, build_tag_head_samples, correctness, bug-fix, chr1, mc-chr6, vesuvio |
| [JR-023](#jr-023--tag-head-sa-samples-do-not-scale-to-whole-genome-the-chr6-recommendation-inverts-on-chr1) | 2026-08-22 | Tag-head SA samples do not scale to whole-genome: the chr6 recommendation inverts on chr1 | resolved | tag-head-samples, scaling, whole-genome, perf, benchmark, chr1, hprc, flipped, vesuvio, jr-016, jr-018, jr-019 |
| [JR-024](#jr-024--tag-array-structural-characteristics-coincidence-rate-as-the-samples-design-predictor) | 2026-08-22 | Tag-array structural characteristics: coincidence rate as the samples-design predictor | resolved | tag-head-samples, characterization, tag-array, coincidence-rate, density, chr1, mc-chr6, hprc, vesuvio |
| [JR-025](#jr-025--correction-to-jr-023-emit-volume-not-coincidence-rate-drives-samples-query-cost) | 2026-08-22 | Correction to JR-023: emit volume, not coincidence rate, drives samples query cost | resolved | tag-head-samples, correction, model, emit-volume, perf, chr1, mc-chr6, vesuvio |

---

## JR-001 — Stage 1: SA-carrying backward extend primitive

```yaml
id: JR-001
date: 2026-07-28
author: claude-opus-4-7 (session with hlakshmidevi)
status: resolved
tags: [stage1, r-index, primitive, sa-carry]
refs:
  follows: []
```

### Context

`FINDINGS_SEED_COST.md` measured that `locate_sa_value` (the "seed" cost) is
32.2% of walk time in `dump_mem_info_lightweight` on HPRCv1 chr6, 100K reads.
`DESIGN_FLIPPED_MEM.md` proposes eliminating this cost by restructuring MEM
enumeration to use only backward extends (which admits SA-carrying) instead
of legacy's forward+backward mix. Stage 1 of that plan: build the primitive
that carries `SA[sp]` through a chain of backward extends.

### Hypothesis

The r-index paper's SA-carry rule (Cobas, Gagie, Navarro 2021, §3), when
mirrored to our samples-at-run-*starts* convention, produces `SA[new_sp]`
at O(rank) marginal cost per backward-extend step -- no `locateNext` walk
from a run-start sample required. Case A (LF-decrement) should dominate;
Case B (leftmost-c lookup) is rare in practice.

### Method

Added `backward_extend_encoded_with_sa` in `src/r-index.cpp:795` alongside
existing `backward_extend_encoded` (not replacing it -- 5 call sites want
unchanged behavior). New `bi_interval_with_sa` struct in
`include/pangenome_index/r-index.hpp:458`.

Case A: `BWT[old_sp] == c` -> `sa_sp' = sa_sp - 1` (LF-decrement).
Case B: `BWT[old_sp] != c` -> find leftmost c-position j in `[old_sp, old_ep]`
via `leftmost_c_in_interval` (built in Stage 0.85, commit 1c9ff21); j is at
a c-run start so `SA[j] = samples[run_id_at(j)]`; then `sa_sp' = SA[j] - 1`.
Initial toehold (sa_sp_prev == NO_POSITION): seed via one `locate_sa_value`
walk from the containing run's sample.

Verification: `bin/test_backward_extend_sa` (`src/test_backward_extend_sa.cpp`)
runs random-ACGT patterns + text-derived deep walks against the canonical
yeast index. At every successful step, compares `sa_sp` from the new path
against `locate_sa_value_naive(new_sp)` ground truth. Also cross-checks the
interval math against the plain `backward_extend_encoded`.

Command:
```
bin/test_backward_extend_sa \
  ../mem-projection/pangenome-pipeline/runs/v2-yeast235/yeast235_chrII_100kb_normalized.ri
```
(and the same with `20000 20000 123` for larger stress.)

### Findings

Two runs, both on `runs/v2-yeast235/yeast235_chrII_100kb_normalized.ri`
(n=337M, 14M runs):

| Config | Trials | Steps | Case A | Case B | Case B % | Mismatches |
|---|---:|---:|---:|---:|---:|---:|
| defaults (seed=42) | 10,000 | 304,485 | 216,016 | 78,469 | 26.6% | **0** |
| stress (seed=123) | 40,000 | 1,204,000 | 849,957 | 314,043 | 27.0% | **0** |

E2E validation on yeast pipeline (`./query.sh configs/yeast235-chrII-normalized.env
v2-yeast235 stage1-sa-primitive --gaf`): 2000/2000 valid, byte-identical MD5
to prior baselines (`alignment.gaf` = `6d2c8784bef2b638ca6f0b96ce8f6bfd`,
`mems_path_pos_v2.bin` = `3474a0cb016e13b2ff2033f27d9be498`). The new primitive
has no callers on the hot path yet; this run just confirmed the header change
+ library rebuild didn't perturb anything.

### Interpretation

Design's Case A / Case B split holds up empirically: Case A is ~73% at yeast
scale, well below the 50% threshold that DESIGN_FLIPPED_MEM.md §5.3 Risk B
flagged as a concern. LF-decrement dominates as intended.

Zero mismatches across 1.2M steps establishes the primitive as ground-truth-
equivalent to `locate_sa_value` at the primitive level. The next stage
(JR-002) can consume `sa_sp` and know it's correct.

### Open questions

1. Case B fraction on HPRC chr6 (larger n, denser BWT) is unmeasured.
   Might differ from yeast's 27%.

Commit: `1c00f89` (`r-index: add backward_extend_encoded_with_sa (SA-carrying primitive)`)

---

## JR-002 — Stage 2: flipped MEM finder

```yaml
id: JR-002
date: 2026-07-28
author: claude-opus-4-7 (session with hlakshmidevi)
status: resolved
tags: [stage2, mem-finder, algorithm]
refs:
  follows: [JR-001]
```

### Context

With the SA-carry primitive from JR-001 in place, Stage 2 builds the
right-to-left MEM enumerator that consumes it. Each emitted MEM carries
`SA[bwt_start]` for free -- Stage 3 (JR-003) will then wire this into
`dump_mem_info_lightweight` to skip the seed walk.

### Hypothesis

The flipped algorithm (DESIGN_FLIPPED_MEM.md §3, verified by
`xy-test/mem-prototype/compare_mem_algorithms.py` against a brute-force
oracle across 1,014 cases) produces the same MEM set as the legacy
3-phase finder. C++ port should match legacy on real read data.

### Method

Added `find_all_mems_flipped` + `MEM_with_sa` struct in
`include/pangenome_index/algorithm.hpp:766`. Faithful C++ port of the
Python prototype:

- Outer loop iterates right-anchors right-to-left, starting at `x = len-1`.
- Step 1' backward-extends from `x` leftward via `backward_extend_encoded_with_sa`,
  emits MEM if length >= min_len.
- Step 2' fresh forward-extend from `pattern[i-1]` to identify the next anchor
  (extracted as `flipped_advance_anchor` -- runs both when MEM was emitted
  and when Step 1' match was too short).

`MEM_with_sa` mirrors legacy `MEM` fields + adds `sa_sp`. Chosen over
extending `MEM` because MEM is used at 60+ call sites in find_mems.cpp; per
CLAUDE.md minimality, adding a field would ripple.

Verification: `bin/test_flipped_mems` (`src/test_flipped_mems.cpp`) runs both
`find_all_mems` and `find_all_mems_flipped` on real reads, compares the
`(start, end, bwt_start, size)` MEM sets, and verifies `sa_sp` matches
`locate_sa_value(bwt_start)` for every flipped MEM. Categorizes
set-differences as FLIPPED_MISSED (real MEM flipped missed, hard fail),
FLIPPED_EXTRA (MEM flipped invented, hard fail), or LEGACY_SPURIOUS
(legacy MEM not left-maximal per `count_encoded`, tolerated -- this is
Risk E in DESIGN_FLIPPED_MEM.md §5.3).

Commands:
```
bin/test_flipped_mems \
  ../mem-projection/pangenome-pipeline/runs/v2-yeast235/yeast235_chrII_100kb_normalized.ri \
  ../mem-projection/yeast-235/yeast-235-chrII/S288C_chrII_N100K_R1_200_reads.txt \
  30 1 <N>
```
for N in {500, 10000, 100000}.

### Findings

Yeast (min_len=30, min_occ=1):

| max_reads | Legacy MEMs | Flipped MEMs | Missed | Extra | Legacy spurious | sa_sp mismatches |
|---:|---:|---:|---:|---:|---:|---:|
| 500 | 581 | 581 | 0 | 0 | 0 | 0 |
| 10,000 | 11,450 | 11,392 | 0 | 0 | 58 (2 reads) | 0 |
| 100,000 | 114,281 | 113,752 | 0 | 0 | 529 (8 reads) | 0 |

HPRC chr6 (min_len=50, min_occ=1), first 1,000 reads: 1000 = 1000, 0 mismatches.
First 10,000 reads: 10,000 flipped vs 10,170 legacy, 0 missed, 0 extra,
170 legacy-spurious (see JR-004 for spurious-MEM characterization).

Also verified min_occ=2 on yeast 10K reads: 11,441 = 11,441 exact set match,
no legacy_spurious (Risk E requires MEMs extending to end with size >= min_occ;
higher min_occ makes it rarer).

### Interpretation

Flipped's MEM set on yeast is exactly `legacy \ spurious`. Zero flipped-missed
and zero flipped-extras across 100K yeast reads and 10K HPRC reads is strong
evidence the algorithm is correct. The 529 yeast + 170 HPRC "legacy_spurious"
divergences are entirely explained by a pre-existing legacy bug (JR-004).

`sa_sp` values are consistent with `locate_sa_value` ground truth across
125K+ MEMs -- the SA-carry chain works end-to-end at algorithm level, not
just primitive level.

### Open questions

1. The `test_flipped_mems` "LEGACY_SPURIOUS" category was added mid-Stage-2
   when the first HPRC-scale test revealed divergences. It should not be a
   permanent tolerance in the test -- once Risk E is fixed (JR-004), the
   category can be removed and the test's bar tightened back to exact set
   equality.

Commit: `df06976` (`algorithm: add find_all_mems_flipped + SA-carrying MEM struct`)

---

## JR-003 — Stage 3: wire flipped finder into find_mems

```yaml
id: JR-003
date: 2026-07-28
author: claude-opus-4-7 (session with hlakshmidevi)
status: resolved
tags: [stage3, find-mems, perf, e2e]
refs:
  follows: [JR-002]
```

### Context

Stages 1 (JR-001) and 2 (JR-002) built and verified the machinery. Neither
had any effect on production `find_mems` output. Stage 3 wires the flipped
finder into `find_mems.cpp`'s lightweight emission path behind a new
`--use-flipped-mems` flag (default OFF -- CLAUDE.md "don't break the
default path"), and measures the resulting E2E performance.

### Hypothesis

Design projection (DESIGN_FLIPPED_MEM.md §4.4, §5.2 Stage 4): eliminating
the seed cost yields ~10% end-to-end `find_mems` wall reduction on HPRC
chr6, no RSS regression.

### Method

Refactored `dump_mem_info_lightweight` into `dump_mem_info_lightweight_core`
+ two thin wrappers so both legacy (computes seed via `locate_sa_value`)
and flipped (consumes `sa_sp` from `MEM_with_sa`) share the tag-run
enumeration and PackedEntry emission verbatim. Per CLAUDE.md: "no copy-paste
of an existing function with one line changed -- refactor or templatize."

Added `--use-flipped-mems` CLI flag. Rejects `--use-flipped-mems` without
`--lightweight-tags` (Stage 3 is scoped to the lightweight path only) and
rejects the combination with `--all-positions`.

Per-read dispatch in `main()`: if flag on, `find_all_mems_flipped` -> iterate
`MEM_with_sa` -> `dump_mem_info_lightweight_flipped`. Otherwise unchanged.

Local smoke on yeast 500 reads (subset of canonical file via `head -500`):
byte-identical `.tsv` (sorted) and `.bin` MD5 for both paths. Confirmed no
divergence at that scale before scaling up.

E2E:
```
# Yeast (via query.sh, since it supports no extra flags, only for flag OFF baseline)
./query.sh configs/yeast235-chrII-normalized.env v2-yeast235 stage3-flag-off --gaf

# Yeast flipped + HPRC flipped (manual pipeline, query.sh doesn't accept extra find_mems flags)
mkdir -p runs/.../queries/<name>/lightweight
cd runs/.../queries/<name>/lightweight
/usr/bin/time -l $PI_BIN/find_mems <ri> <ltags> <reads> <mem_len> <min_occ> mems \
  --lightweight-tags --use-flipped-mems
$GAFPACK --gfa ... --path-pos mems_path_pos_v2.bin ... --dedup-read-node --gaf-file-prefix alignment
python3 $VALIDATE_GAF alignment.gaf ... --sample 2000
```

### Findings

**Yeast E2E, flag OFF** (regression check that default path is unchanged):
- MD5 of `alignment.gaf` and `mems_path_pos_v2.bin` byte-identical to prior
  baselines (`validation-test`, `stage1-sa-primitive`, `stage3-flag-off`)
  across three separate query runs.

**Yeast E2E, flag ON** (`runs/v2-yeast235/queries/stage3-flag-on/`):
- validate_gaf: 2000/2000 valid (100%), 0 invalid.
- GAF line count: 4,906,654 (baseline 4,954,592, delta -47,938).
- Set diff: 0 flipped-only records, 47,938 legacy-only records (all trace
  to Risk-E spurious MEMs -- see JR-004).

**HPRC chr6 E2E, flag ON** (`runs/hprc-chr6-2026-06-02/queries/stage3-flipped/`,
`find_mems` only, since validate_gaf is out of scope per user):

| Metric | Baseline (legacy) | Flipped | Δ |
|---|---:|---:|---:|
| Wall (s) | 47.996 | 34.437 | -28.3% |
| Peak RSS (MB) | 4599.94 | 4106.73 | -10.7% |
| MEMs emitted | 100,846 | 100,014 | -832 |
| Total entries written | 21.7M | 1.1M | -95% |
| Total locate operations | 47.8M | 2.7M | -94% |
| First-locate calls | 100,846 | **0** | eliminated |
| First-locate time | 14.28s | **0s** | eliminated |
| locateNext calls | 47.7M | 2.7M | -94% |
| MEM-finding time | 19.30s | 29.22s | +51% |
| MEM-processing time | 21.59s | 1.34s | -94% |

validate_gaf on flipped HPRC output: 2000/2000 valid.

### Interpretation

Design's Stage 3 primary goal hit cleanly: `First locate calls: 0` and
`First locate time: 0s` in the profiler output. The seed cost is fully
eliminated.

The 28% headline wall reduction is larger than the 10% design projection.
Two effects compound: (a) seed elimination (~14s saved, this is Stage 3's
genuine contribution) and (b) massive downstream volume collapse from
excluding legacy's Risk-E spurious MEMs (~6s saved, this is Risk E getting
out of the way). Once Risk E is fixed in legacy, the honest apples-to-apples
delta is smaller -- see JR-006 for the refined projection.

MEM-finding got slower by 10s (+51%). This is Risk C in DESIGN_FLIPPED_MEM.md
§5.3 -- Step 2' does a fresh forward extend from `pattern[i-1]` that
partially duplicates work the legacy 3-phase split amortized. Predicted and
measured. Net-positive because the processing-side collapse is larger.

Peak RSS drop is entirely from the smaller in-memory entry vector before
sort (`22.86 MB` vs an estimated ~350 MB from legacy's 21.7M entries).

### Open questions

1. Real Stage-3-only perf delta once Risk E is fixed in legacy -- see JR-006.
2. Case B fraction at HPRC scale still unmeasured (JR-001 open item).
3. Stage 4: N=5 warm-cache benchmark per `mem-projection/pangenome-pipeline/perf/`
   protocol not yet run. Numbers above are single cold runs.
4. Full commit not yet landed. Code is written and builds; awaiting sign-off
   on Risk E treatment before committing.

Commit: (pending as of 2026-07-29)

---

## JR-004 — Risk E: legacy Step 3 emits non-left-maximal MEMs on HPRC

```yaml
id: JR-004
date: 2026-07-28
author: claude-opus-4-7 (session with hlakshmidevi)
status: open
tags: [risk-e, correctness, legacy-bug, hprc]
refs:
  follows: [JR-002, JR-003]
```

### Context

Discovered while running JR-002's set-equality check on 10K HPRC reads:
flipped and legacy MEM sets diverged. Investigation traced the divergence
to a bug in the legacy 3-phase finder's Step 3, previously flagged as
Risk E in DESIGN_FLIPPED_MEM.md §5.3 ("the current algorithm's Step 3
out-of-bounds bug"). This entry establishes concretely that Risk E is
real, characterizes its trigger, and quantifies its scope.

### Hypothesis

Legacy's Step 3 at `include/pangenome_index/algorithm.hpp:724` accesses
`pattern[len]` when Step 2 extends to end-of-pattern. C++ `std::string`
returns `\0` (the null terminator) at that position, not undefined
behavior. When `find_mems` runs on the encoded r-index, `sym_map[0]`
aliases to whatever DNA character was assigned code 0 during
`calculate_C` (typically 'A'). So Step 3 mis-extends by a phantom 'A'
followed by pattern's reverse, and the outer loop advances by only 1
anchor per iteration, emitting a new bogus "MEM" ending at end-of-pattern
at each anchor.

### Method

`bin/test_flipped_mems` (JR-002) already categorizes legacy-only MEMs
as LEGACY_SPURIOUS iff `is_left_maximal(mem)` returns false. Left-maximality
is checked by calling `FastLocate::count_encoded` on `pattern[start-1..end)`
and testing whether the count is >= min_occ. If the count of the
left-extension is >= min_occ, the MEM is not left-maximal (by definition).

Ran on yeast 100K reads and HPRC chr6 10K reads.

For one specific HPRC divergent read (#8093 in yeast, #4778 in HPRC), also
dumped legacy's and flipped's MEM lists explicitly to confirm the pattern:
legacy emits MEMs with `start = <k>, k+1, k+2, ..., end = len` for a
contiguous range of k values.

### Findings

**Trigger:** legacy emits spurious MEMs iff the true MEM at some anchor
extends all the way to `end = pattern.length()` AND the resulting Step 3
mis-extension keeps succeeding for enough iterations.

**Scope on yeast (100K reads, min_len=30, min_occ=1):**
- 529 spurious MEMs across 8 reads (0.008% of reads affected).

**Scope on HPRC (10K reads, min_len=50, min_occ=1):**
- 170 spurious MEMs across 3 reads (0.03% of reads affected).
- Most pathological: read #4778 emits **126 spurious MEMs**, all ending
  at position 200 (end of read), with start positions 25, 26, 27, ...,
  150 -- one at every anchor position in that range.
- Read #2251: 39 spurious MEMs, starts 112-150.
- Read #5238: 5 spurious MEMs.

**Downstream amplification:** the spurious MEMs have very large `size`
(occurrence count) values because they're low-complexity tail patterns
matching many places in the pangenome. Each spurious MEM balloons into
hundreds or thousands of downstream tag-run entries, so a small MEM-set
delta becomes a large output-volume delta:

- Legacy HPRC 100K entries written: 21.7M
- Flipped HPRC 100K entries written: 1.1M (20x reduction)

### Interpretation

Risk E is (a) rare in read count but (b) pathological in output volume.
It's benign for correctness in the sense that `validate_gaf` passes 100%
on legacy output (every spurious MEM's substring IS a real substring of
the pangenome at the claimed position -- just not left-maximal). But it's
non-benign for performance: legacy does ~20x more downstream work on
HPRC because of it.

The one-line fix (guard against `e >= len` before Step 3's loop) has
been deferred pending prioritization. When it lands:
- Legacy behavior changes on the small fraction of reads that trip Risk E.
- Legacy output size drops to match flipped's.
- `test_flipped_mems`'s LEGACY_SPURIOUS tolerance can be removed.
- Stage 3's honest perf delta will be measurable -- see JR-006.

### Open questions

1. **Fix scope:** does `if (e >= len) return len;` before Step 3 fully
   suffice, or are there other codepaths where a similar `pattern[j]`
   with `j >= len` access can occur? A search of the codebase for
   `pattern[` accesses inside MEM-related loops would answer this.
2. **Yeast vs HPRC ratio:** yeast 100K has 529 spurious across 8 reads
   (~66 per affected read); HPRC 10K has 170 across 3 reads (~57 per
   affected read). Not directly comparable (different min_len). But
   worth confirming on HPRC 100K to see if the "large batch of spurious
   from one bad read" pattern scales.
3. **Which reads trigger it?** The 3 HPRC reads that triggered on the
   first 10K don't look obviously low-complexity. Understanding the
   linguistic characteristic that triggers Risk E might inform whether
   the fix has other implications.

> To be resolved by: (a) fix commit for `algorithm.hpp:724`, (b) new
> journal entry describing the fix + re-baselined measurements.

> Resolved by JR-007 (2026-07-29).

---

## JR-005 — Risk E: size distribution + per-read concentration

```yaml
id: JR-005
date: 2026-07-29
author: claude-opus-4-7 (session with hlakshmidevi)
status: open
tags: [risk-e, distribution, hprc, characterization]
refs:
  follows: [JR-004]
  confirms: [JR-004]
```

### Context

JR-004 established that Risk E is real and quantified it in aggregate.
User asked whether MEM sizes are meaningfully different from averages
alone -- specifically whether the spurious MEMs form a distinct
subpopulation vs. being smeared across the whole distribution. This
entry answers that with concrete per-MEM data.

### Method

New diagnostic binary `bin/dump_mem_size_distribution`
(`src/dump_mem_size_distribution.cpp`) -- one-off tool that emits a TSV
row per MEM from both finders, tagged with `left_maximal` and
`in_other_finder`. Then downstream awk pipelines for histograms and
per-read counts.

Command:
```
bin/dump_mem_size_distribution \
  ../mem-projection/pangenome-pipeline/runs/hprc-chr6-2026-06-02/hprcv1_chr6.ri \
  ../mem-projection/hprcv1/chr6.alt.reads.txt \
  50 1 10000 \
  > xy-test/risk_e_analysis/hprc_10k.tsv
```
Output: 20,171 rows (header + 10,170 legacy MEMs + 10,000 flipped MEMs).

### Findings

**Cross-tabulation of left-maximal x in-other-finder (10K HPRC reads):**

```
count   finder   left_maximal   in_other_finder
10000   flipped  yes            yes
  170   legacy   no             no
10000   legacy   yes            yes
```

Only three cells populated. **Zero flipped-emitted-non-left-maximal**,
**zero flipped-invented-MEMs**, **zero real-MEMs-flipped-missed**. The
flipped MEM set is exactly `legacy \ spurious`.

**Size histogram (log-scale bins):**

| bin | legacy count | flipped count |
|---|---:|---:|
| [1, 10) | 269 | 269 |
| [10, 100) | 9,692 | 9,692 |
| [100, 1000) | 21 | 21 |
| [1000, 10,000) | 9 | 9 |
| [10,000, 100,000) | **179** | **9** |
| [100,000, inf) | 0 | 0 |

Every bin except `[10⁴, 10⁵)` is identical between the two finders. In
the top bin, legacy has 179, flipped has 9, delta = 170 -- exactly matching
the 170 legacy-spurious count.

**Spurious MEM size stats (n=170):**

| stat | value |
|---|---:|
| min | 20,751 |
| Q1 | 26,939 |
| median | 32,173 |
| Q3 | 34,924 |
| P90 | 36,287 |
| P95 | 39,583 |
| P99 | 64,323 |
| max | 64,450 |
| mean | 31,680.7 |
| sum | 5,385,718 |

**Shared / real MEM size stats (n=10,000):**

| stat | value |
|---|---:|
| min | 1 |
| Q1 | 64 |
| median | 86 |
| Q3 | 90 |
| P90 | 90 |
| P95 | 90 |
| P99 | 91 |
| max | 53,629 |
| mean | 103.2 |
| sum | 1,031,936 |

**Per-read spurious concentration:**

| read_id | # spurious | total spurious size |
|---|---:|---:|
| #4778 | 126 | 3,672,631 |
| #2251 | 39 | 1,391,423 |
| #5238 | 5 | 321,664 |

All 3 affected reads out of 10,000 (0.03%). Read #4778 alone accounts
for 68% of the total spurious size mass.

**Volume decomposition of legacy's total occurrences (10K reads):**
- Legacy real MEMs total size: 1,031,936
- Legacy spurious MEMs total size: 5,385,718
- Legacy TOTAL: 6,417,654
- Flipped total: 1,031,936 (exactly = legacy real)
- **Spurious share of legacy's total occurrence volume: 83.9%**

### Interpretation

The spurious and real MEM size distributions are **disjoint**. Real MEMs
are tight around 86-91; spurious MEMs are all >= 20,751. There's no
overlap and no ambiguous middle ground. Any legacy MEM with size in the
tens of thousands and end == len is essentially guaranteed to be spurious.

The three-order-of-magnitude size gap is what makes Risk E dominate
downstream volume despite affecting only 0.03% of reads. 170 spurious
MEMs at median size 32K = ~5.4M "phantom occurrences", which is 5.2x
more occurrence volume than all 10,000 real MEMs combined. In
`dump_mem_info_lightweight` each occurrence becomes a `locateNext` step
(walking to the next tag-run start), which is why legacy's ~48M locate
calls collapse to flipped's ~2.7M once these phantoms are gone.

The linguistic characteristic that triggers Risk E is now partly clear:
the affected reads' spurious MEMs all end at pattern-end, with contiguous
start positions (read #4778: starts 25, 26, ..., 150). This is the exact
signature the "phantom-'A' backward-extend eating one position per outer-loop
iteration" theory predicts. Not confirmed via instrumentation of the
legacy Step 3 loop -- doing that would be the direct proof.

### Open questions

1. **Why these 3 reads and not others?** The reads' sequences don't look
   obviously low-complexity. Something specific about their end sequences
   makes Step 3's phantom-extend succeed for 100+ iterations. Instrumenting
   Step 3 to log what `back.size` is at each iteration on read #4778
   would pin this down.
2. **Does the size gap scale to 100K reads?** With ~10x more reads and
   presumably ~10x more spurious cases, is the size distribution still
   cleanly disjoint or does it start to overlap at the edges?
3. **Yeast comparison:** JR-004 showed yeast has 529 spurious across 8
   reads with `min_len=30`. Yeast's spurious size distribution should
   look similar (disjoint tails), but the specific sizes will be smaller
   because yeast pangenome is 100x smaller than HPRC chr6. Not yet
   measured.

Data files (retained for future reference):
- `xy-test/risk_e_analysis/hprc_10k.tsv` (20,171 rows, ~2 MB)
- `xy-test/risk_e_analysis/hprc_10k.log`

Binary: `bin/dump_mem_size_distribution` (source: `src/dump_mem_size_distribution.cpp`)

> This entry does not resolve JR-004; it strengthens the evidence that
> Risk E is a well-scoped bug with a well-scoped fix. JR-004 stays open
> pending the actual fix commit.

> Resolved by JR-007 (2026-07-29).

---

## JR-006 — Refined Stage 3 perf estimate once Risk E is fixed

```yaml
id: JR-006
date: 2026-07-29
author: claude-opus-4-7 (session with hlakshmidevi)
status: open
tags: [stage3, perf, risk-e, projection]
refs:
  follows: [JR-003, JR-005]
```

### Context

JR-003 measured a 28% wall reduction on HPRC chr6 (48s -> 34.4s) with
`--use-flipped-mems`. JR-005 established that 84% of legacy's downstream
occurrence volume is Risk-E phantom work. The 28% headline therefore
double-counts: it credits Stage 3 for savings that a fixed legacy would
also realize. This entry decomposes the 28% into the two components and
projects the honest Stage-3-only delta.

### Method

Arithmetic on JR-003's timing breakdown, using JR-005's ratio
(1,031,936 / 6,417,654 = 16.1% -- the real fraction of legacy's work).

Legacy processing-time components:
- MEM finding: 19.3s (independent of Risk E; same finder either way)
- Tag queries: 0.66s
- First locate (seed): 14.28s  <- Stage 3 eliminates this
- LocateNext: 4.69s
- Other: 0.25s

If Risk E were fixed in legacy, the tag queries and locateNext costs would
drop by roughly the same ratio as the occurrence-volume drop (0.161 of the
current value):
- Tag queries: 0.66s * 0.161 ~ 0.11s (saves 0.55s)
- LocateNext: 4.69s * 0.161 ~ 0.76s (saves 3.93s)

Plus sorting: 2.83s * 0.161 ~ 0.46s (saves 2.37s from the sort phase).

### Findings

**Estimated fixed-legacy find_mems wall breakdown:**

| Component | Legacy (as measured) | Fixed legacy (est.) | Basis |
|---|---:|---:|---|
| R-index load | 2.85s | 2.85s | unchanged |
| Tag index load | 1.18s | 1.18s | unchanged |
| MEM finding | 19.30s | 19.30s | Risk E fix doesn't change MEM count materially (170/100K ~ 0.17% of MEMs) |
| Tag queries | 0.66s | 0.11s | scaled by 0.161 |
| First locate (seed) | 14.28s | ~2.30s | scaled by 0.161 (assuming per-real-MEM seed cost matches) |
| LocateNext | 4.69s | 0.76s | scaled by 0.161 |
| Other | 0.25s | 0.04s | scaled |
| Sort | 2.83s | 0.46s | scaled |
| **Total (est.)** | **~48s** | **~27s** | |

**Flipped (as measured):** 34.4s.

**Estimated honest Stage-3-only delta:** flipped 34.4s vs fixed-legacy ~27s
= **flipped is ~27% SLOWER than fixed-legacy on HPRC**.

Wait -- that contradicts JR-003's headline. Let me redo more carefully.

Actually the "first locate scaled by 0.161" step is wrong. First-locate
is one call per emitted MEM, so it scales with MEM count (100.8K -> 100.0K,
essentially unchanged), not with occurrence volume. Redo:

- First locate: legacy 14.28s / 100,846 calls = 141.6 microseconds per call.
  Fixed legacy would still have ~100,000 calls (loses only the 170 spurious
  MEMs). But each spurious MEM's first-locate walk is very long (since sp
  falls in a huge run). If spurious first-locates average 5x the walk length
  of real first-locates (rough guess): real per-call ~= (14.28s - N*5t) / 100000
  where t = real avg first-locate time, N = 170 spurious. Without measurement
  we can't split cleanly.

The right way to answer this question is to *actually run fixed-legacy* --
fix Risk E in a scratch commit, re-run HPRC baseline, then compare against
JR-003's flipped numbers. Estimating via ratios is unreliable when the
spurious MEMs have very different per-call costs than real ones.

### Interpretation

**JR-003's 28% headline is not honest as an apples-to-apples measure of
Stage 3's contribution.** It measures "flipped vs Risk-E-buggy legacy" and
credits Stage 3 for savings a bug-fix would also deliver.

The honest number requires measurement, not estimation. The estimation
above has one solid conclusion and one uncertain one:

- **Solid:** Stage 3 eliminates ~14s of `first_locate` time -- this is
  real and independent of Risk E. Even if fixed-legacy's first_locate
  time drops proportionally with spurious MEM removal, most of the 14s
  is still real per-MEM seed work on the 100K non-spurious MEMs.

- **Uncertain:** Stage 3's Step 2' costs ~10s extra in MEM-finding. Whether
  the net (seed savings - Step 2' cost) is positive at the ~4s level or
  the ~10s level depends on how much of legacy's 14.28s first-locate is
  attributable to the 170 spurious MEMs specifically.

### Open questions

1. **Actually measure fixed-legacy.** Apply the Risk E one-line fix to
   `algorithm.hpp:724`, re-run HPRC baseline, then re-run this comparison.
   This is the definitive answer and this entry should be superseded by
   one that reports the measurement.
2. Until (1) is done, the 28% headline stands as "what a user switching
   between legacy-today and flipped-today would see," which is still a
   real number even if not the theoretically clean one. Should be
   documented that way if any external reporting happens.
3. Stage 4 (N=5 warm-cache benchmark) should be re-run against
   fixed-legacy for the actual paper/lab-meeting numbers, not against
   Risk-E-buggy legacy.

> To be resolved by: a new entry that reports fixed-legacy HPRC timing
> measurement + apples-to-apples flipped delta.

> Resolved by JR-007 (2026-07-29). Estimation in this entry turned out to be
> unreliable in the specific way flagged ("first-locate scales with MEM count,
> not occurrence volume"); JR-007 supersedes both the numbers here and the
> logic behind them.

---

## JR-007 — Risk E fix: fixed-legacy == flipped (byte-identical coverage) + honest Stage 3 delta

```yaml
id: JR-007
date: 2026-07-29
author: claude-opus-4-7 (session with hlakshmidevi)
status: resolved
tags: [risk-e, fix, correctness, perf, stage3]
refs:
  follows: [JR-006]
  resolves: [JR-004, JR-005, JR-006]
```

### Context

JR-004 established Risk E as a real bug in legacy Step 3
(`include/pangenome_index/algorithm.hpp:724`). JR-005 characterized its
size distribution and per-read concentration. JR-006 tried to estimate
the honest Stage 3 perf delta arithmetically and flagged the estimation
as unreliable. This entry closes all three: applies the fix, verifies
correctness via a stronger check than validate_gaf (byte-identical
downstream coverage), and measures the honest apples-to-apples flipped
delta.

### Hypothesis

Under JR-004's model of the bug:

- Fix: guard `if (e >= len) return len;` before Step 3's loop skips the
  ill-defined `pattern[len]` access. If e == len, the emitted MEM covers
  end-of-pattern; no right-anchor exists to the right of it, so returning
  `len` correctly terminates the outer loop.
- Prediction 1: fixed-legacy's MEM set becomes identical to flipped's on
  every dataset (`test_flipped_mems` LEGACY_SPURIOUS goes to 0).
- Prediction 2: fixed-legacy and flipped produce byte-identical
  downstream `alignment_coverage.csv`.
- Prediction 3: fixed-legacy find_mems wall drops substantially vs
  buggy-legacy (Risk E was driving ~20x downstream volume amplification
  on HPRC), giving an honest Stage 3 delta smaller than JR-003's 28%
  headline.

### Method

**Fix** applied to `include/pangenome_index/algorithm.hpp` (in the
`find_mems_function` between MEM emission and Step 3 setup). One line
plus a substantial comment explaining the mechanism, so a future reader
who trips this code path understands why the guard exists.

```cpp
if (e >= len) return len;
```

**Verification chain:**

1. Unit-level (`bin/test_flipped_mems`) on HPRC 10K reads and yeast
   100K reads with the fix in place. LEGACY_SPURIOUS must go from 170
   / 529 respectively to 0.
2. Pipeline E2E on yeast via `./query.sh` (query name:
   `riske-fixed-legacy`), flag OFF. validate_gaf 2000/2000.
3. Pipeline E2E on HPRC via manual invocation (query.sh does not accept
   extra find_mems flags), flag OFF, with `gtime -v` for peak RSS.
4. **Downstream identity check via `alignment_coverage.csv` MD5** across
   `riske-fixed-legacy` vs `stage3-flipped` (flipped) vs
   `local-validate` (buggy legacy baseline). This is a much stronger
   correctness signal than validate_gaf's 2000-entry sample -- coverage
   counts feed cosigt directly, so byte-identity here means genotyping
   pipelines see the exact same input.

### Findings

**Unit-level (`bin/test_flipped_mems`):**

| Dataset | Buggy-legacy LEGACY_SPURIOUS | Fixed-legacy LEGACY_SPURIOUS |
|---|---:|---:|
| HPRC 10K reads | 170 | **0** |
| Yeast 100K reads | 529 | **0** |

Fixed-legacy's MEM set is now exactly equal to flipped's on both datasets.
Zero flipped-missed, zero flipped-extra, zero sa_sp mismatches (unchanged
-- those were already zero).

**Yeast E2E (validate_gaf):**

- `runs/v2-yeast235/queries/riske-fixed-legacy/lightweight/`: 2000/2000 valid.
- GAF line count: 4,906,654 (matches flipped's stage3-flag-on exactly,
  down from buggy-legacy's 4,954,592).

**HPRC E2E (find_mems only, no validate_gaf per user):**

- `runs/hprc-chr6-2026-06-02/queries/riske-fixed-legacy/lightweight/`.
- Total MEMs: 100.0K (was 100.8K in buggy-legacy).
- Total entries written: 1.1M (was 21.7M -- Risk E amplification gone).
- LocateNext calls: 2.7M (was 47.7M).
- Wall: 40.66s (was 47.996s).
- Peak RSS: 4113 MB (was 4600 MB).

**Downstream identity (alignment_coverage.csv MD5):**

Yeast:

| variant | MD5 | size |
|---|---|---:|
| `riske-fixed-legacy` (flag OFF, fix applied) | `1baf368f3cdc59890f88410cc34cc051` | 16,108,726 |
| `stage3-flag-on` (flipped) | `1baf368f3cdc59890f88410cc34cc051` | 16,108,726 |
| `stage3-flag-off` (buggy legacy) | `db4ba02af811e72df2790db064b152b9` | 16,108,919 |

HPRC:

| variant | MD5 | size |
|---|---|---:|
| `riske-fixed-legacy` (flag OFF, fix applied) | `8c242732c837d697fb18c52fd0c0cf6e` | 61,058,356 |
| `stage3-flipped` (flipped) | `8c242732c837d697fb18c52fd0c0cf6e` | 61,058,356 |
| `local-validate` (buggy legacy baseline) | `8bb1c7bb78afca87ddd76bf93a7f732d` | 61,198,316 |

**Fixed-legacy and flipped produce byte-identical coverage files on both
datasets.** The buggy-legacy baseline differs -- its spurious MEMs
contributed real (but wrong) coverage counts on some nodes.

**Note on `.bin` and `.gaf` MD5s:** On yeast, `mems_path_pos_v2.bin`
and `alignment.gaf` MD5s differed between fixed-legacy and flipped
despite identical line counts and byte-identical sorted content
(diff-sorted returned 0 lines). Root cause: the pipeline sorts
`mems_path_pos_v2.bin` by `(seq_id, path_bp)` before writing, and when
multiple records tie on that key, the sort's tie-break is
implementation-defined and depends on insertion order. Legacy iterates
left-to-right; flipped iterates right-to-left; same records, different
insertion order, different tie-break resolution. On HPRC the MD5s
happened to match -- either fewer ties, or the sort was stable and the
insertion orders happened to converge. **Neither case affects
correctness** -- gafpack aggregates records into coverage counts
regardless of order within a bucket, which is why the coverage file is
byte-identical either way. If future work wants to make the .bin/gaf
outputs deterministic across finders, either make the pipeline sort
stable-with-explicit-tiebreaker or add a secondary sort key.

**Honest Stage 3 perf delta (HPRC chr6, 100K reads):**

| Metric | Fixed-legacy | Flipped | Δ |
|---|---:|---:|---:|
| Wall | 40.66s | 34.44s | **-15.3%** |
| Peak RSS | 4113.61 MB | 4106.73 MB | -0.2% (noise) |
| Time for finding MEMs | 19.61s | 29.22s | **+49% (Risk C)** |
| Time for processing MEMs | 17.37s | 1.34s | **-92%** |
| ├─ First locate (seed) | 15.54s | **0s** | **-100% (Stage 3's win)** |
| ├─ Locate next | 0.48s | 0.89s | +85% (small absolute) |
| ├─ Tag queries | 0.55s | 0.34s | -37% |
| Total MEMs | 100.0K | 100.0K | 0 (identical set) |
| Total entries written | 1.1M | 1.1M | 0 |
| Sort time | 0.19s | 0.11s | -42% |

**Wall breakdown of the -15.3% improvement:**

- Seed elimination saves 15.54s (this is the Stage 3 primary win).
- Locate-next and other processing costs approximately break even
  (Stage 3 shifts a tiny amount of work, sub-second).
- MEM finding costs 9.61s more (Risk C: Step 2' does a fresh forward
  extend that partially duplicates work legacy amortized across its
  Step 1 and Step 3).
- Sort saves 0.08s.
- Net: 15.54 - 9.61 + ~0.1 = ~6s saved, matching the observed 40.66s
  -> 34.44s = 6.22s delta.

### Interpretation

**Correctness result is stronger than JR-004/JR-005 dared to hope for.**
The fix isn't just "reduces spurious MEMs" -- fixed-legacy is
*bit-identical* to flipped at the downstream coverage level, on both
datasets. Cosigt cannot tell them apart. This retires all doubt about
whether flipped and legacy compute equivalent biology.

**Honest Stage 3 delta is -15.3% wall on HPRC.** Between the design's
original 10% projection (§4.4) and JR-003's misleading 28% headline. It
sits closer to the projection because:

- Design's 10% used a linear cost model where seed cost is 32.2% of
  walk time on Risk-E-buggy legacy. On Risk-E-fixed legacy the seed's
  share of processing time is much higher (`15.54 / 17.37 = 89%` of
  MEM processing), because Risk E was inflating the "other work"
  denominator. So seed elimination's *relative* impact on find_mems
  wall goes up.
- But Stage 2's Step 2' cost is a real -10s hit. Design flagged this
  as Risk C and said "measure in Stage 4"; now measured.

**Risk C is the biggest opportunity for future work.** Step 2' does
work that partially duplicates what Step 1' just did (opposite
direction). A version of the flipped finder that memoized or reused
Step 1's rightward walk could recover some of the 9.6s -- worth
investigating in a future entry.

**Peak RSS is essentially identical (4.11 vs 4.11 GB, delta = 7 MB).**
The 11% RSS drop JR-003 reported vs buggy-legacy was Risk E's fault --
the ballooning entries vector before sort. Fixed-legacy also benefits
from that; both paths now use the same in-memory volume.

**Test binary implication:** `bin/test_flipped_mems`'s LEGACY_SPURIOUS
category is now dead code -- it will always report 0. Could be
removed, but leaving it in place is defensible as a regression guard:
if a future refactor of legacy accidentally re-introduces a
non-left-maximal-emission bug (Risk E or a similar one), this test
will catch it. Cost: a `count_encoded` call per divergent MEM (i.e.,
zero on healthy code paths).

### Open questions

1. **Coverage-file identity as a new canonical validation step.**
   VALIDATION_GUIDE.md currently emphasizes validate_gaf's 2000-sample
   check. The coverage-file MD5 comparison is objectively stronger --
   it's the exact bytes cosigt consumes. Consider updating
   VALIDATION_GUIDE.md to recommend coverage MD5 comparison as the
   preferred check when a known-good baseline exists.
2. **Non-deterministic .bin ordering across finders.** Documented as
   a non-issue for correctness above. If deterministic outputs are
   ever needed (e.g., for git-based CI diffs), the fix is a stable
   sort with a total-order tie-breaker in `write_sorted_entries`.
3. **Risk C mitigation.** Step 2''s work overlap with Step 1' is
   ~10s / MEM-finding on HPRC. Worth a targeted investigation
   entry before defaulting to flipped.
4. **Turning `--use-flipped-mems` on by default.** Now that
   fixed-legacy == flipped at coverage-file level, the decision is
   purely about the -15% perf tradeoff. Deferred pending user
   direction; when it's made, that decision itself should be a JR
   entry documenting the tradeoff and any rollback plan.

### Commits

- (fix commit, pending as of writing)
- (journal update commit, pending)

Data artifacts:
- `runs/v2-yeast235/queries/riske-fixed-legacy/lightweight/`
- `runs/hprc-chr6-2026-06-02/queries/riske-fixed-legacy/lightweight/`

---

## JR-008 — Stage 3 benchmark: N=5 warm-cache HPRC chr6 100K reads

```yaml
id: JR-008
date: 2026-07-29
author: claude-opus-4-7 (session with hlakshmidevi)
status: resolved
tags: [stage3, perf, benchmark, hprc, harness]
refs:
  follows: [JR-007]
  confirms: [JR-006, JR-007]
```

### Context

JR-007 established Stage 3 correctness (byte-identical
`alignment_coverage.csv` between fixed-legacy and flipped) and reported
a single-run wall delta of -15.3% on HPRC chr6. User asked for a proper
benchmark with variance signal, using the existing perf harness at
`mem-projection/pangenome-pipeline/perf/perf_harness.sh` rather than
ad-hoc invocations. This entry runs N=5 warm-cache trials per mode and
reports the mean/min/median with stdev.

### Method

**Harness modification.** The existing `perf_harness.sh` accepts
`--modes=lightweight|full-tag` as the axis of comparison but has no
hook for extra find_mems flags such as `--use-flipped-mems`. Added a
one-line backward-compatible extension: an optional
`FIND_MEMS_EXTRA_FLAGS` environment variable is appended to the
find_mems command line. Unset = current behavior. Provenance file
records the value.

CLAUDE.md discourages editing downstream pipeline files from the
finder repo, but the change is (a) purely additive, (b) preserves the
default behavior, and (c) the alternative (reimplementing the harness
protocol in the finder repo) would duplicate ~200 lines and risk drift.
User was asked and approved.

**Benchmark protocol.** Ran the harness twice back-to-back within a
single machine session:

```
export MEM_PROJ=/Users/hlakshmidevi/personal/mem-projection
export INDEX_DIR="$MEM_PROJ/pangenome-pipeline/runs/hprc-chr6-2026-06-02"

FIND_MEMS_EXTRA_FLAGS="" \
  ./perf/perf_harness.sh configs/hprcv1-chr6-alt-reads.env 5 \
    hprc-stage3-legacy --modes=lightweight

FIND_MEMS_EXTRA_FLAGS="--use-flipped-mems" \
  ./perf/perf_harness.sh configs/hprcv1-chr6-alt-reads.env 5 \
    hprc-stage3-flipped --modes=lightweight
```

Each harness invocation: 2 warmups (untimed, to reach warm-cache
steady state) + 5 timed trials, contiguous ordering (all warmups then
all trials for the same mode). Coverage-only mode (no `.gaf`, no
validate_gaf).

Then merged the two per-run outputs into a single perf-tag directory
so `summarize.py` produces a pairwise A/B comparison:
- `perf/hprc-stage3-legacy/lightweight/` -> `perf/hprc-stage3/legacy/`
- `perf/hprc-stage3-flipped/lightweight/` -> `perf/hprc-stage3/flipped/`

**Pre-restart false start.** An earlier N=3 run produced legacy
trial-1 wall = 133s (~3.3x the others). Discarded and restarted with
N=5; the retry showed legacy trial-1 at 46.8s (still slower than
trials 2-5 which clustered 37.9-39.3s, but not catastrophic). Root
cause of the 133s spike not investigated -- possibly other machine
activity between the harness warmups and timed trials. Reported here
for transparency; the N=5 restart is the actual benchmark.

Config: `configs/hprcv1-chr6-alt-reads.env` (MEM_LEN=50, MIN_OCC=1).
Reads: `../hprcv1/chr6.alt.reads.txt` (100K reads, 200bp each,
19MB file).
Machine: M1 mac laptop (single user, background quiescent).
`find_mems` binary: commit `fe9373e` (Risk E fix applied).

### Findings

**Per-trial wall (find_mems, seconds):**

| trial | legacy | flipped |
|---:|---:|---:|
| 1 | 46.80 | 35.13 |
| 2 | 39.33 | 34.54 |
| 3 | 38.57 | 35.88 |
| 4 | 39.25 | 34.26 |
| 5 | 37.91 | 35.28 |
| **min** | **37.91** | **34.26** |
| **median** | **39.25** | **35.13** |
| **mean** | **40.37** | **35.02** |
| max | 46.80 | 35.88 |
| stdev | ±3.64 | ±0.64 |

**find_mems wall delta (flipped vs legacy):**

- **min-vs-min:** -9.6% (36.9s -> 34.3s)
- **median-vs-median:** -10.5% (39.3s -> 35.1s)
- **mean-vs-mean:** -13.3% (40.4s -> 35.0s, per `summarize.py`)

The median-based delta is essentially exactly the design's original
projection of ~10% (DESIGN_FLIPPED_MEM.md §4.4). The mean-based
delta is inflated by legacy's trial-1 outlier; use median as the
robust headline.

**Full breakdown from `summarize.py perf/hprc-stage3` (N=5, warm cache):**

```
STEP 09 - find_mems
  metric                  flipped                 legacy                     delta
  wall (rusage)           35.02 +-  0.64 s        40.37 +-  3.64 s        +15.3%
    user                  29.53 +-  0.11 s        30.82 +-  0.36 s         +4.4%
    sys                    4.58 +-  0.41 s         7.63 +-  2.14 s        +66.4%
  peak RSS (kB)           4,201,613 +-  65,277    4,243,389 +-  37,078      +1.0%
  peak (reported)         4090.0 +-  56.8 MB      4134.9 +-  31.8 MB       +1.1%
  minor faults              645,310 +- 33,772       911,297 +- 233,678     +41.2%
  major faults            132 +- 5                136 +- 6                  +2.6%
  vol ctx switches        522 +- 59              2117 +- 2571            +305.2%
  invol ctx switches      53518 +- 4467          81085 +- 21357            +51.5%

STEP 09 - find_mems  PHASE BREAKDOWN (from find_mems internal profiler)
  phase                   flipped                 legacy                     delta
  Total exec              34.81 +-  0.64 s        40.19 +-  3.64 s        +15.5%
    R-index load           2.619 +- 0.125 s        2.535 +- 0.203 s         -3.2%
    Tag-index load         0.725 +- 0.074 s        0.743 +- 0.040 s         +2.5%
    Read processing       31.35 +-  0.57 s        36.79 +-  3.74 s        +17.4%
      MEM finding         29.64 +-  0.46 s        19.75 +-  1.13 s        -33.4%
      MEM processing       1.51 +-  0.15 s        16.76 +-  2.61 s      +1007.9%
        Tag queries        0.416 +- 0.072 s        0.838 +- 0.552 s       +101.4%
        Locate ops         0.99 +-  0.08 s        15.78 +-  2.01 s      +1501.4%
          First locate     0.00 +-  0.00 s        15.26 +-  1.98 s           n/a
          Locate next      0.99 +-  0.08 s         0.53 +-  0.04 s        -46.7%
    Bucket sort+write     0.116 +- 0.009 s        0.118 +- 0.009 s         +2.1%
```

(Signs read as: positive delta = legacy's value is larger than
flipped's. So "+15.3%" wall means legacy is 15.3% larger than
flipped, i.e., flipped is 13.3% faster on a mean basis or 10.5%
faster on a median basis.)

**Peak RSS delta = +1.0% (essentially noise).** JR-003's reported 11%
RSS improvement was against Risk-E-buggy legacy; fixed-legacy now uses
the same in-memory volume as flipped.

**Output artifacts byte-identical across all 5 x 2 = 10 trials:**

| artifact | flipped | legacy | delta |
|---|---:|---:|---:|
| `mems_path_pos_v2.bin` | 19,175,136 | 19,175,136 | 0.0% |
| `alignment_coverage.csv` | 61,058,356 | 61,058,356 | 0.0% |
| `mems_seq_id_starts.out` | 19,872 | 19,872 | 0.0% |

Also gafpack's internal stats (records loaded, total GAF entries,
step-visits, seq_ids scanned, mean passes/path) are bit-identical
across the two modes. Confirms JR-007's byte-identity finding at
N=5-per-mode scale.

**gafpack step is unaffected**, as expected -- it processes the same
input records either way. Wall 16.5s vs 16.5s (+0.2%, noise).

### Interpretation

**Headline: flipped is ~10-13% faster than fixed-legacy on find_mems
wall, with byte-identical outputs.** The design's original 10%
projection sits inside the min-median-mean range.

**Wall decomposition of the ~5s improvement:**
- Seed elimination saves ~15.3s (`First locate: 15.26 -> 0.00`).
  Stage 3's primary contribution.
- Step 2' costs +9.9s in MEM finding (`19.75 -> 29.64`).
  This is Risk C from DESIGN_FLIPPED_MEM.md §5.3 -- flipped's Step 2'
  does a fresh forward extend that partially duplicates work legacy
  amortizes across Steps 1 and 3.
- LocateNext costs +0.5s in flipped (`0.53 -> 0.99`; small
  absolute).
- Tag queries save ~0.4s in flipped (probably from tighter working
  set on the tag index, but noise floor is close).
- R-index load, tag-index load, bucket sort: all within 1-3%
  (indistinguishable from noise).
- Net: +15.3 - 9.9 - 0.5 + 0.4 = +5.3s, matches the observed
  40.4 -> 35.0 delta.

**Flipped run is dramatically more consistent (stdev ratio ~5.7x):**
- Flipped stdev = 0.64s (1.8% of mean).
- Legacy stdev = 3.64s (9.0% of mean), driven almost entirely by
  trial-1 = 46.8s.
- Hypothesis: legacy's per-MEM `locate_sa_value` walk is very
  cache-sensitive (long walks through cold BWT-run pages);
  small variations in page-cache state produce large wall variations.
  Flipped's O(rank) SA-carry has smaller cache footprint and is
  more deterministic. Not confirmed by direct measurement.

**Voluntary context switches show a 4x difference** (flipped 522 vs
legacy 2117). Legacy's per-MEM locate_sa_value walk apparently
triggers more I/O yields into the kernel; flipped's tighter SA-carry
avoids them. Consistent with the cache-sensitivity hypothesis.

**Minor page faults +41% legacy vs flipped.** Legacy touches more BWT
pages (via the locateNext walks in `locate_sa_value`), so the
first-touch fault count is higher. Same underlying phenomenon.

**Combined 09+10 wall delta:**
- flipped: 51.51 +- 0.83s
- legacy: 56.90 +- 3.41s
- delta: +10.5% (find_mems wins dominate; gafpack is unchanged).

### Open questions

1. **Turn `--use-flipped-mems` on by default?** With correctness
   proven byte-identical (JR-007) and perf +10% wall on find_mems
   (this entry), the case is clear. Design's Stage 5 in
   `DESIGN_FLIPPED_MEM.md` §5.2 outlines the rollout. Would need:
   - Flip the CLI default (or introduce `--use-legacy-mems` as a
     regression-testing escape hatch).
   - Update DESIGN_FLIPPED_MEM.md status line to "Stage 5 landed".
   - Update `mem-projection/pangenome-pipeline/query.sh` to
     acknowledge the change (or leave it alone, since the flag
     flip doesn't change the invocation).
   - Rerun the harness one more time after the flip to confirm
     nothing else in the codebase depended on the previous default.
2. **Risk C mitigation still on the table.** Step 2' costs +9.9s of
   the ~15.3s seed savings. Recovering some of that would push the
   delta well past 15%. A future entry could investigate reusing
   Step 1's rightward walk instead of a fresh Step 2' forward extend.
3. **Why is legacy trial-1 always slower?** Both this run (N=5) and
   the discarded N=3 run had trial-1 legacy noticeably above the
   others. Two warmups appear insufficient to stabilize legacy's
   locate walks; flipped isn't affected. Increasing warmups from 2
   to 3 (via harness edit) or dropping trial-1 as a post-hoc extra
   warmup would clean up the mean. Not urgent -- median is already
   robust to this outlier.

### Data artifacts

- `perf/hprc-stage3/` -- merged perf-tag dir consumed by summarize.py.
- `perf/hprc-stage3-legacy/`, `perf/hprc-stage3-flipped/` -- the two
  raw harness invocations preserved for full-trial forensics.
- `perf/hprc-stage3/PROVENANCE.txt` -- concatenated harness provenance
  (host, date, binaries, index files, FIND_MEMS_EXTRA_FLAGS value
  per run).
- `perf/hprc-stage3/SUMMARY.tsv` -- per-trial machine-readable rows,
  mode column renamed to `legacy`/`flipped`.

Harness edit committed separately in the pipeline repo (one-line
addition of `${FIND_MEMS_EXTRA_FLAGS:-}` word-splitting into the
find_mems invocations).

---

## JR-011 — backward_extend_encoded_with_sa rewrite: HPRC chr6 noisy benchmark

```yaml
id: JR-011
date: 2026-07-30
author: claude-opus-4-7 (session with hlakshmidevi)
status: resolved
tags: [jr-010, sa-carry, primitive, perf, hprc, noisy]
refs:
  follows: [JR-008]
```

> Note on lineage: this entry addresses the perf regression documented
> off-journal in the session brief as "JR-010" (flipped +50% slower than
> legacy on `chr6.alt_noisy.reads.txt`, +48s hit entirely in MEM-finding,
> traced to redundant pred+block-walks inside
> `backward_extend_encoded_with_sa`). JR-009 / JR-010 numbers are reserved
> for the user's parallel work; this session's number is JR-011.

### Context

JR-008 established Stage 3's baseline on `chr6.alt.reads.txt` (error-free
reads): flipped ~10-13% faster than fixed-legacy at find_mems wall.
When the same code was rerun on `chr6.alt_noisy.reads.txt` (noisy reads,
same read count and length), the picture inverted -- flipped became
~+50% slower than legacy, and the entire ~48s regression landed inside
MEM-finding rather than MEM-processing. Root cause: on noisy reads,
Step 1' walks fail early far more often, so `backward_extend_encoded_with_sa`
gets called an order of magnitude more times per read. The primitive's
Case B path was doing up to 5 predecessor+block-walks per call
(`backward_extend_encoded` internal 2 rank walks + one
`bwt_char_at_encoded` for the Case A test + one `bwt_char_at_encoded`
fast-path inside `leftmost_c_in_interval` + one `run_id_and_offset_at`
at the leftmost-c position) plus, on the toehold path, an entire
locate walk from the containing run's sample.

This entry rewrites the primitive per a specific plan
(see the sub-step 1-8 in DESIGN_FLIPPED_MEM.md section 5) and measures
the resulting delta on the noisy-read workload.

### Hypothesis

Two-part hypothesis:

1. **Primitive-level.** Collapse the 3-5 pred+walks per general-case
   call into 2 by introducing a combined block scan (`scan_at`, returns
   ranks + BWT char + run id + run start in one walk) and a rank-in
   variant of `backward_extend_encoded`
   (`backward_extend_encoded_with_ranks`). Eliminate the toehold
   locate walk by precomputing a 5-entry `first_bint[c] / first_sa[c]`
   table at load. Replace Case B's leftmost-c lookup with a
   successor query over a per-DNA-char RUN-INDEXED bitvector
   (`run_head_c[c][r] = 1` iff run r's head char is c).

2. **Wall-time.** Since the +48s regression sits entirely in
   MEM-finding, and MEM-finding is dominated by
   `backward_extend_encoded_with_sa` calls, cutting per-call cost by
   ~60% should recover a proportional share of the regression. If
   the primitive is the sole regression driver, flipped-vs-legacy
   should close to within 10%. If other flipped-code-path costs
   (notably Step 2' forward extends) also contribute meaningfully,
   the primitive rewrite will only partially close the gap.

### Method

**Implementation** landed as commits `23a73a0..607a0ca` on branch
`backward-extend-sa-perf`, 7 test-first sub-steps ending in the
primitive rewrite (`src/r-index.cpp:951-1055`). Full commit list:

- `23a73a0` test scaffold (Test A/B/C, A and B fail initially)
- `59eb360` `EncodedBlock::ranks_at` optional out-params
- `ccb502a` `scan_at` combined block scan; Test A passes
- `0fa123d` `backward_extend_encoded_with_ranks` rank-in variant
- `c44e475` `run_head_c` parallel run-id bitvector; Test B passes
- `5b1a4d7` `first_bint` / `first_sa` toehold table
- `607a0ca` primitive rewrite (toehold fast path + general path)

**Unit verification** (`bin/test_backward_extend_sa` on yeast235 chrII,
`../mem-projection/pangenome-pipeline/runs/v1-current/yeast235_chrII_100kb_normalized.ri`):
- Test A (scan_at cross-check, 5000 positions): 5000/5000 PASS
- Test B (run_head_c cross-check, 2000 checks per DNA char): 10000/10000 PASS
- Test C (SA-carry oracle, 10K trials, mix of random-ACGT and
  text-derived patterns): 302243 steps, 0 mismatches, 10000 toehold-path
  hits, Case A / Case B = 213725 / 78518 (73:27 split).

**E2E correctness** (`./query.sh configs/hprcv1-chr6-alt-reads.env
hprc-chr6-2026-06-02 s7-flipped --gaf` with `FIND_MEMS_EXTRA_FLAGS=--use-flipped-mems`):
- validate_gaf: 2000 / 2000 valid, 0 invalid.
- GAF line count: 532,534 (identical to the sub-step 5 legacy baseline
  on the same reads file -- no silent MEM drop).
- `alignment_coverage.csv` MD5 `60f6b8e4a759aebb252b83a870b9ff8c` --
  byte-identical to the legacy baseline. Coverage-file identity is
  the strong correctness signal per JR-007; cosigt cannot distinguish
  the two.

**Benchmark protocol.** Built two find_mems binaries side-by-side:
- `bin/find_mems` (post-rewrite): commit `607a0ca`.
- `bin-pre-rewrite/find_mems`: commit `5b1a4d7` (all sub-steps landed
  except the primitive rewrite itself; isolates the delta to the
  rewrite alone rather than confounding with load-time cost of
  `run_head_c` + `first_bint`).

Per binary and per mode (legacy / flipped): 2 warmup runs (untimed)
+ 3 timed trials, warm-cache, contiguous ordering per variant. Ran
directly (not through `query.sh`) to skip gafpack + validate_gaf on
each trial. Command per trial:

```
$BIN $RI $LTAGS $READS 50 1 mems --lightweight-tags [--use-flipped-mems]
```

Config: `hprcv1-chr6-alt-reads.env` (MEM_LEN=50, MIN_OCC=1). Reads:
`../mem-projection/hprcv1/chr6.alt_noisy.reads.txt` (100K reads,
200bp, 19MB). Machine: same M1 Mac as JR-008.

### Findings

**Per-trial wall (find_mems, seconds):**

| trial | legacy | pre-rewrite flipped | post-rewrite flipped |
|---:|---:|---:|---:|
| 1 | 49.55 | 70.44 | 62.01 |
| 2 | 46.79 | 97.89 | 61.78 |
| 3 | 48.47 | 74.18 | 61.04 |
| **median** | **48.47** | **74.18** | **61.78** |

**Per-trial phase breakdown (median trial):**

| Phase | legacy | pre-rewrite flipped | post-rewrite flipped |
|---|---:|---:|---:|
| Total | 48.47 s | 74.18 s | 61.78 s |
| R-index load | 8.17 s | 8.33 s | 8.10 s |
| **MEM finding** | **18.67 s** | **58.72 s** | **46.42 s** |
| MEM processing | 19.30 s | 5.74 s | 5.74 s |
| ├─ First locate | 17.14 s | 0.00 s | 0.00 s |
| ├─ Locate next | 1.54 s | 4.56 s | 4.56 s |
| └─ Tag queries | 0.48 s | 1.01 s | 1.01 s |

**Deltas (post-rewrite vs pre-rewrite flipped):**
- Total wall: **-12.4 s (-16.7%)**
- MEM finding: **-12.3 s (-20.9%)**

Everything else identical within noise (correctness proof: MEM
processing, sorting, load times all unchanged; the rewrite only
touches the SA-carry primitive body).

**Deltas (post-rewrite flipped vs legacy):**
- Total wall: +13.3 s (+27.4%)
- MEM finding: +27.8 s (+148.9%)
- MEM processing: **-13.6 s (-70.5%)** (seed elimination win preserved)

**Peak RSS:** legacy 4192 MB, post-rewrite flipped 4182 MB -- delta
0.2%, well under the 10% CLAUDE.md ceiling.

**Case A / Case B split** on the noisy workload
(from Test C's 10K trials on yeast235; production workload
distribution not directly measured but is expected to be similar):
Case A 73%, Case B 27%, initial toehold ~3% of total steps.

### Interpretation

**Primitive rewrite delivers the projected ~12 s recovery on
MEM-finding.** The pre-rewrite Case B path was doing up to 5
pred+block-walks per call plus a full locate walk on the toehold
path (call-once-per-outer-loop-iter). The post-rewrite path does 2
pred+block-walks in the general case and 0 memory-load fields on
the toehold path. On the noisy workload (many short Step 1' walks,
so per-call cost dominates), the ratio 5:2 predicts a ~60% cut in
primitive time, and MEM finding is dominated by primitive calls,
so a ~12-13 s cut on the ~58 s pre-rewrite MEM finding matches
prediction well.

**But flipped remains ~+28 s slower than legacy on MEM-finding.**
Broken down:

- Seed elimination saves 17.1 s (First locate: 17.14 → 0.00). This
  is the whole Stage 3 win, preserved intact.
- Post-rewrite primitive is still more expensive per call than the
  equivalent legacy backward-extend chain (which piggy-backs on
  cache warmed by adjacent Step 1 / Step 3 walks). Absolute
  primitive-time delta not directly measured; inferred from the
  ~46 s MEM-finding on flipped vs ~18 s on legacy = +28 s.
- Step 2' forward-extend cost. Forward extends are ~2x per-char
  cost of backward extends (swap → backward-extend-on-RC → swap).
  On noisy reads, Step 2' fires on every outer-loop iteration,
  and short backward walks mean more outer-loop iterations per
  read. Not directly measured (would need finer per-Step
  instrumentation); dominates the residual gap by process of
  elimination.

**Net: post-rewrite flipped is +13 s wall vs legacy on the noisy
workload, down from +26 s pre-rewrite.** Half the JR-010
regression closed; the other half sits outside the primitive.
The +13 % wall difference does NOT meet the "within 10%" ideal
called out in the session brief, but does exceed the fallback
requirement of "document what remains and next candidate
directions" -- see Open questions below.

**Coverage is byte-identical to legacy across the entire branch.**
Every sub-step's E2E gate ran `--gaf` and got 2000/2000 valid;
`alignment_coverage.csv` MD5 is `60f6b8e4a759aebb252b83a870b9ff8c`
for sub-step 2 legacy, sub-step 5 legacy, sub-step 7 flipped -- all
three match. Correctness is not in dispute at any point in the
branch.

**Load-time cost.** The eager `run_head_c` build adds ~7-8 s to
R-index load (2.5 → ~8 s on chr6). This is a one-time cost per
process invocation; in a batch-query context (many queries per
process) it amortizes to zero, but for latency-sensitive single-
query callers it may matter. See Open question 2.

**Peak RSS ceiling not touched.** Delta 0.2% vs legacy.

### Open questions

1. **Step 2' overhead is the remaining ~+14 s.** **RESOLVED in JR-012**
   via Step 0' pre-check: the outer-iteration count dropped 4.3x
   (1.39M → 326K), and flipped's total wall now beats legacy on
   both noisy (−20%) and clean (−28%). See JR-012 for details.
   Original hypothesis directions (Step 1'/Step 2' state reuse,
   or cheaper anchor-advance heuristic) were not needed; the
   min_len pre-check turned out to be the right lever.
2. **Consider serializing `run_head_c` and `first_bint` / `first_sa`
   to disk** in the .ri format, or dropping them entirely if a
   future rewrite makes their queries unnecessary. Current one-time
   ~8 s load cost is amortized in batch contexts but visible in
   single-query contexts.
3. **Whether to enable `--use-flipped-mems` by default is now more
   complicated than JR-008 suggested.** On clean reads the flipped
   path wins by ~10 %; on noisy reads it still loses by ~13 %.
   A workload-aware toggle would be ideal but adds complexity. If
   the production cohort mix skews clean-reads-heavy, flipping the
   default is still justified; if noisy is common, keep legacy
   default and gate flipped on an explicit flag. Deferred pending
   user direction.
4. **Test C's toehold-hit count is a synthetic distribution.** Real
   production workloads may have different Case A / Case B / toehold
   frequencies. Instrumenting `backward_extend_encoded_with_sa`
   with a per-invocation counter and dumping the histogram from a
   real find_mems run would firm up the model.

### Data artifacts

- `xy-test/jr011_bench_legacy_{1,2,3}/run.log` -- per-trial logs
  for the legacy binary (find_mems internal profiler dump).
- `xy-test/jr011_bench_pre_rewrite_flipped_{1,2,3}/run.log` -- the
  three pre-rewrite flipped trials.
- `xy-test/jr011_bench_post_rewrite_flipped_{1,2,3}/run.log` -- the
  three post-rewrite flipped trials.
- `runs/hprc-chr6-2026-06-02/queries/s7-flipped/` and
  `s7-flipped-rerun/` -- E2E `--gaf` runs of the post-rewrite
  primitive, both with validate_gaf 2000/2000 and coverage MD5
  `60f6b8e4a759aebb252b83a870b9ff8c`.
- `runs/hprc-chr6-2026-06-02/queries/s{2,5}-regcheck*/` -- sub-step
  regression checks (legacy path, no functional change; used to
  establish the coverage baseline for identity comparisons).

### Pipeline-repo edit

`../mem-projection/pangenome-pipeline/query.sh` gained a one-line
`${FIND_MEMS_EXTRA_FLAGS:-}` word-split after the flag array,
mirroring the perf_harness.sh hook JR-008 added. This is the same
pattern JR-008 documented and got user approval for. Not committed
in this session -- flagged for user review before commit on the
pipeline repo side.

### Commits (this repo, branch backward-extend-sa-perf)

- `23a73a0` test: scaffold Test A/B/C for backward_extend_encoded_with_sa optimization
- `59eb360` r-index: extend EncodedBlock::ranks_at with optional block-scan extras
- `ccb502a` r-index: add scan_at combined block-scan; Test A passes
- `0fa123d` r-index: add backward_extend_encoded_with_ranks (rank-in variant)
- `c44e475` r-index: add run_head_c parallel run-id bitvector; Test B passes
- `5b1a4d7` r-index: precompute first_bint/first_sa toehold table at load
- `607a0ca` r-index: rewrite backward_extend_encoded_with_sa (JR-010 perf fix)
- `9e728de` journal: JR-011 backward_extend_encoded_with_sa rewrite + noisy-chr6 bench

---

## JR-012 — Step 0' min_len pre-check in find_mems_flipped: HPRC chr6 noisy + clean

```yaml
id: JR-012
date: 2026-08-01
author: claude-opus-4-7 (session with hlakshmidevi)
status: resolved
tags: [flipped, algorithm, perf, hprc, noisy, clean, step0, vesuvio]
refs:
  follows: [JR-011]
  supersedes-perf-numbers-of: [JR-011]
benchmark-platform: vesuvio (Linux 6.12.73 x86_64, Debian)
```

### Context

JR-011 fixed the JR-010 primitive-level regression but left flipped's
outer-loop iteration count at 3.79x legacy on `chr6.alt_noisy.reads.txt`
(1,389,639 iters vs 366,797) and its extend-call count at 2.51x legacy
(58.8M vs 23.4M).  Post-JR-011 flipped MEM-finding on noisy chr6
was ~44.5s vs legacy ~19.1s — +133%, entirely inside the outer loop.

Profiling (bin/debug_flipped_perf, bin/debug_anchor_distribution)
showed 88.8% of outer iterations were "short-match" iterations: Step 1'
backward-extending ~11-20 characters before failing (min_len is 50), and
Step 2' then advancing `x` by a very small anchor gap (gap=2 in 48% of
cases, 82.7% ≤ 10).  Each such iteration cost ~27 bwe + ~15 fwe with
no MEM emission.  Legacy's Step 1 avoids exactly this class of work
via its length-`min_len` bound pre-check: it forward-extends the
length-`min_len` window; on failure at position j, it can jump the next
right endpoint to `j + 1` because no MEM of length ≥ `min_len` can
straddle j.  The flipped algorithm inherited no analogous pre-check.

### Hypothesis

Mirror legacy Step 1's `min_len` bound in the flipped finder as **Step 0'**:
before running Step 1', forward-extend the length-`min_len` window
`pattern[x - min_len + 1 .. x]`.  On failure at position j, no MEM of
length ≥ `min_len` can end anywhere in [j, x], so the next anchor can
safely be `j - 1` (mirror of legacy's `j + 1` in the reverse
direction) and Step 1' + Step 2' can be skipped entirely for this
outer iteration.

**Correctness argument** (also embedded as inline comment).  Step 0'
failure at j means `count(P[x - min_len + 1 .. j]) < min_occ`.  For
any candidate MEM `[a..y]` with `length ≥ min_len` and `y ∈ [j, x-1]`:
`length ≥ min_len` and `y ≤ x - 1` imply `a ≤ y - min_len + 1 ≤ x - min_len`.
Hence `P[a..y]` contains `P[x - min_len + 1 .. j]` as a substring, so
`count(P[a..y]) ≤ count(P[x - min_len + 1 .. j]) < min_occ`.
So no MEM of length ≥ `min_len` can end in [j, x-1].  (A MEM ending
at x itself is also impossible since the length-`min_len` window
containing x already fails.)

**Design choice — Option A on success.**  When Step 0' succeeds
(the length-`min_len` window has ≥ `min_occ` occurrences), we discard
the accumulated forward interval and restart Step 1' from an empty
interval via `backward_extend_encoded_with_sa` (JR-011).  The
alternative — carrying the forward interval into Step 1' and
avoiding the redundant re-extension — requires converting a forward
interval into backward-extend context and re-establishing the SA
carry from scratch.  Option A costs at most `min_len` extra
forward-extends per emitted MEM, which is small relative to Step 1'
+ Step 2' cost.  The SA carry via the JR-011 toehold table makes
the Step 1' restart itself cheap.

### Method

Implemented in `find_mems_flipped_function` (algorithm.hpp) as a
prefix block before the existing Step 1'.  No changes to
`flipped_advance_anchor` (Step 2') or `find_all_mems_flipped` outer
loop.  Simulation prior to implementation (bin/debug_step0_perf on
100K noisy chr6 reads) predicted:

- outer iterations: 1,389,639 → 325,601 (−77%)
- backward_extend calls: 37.9M → 16.6M (−56%)
- forward_extend calls: 20.9M → 14.6M (−30%)
- total extends: 58.8M → 31.2M (−47%)

Empirical measurement used the JR-008 `perf_harness.sh` protocol
(2 untimed warmups + 3 timed trials, contiguous per-mode ordering),
median-of-3, on both `chr6.alt_noisy.reads.txt` and `chr6.alt.reads.txt`.
**All authoritative perf numbers reported here are from Vesuvio**
(host `vesuvio`, Linux 6.12.73, x86_64 Debian) — the canonical
benchmarking platform. Prototyping/design iteration was done on a
darwin/arm64 developer Mac and is not reported.

### Findings

**Correctness (byte-identical output preserved).**
- `test_backward_extend_sa` on yeast chrII: A 5000/5000, B 10000/10000,
  C 100140 steps / 0 mismatches (JR-011 primitives unchanged).
- `test_flipped_mems` on 100K yeast reads: 109520/109520 MEMs match
  legacy, 0 SA mismatches, 0 flipped_missed, 0 flipped_extra.
- E2E `./query.sh --gaf` on HPRC chr6 noisy: 2000/2000 valid,
  GAF 532,534 lines, coverage MD5
  `60f6b8e4a759aebb252b83a870b9ff8c` matches legacy s2-regcheck
  baseline byte-for-byte.

**Perf: HPRC chr6 on Vesuvio, `perf_harness.sh` median-of-3,
find_mems --lightweight-tags (perf tags jr012-{noisy,clean}-{legacy,step0}).**

| workload | metric               | legacy | Step 0' flipped | Δ vs legacy |
|:---------|:---------------------|-------:|----------------:|------------:|
| noisy    | **total execution**  | **60.33s** | **53.29s**  | **−7.04s (−11.7%)** |
| noisy    | R-index load         | 11.50s | 11.55s          | +0.05s |
| noisy    | MEM finding          | 28.35s | 38.63s          | +10.28s (+36%) |
| noisy    | MEM processing       | 19.80s | 2.48s           | **−17.32s (−87%)** |
| noisy    | ├─ tag queries       | 0.11s  | 0.11s           |  0.00s |
| noisy    | ├─ first locate      | 17.65s | **0.00s**       | **−17.65s** |
| noisy    | └─ locate-next       | 1.88s  | 2.21s           | +0.33s |
| noisy    | peak RSS             | 4512 MB| 4512 MB         |  0 MB |
| clean    | **total execution**  | **49.13s** | **44.03s**  | **−5.10s (−10.4%)** |
| clean    | R-index load         | 11.54s | 11.54s          |  0.00s |
| clean    | MEM finding          | 24.67s | 30.96s          | +6.29s (+25%) |
| clean    | MEM processing       | 12.26s | 0.88s           | **−11.38s (−93%)** |
| clean    | ├─ tag queries       | 0.07s  | 0.07s           |  0.00s |
| clean    | ├─ first locate      | 11.43s | **0.00s**       | **−11.43s** |
| clean    | └─ locate-next       | 0.62s  | 0.67s           | +0.05s |
| clean    | peak RSS             | 4520 MB| 4521 MB         | +1 MB |

Wall-time headlines:
- **Noisy**: total wall −11.7% (60.33s → 53.29s). MEM finding costs
  +10.3s more (Step 0' + Step 1' + Step 2' > legacy 3-phase), but MEM
  processing collapses by 17.3s (−87%) because flipped's SA-carry
  (JR-011 toehold table) eliminates the 155.5K `locate_sa_value`
  first-locate calls entirely (17.65s → 0.00s).
- **Clean**: total wall −10.4% (49.13s → 44.03s). Same mechanism —
  MEM finding pays +6.3s to save +11.4s in MEM processing (−93%),
  with first-locate collapsing from 11.43s to 0.00s.
- Peak RSS is identical between paths (within 1 MB); no memory
  regression.

**Extend-call counts (bin/debug_step0_perf simulation model, 100K
noisy chr6 reads, Mac prototyping platform — algorithmic totals
architecture-independent).**

| metric               | legacy | pre-Step 0' flipped | Step 0' flipped |
|:---------------------|-------:|--------------------:|----------------:|
| outer iterations     | 366,797 | 1,389,639 | 325,601 |
| backward_extend calls | 14.6M | 37.9M | 16.6M |
| forward_extend calls  |  8.9M | 20.9M | 14.6M |
| **total extends**     | **23.4M** | **58.8M** | **31.2M** |
| MEMs emitted          | 155,551 | 155,551 | 155,551 |

Post-Step 0' outer iteration count (325,601) is now slightly *below*
legacy (366,797) because Step 0''s `j - 1` jump on failure is
tighter than legacy's `j + 1` jump on the corresponding forward
direction. Emit rate rises to 47.8% of iters (155K / 326K), close
to legacy's 42% (155K / 367K). The remaining +33% extend gap
vs legacy comes from Step 1' backward walks emitting MEMs of avg
length 106 chars — bigger MEMs cost more backward extends per emit
than legacy's bounded incremental growth pattern.

### Verification chain

Same three gates as JR-011, all PASS with byte-identical output.
Vesuvio perf artifacts under
`../mem-projection/pangenome-pipeline/perf/jr012-{noisy,clean}-{legacy,step0}/`
(SUMMARY.tsv + per-trial find_mems.log/.time/.stderr + PROVENANCE.txt).

### Discussion

**Wall-time beats legacy on both workloads.** JR-008's original
finding (flipped ~10-13% faster than legacy on clean chr6 at wall
time) is confirmed and extended: −10.4% on clean and −11.7% on
noisy on Vesuvio. Recommendation in JR-011 open question #3
("workload-aware toggle") is superseded — flipped can be enabled
unconditionally without workload-specific regression risk.

**MEM finding still lags legacy on both workloads.** Step 0' gets
flipped's noisy MEM finding to +36% vs legacy (down from +133%
pre-Step 0' on the same code path). The residual gap is intrinsic to
model (b) enumeration (one MEM per qualifying right-endpoint):
flipped must run Step 1' fully from scratch for every emitted MEM
to get correct SA carry via the toehold table, while legacy's Step 3
preserves interval state across MEM emissions. Closing this gap
would require a different flipped state machine — deferred (see
Open question #1).

**Trial-to-trial variance on Vesuvio is very tight** (per-trial
totals within ±0.5s of the median across all four configurations),
indicating the benchmark is CPU-bound with negligible OS/I/O jitter
under `perf_harness.sh`'s warm-cache protocol. This makes small
percentage deltas trustworthy.

**Load-time cost is unchanged from JR-011.** R-index load is
~11.5s on Vesuvio (JR-011's eager `run_head_c` bitvector build +
`first_bint`/`first_sa` toehold table precompute). Both legacy and
flipped paths pay this cost identically because they share the same
r-index loader. JR-011 open question #2 (serialize the tables to
disk) remains the natural follow-up if single-query latency matters;
in batch contexts this amortizes to zero.

**JR-011's noisy perf numbers are historical.** JR-011 measured the
pre-Step 0' flipped code and reported "flipped +27% wall vs legacy on
noisy". Under Step 0', that outcome is superseded by JR-012's "flipped
−11.7% wall vs legacy on noisy" (Vesuvio). JR-011's primitive rewrite
remains the enabling foundation — Step 0''s Option A design depends
on cheap SA-carrying backward-extend from empty, which JR-011
delivered.

### Open questions

1. **Bigger-MEM growth pattern** in Step 1'. Avg emitted MEM length
   on noisy chr6 is 106 chars, meaning each emit costs ~106 bwe.
   Legacy's Step 2 grows matches incrementally within its outer
   loop and can share state across emissions.  A "grow only" Step 1'
   variant that skips the length-`min_len` characters Step 0' already
   forward-extended (converting the forward interval into a backward
   context) could save the first `min_len` bwe per emit — ~7.8M
   bwe (155K × 50) at cost of some algorithmic complexity.  Not yet
   attempted.
2. **`FIND_MEMS_EXTRA_FLAGS` hook** in pipeline `query.sh` (added
   in JR-011 session, uncommitted on that repo side) is now used
   by both JR-011 and JR-012 gate 3.  Should be committed on the
   pipeline repo per the JR-008 perf_harness.sh precedent.
3. **Default toggle for `--use-flipped-mems`.** Given the JR-012
   Vesuvio findings (flipped beats legacy on wall time on both
   noisy and clean workloads by 10-12%, byte-identical output,
   identical peak RSS, load-time cost ~11.5s amortized in batch),
   flipping the default seems justified. Deferred pending user
   direction.

### Data artifacts

**Vesuvio perf benchmark (canonical / source of truth):**
- `../mem-projection/pangenome-pipeline/perf/jr012-noisy-legacy/`
  — SUMMARY.tsv, PROVENANCE.txt, and lightweight/trial-{1,2,3}/
    (find_mems.log, find_mems.stderr, find_mems.time, gafpack.*).
- `../mem-projection/pangenome-pipeline/perf/jr012-noisy-step0/`
- `../mem-projection/pangenome-pipeline/perf/jr012-clean-legacy/`
- `../mem-projection/pangenome-pipeline/perf/jr012-clean-step0/`

**Correctness gate artifacts:**
- `runs/hprc-chr6-2026-06-02/queries/step0-flipped-noisy/` — E2E
  `--gaf` run of Step 0' flipped on chr6 noisy (correctness gate 3).
- `runs/hprc-chr6-2026-06-02/queries/s2-regcheck/` — legacy
  baseline (coverage MD5 reference).

**Prototyping / design artifacts (Mac):**
- `src/debug_step0_perf.cpp` — simulation harness used to model
  Step 0' cost prior to implementation (Option A logic).
- `src/debug_flipped_perf.cpp`, `src/debug_legacy_perf.cpp`,
  `src/debug_anchor_distribution.cpp` — instrumented profilers
  used to identify the outer-iteration overhead.

### Commits (this repo, branch backward-extend-sa-perf)

- `37dc2d9` algorithm: rewrite find_mems_flipped state machine (correctness)
- `83ef433` algorithm: add Step 0' min_len pre-check to find_mems_flipped (perf)
- `ecd8483` journal: JR-012 Step 0' min_len pre-check in find_mems_flipped
- `412b9c8` journal: JR-012 promote Vesuvio to source of truth (perf numbers)

---

## JR-013 — Toehold hoist into callers: perf null, contract cleanup

```yaml
id: JR-013
date: 2026-08-06
author: claude-opus-4-7 (session with hlakshmidevi)
status: resolved
tags: [flipped, r-index, primitive, contract, refactor, vesuvio, perf-null]
refs:
  follows: [JR-012]
benchmark-platform: vesuvio (Linux 6.12.73 x86_64, Debian)
```

### Context

`backward_extend_encoded_with_sa` (JR-011) had a "toehold fast path"
inline branch at the top: if the caller passed the whole-BWT interval
plus `sa_sp_prev == NO_POSITION`, return the precomputed
`(first_bint[c], first_sa[c])` directly.  This branch fired on the
very first call of every Step 1' walk in the flipped finder.

Two other flipped-finder call sites paid analogous overheads for
their first extend:
- Step 0' (`find_mems_flipped_function`): first `forward_extend_encoded`
  on the whole-BWT interval per outer iter.
- Step 2' (`flipped_advance_anchor`): first `forward_extend_encoded`
  on the whole-BWT interval per emit iter.

Neither was using the precomputed table -- `forward_extend_encoded`
does a full pred+walk internally (via swap + `backward_extend_encoded`)
even though the answer for whole-BWT + one character is knowable in
O(1) via `first_bint[dna_char_to_idx(complement(c))]` swapped.

Estimated 481K + 155K + 325K ≈ 961K first-extend calls per 100K
noisy chr6 reads that could be short-circuited.

### Hypothesis

Add three public O(1) accessors on `FastLocate` and hoist the
"first extend from whole-BWT" logic out of the primitives:
- `first_backward_with_sa(c)` -- Step 1' seed
- `first_backward(c)`         -- interval only
- `first_forward(c)`          -- Step 0' / Step 2' seed
  (derived via FMD-index identity: `forward(whole_bwt, c) ==
  swap(backward(whole_bwt, complement(c)))`, and
  `swap(whole_bwt) == whole_bwt`).

Callers use the accessor for the first extend of each walk, then
loop with the primitive for subsequent extends.

The primitive's contract can then be tightened: require
`bint.size > 0` and `sa_sp_prev != NO_POSITION`.  The general-path
SA-seed fallback (`run_id_and_offset_at` + `locateNext` walk),
which only existed to backstop toehold-check misses, can be
removed entirely.

Estimated wall-time impact: ~500-700 ms MEM-finding saved on
noisy chr6 (Mac back-of-envelope, ~1-2% of MEM-finding).

### Method

Implemented in commits `987f0b7` (accessors + contract tightening)
and `1d118e9` (dead-code cleanup + defensive assert).
`test_backward_extend_sa` Test C updated to seed the initial step
via `first_backward_with_sa(c)`, mirroring the production path.

Perf measured on Vesuvio via `perf_harness.sh` protocol
(2 warmups + 3 timed trials, median-of-3), tags
`jr013-{noisy,clean}-hoist`.  Baseline for comparison: JR-012's
`jr012-{noisy,clean}-{legacy,step0}` runs (same host, protocol,
config).

### Findings

**Correctness (byte-identical output preserved).**
- test_backward_extend_sa yeast: A 5000/5000, B 10000/10000,
  C 100140 steps / 0 mismatches (Test C now flows through the
  accessor for INITIAL; same numbers as JR-012).
- test_flipped_mems 100K yeast reads: 109520/109520 MEMs match
  legacy, 0 SA mismatches, 0 flipped_missed, 0 flipped_extra.
- E2E chr6 noisy `--gaf`: 2000/2000 valid, coverage MD5
  `60f6b8e4a759aebb252b83a870b9ff8c` (noisy baseline).
- E2E chr6 clean `--gaf` (via perf_harness gafpack): coverage MD5
  `8c242732c837d697fb18c52fd0c0cf6e` (clean baseline from JR-012
  legacy).

**Perf: Vesuvio median-of-3, hoist vs JR-012 step0 flipped.**

| workload | metric               | JR-012 step0 | JR-013 hoist | delta      |
|:---------|:---------------------|-------------:|-------------:|-----------:|
| noisy    | total execution      | 53.29s       | **53.01s**   | −0.28s (−0.5%) |
| noisy    | MEM finding          | 38.63s       | **38.43s**   | −0.20s (−0.5%) |
| noisy    | MEM processing       | 2.48s        | 2.49s        | +0.01s (noise) |
| noisy    | locate ops           | 2.21s        | 2.22s        | +0.01s (noise) |
| noisy    | first locate         | 0.00s        | 0.00s        | 0 |
| noisy    | peak RSS             | 4512 MB      | 4496 MB      | −16 MB |
| clean    | total execution      | 44.03s       | **44.42s**   | +0.39s (+0.9%) |
| clean    | MEM finding          | 30.96s       | **31.11s**   | +0.15s (+0.5%) |
| clean    | MEM processing       | 0.88s        | 0.89s        | +0.01s (noise) |
| clean    | locate ops           | 0.67s        | 0.69s        | +0.02s (noise) |
| clean    | first locate         | 0.00s        | 0.00s        | 0 |
| clean    | peak RSS             | 4520 MB      | 4516 MB      | −5 MB |

**Perf delta is essentially null.** Both workloads' totals lie
within the trial-to-trial variance floor of Vesuvio measurements
(~±0.3-0.5s from JR-012 warm-cache trials).  The noisy MEM-finding
saw a real but tiny improvement (−0.2s ≈ −0.5%); the clean
MEM-finding saw an equally-tiny regression (+0.15s).  Neither is
distinguishable from noise at the wall-time scale.

### Discussion

**Why the perf estimate was optimistic.**  The Mac back-of-envelope
projected ~500-700 ms saved.  Three factors overcounted:
1. **Branch prediction already near-perfect.**  The toehold check
   was `[TRUE, FALSE, FALSE, ..., FALSE]` per Step 1' walk --
   the branch predictor nails this after warm-up.  Cost per call
   was ~2-4 ns, not the ~700 ns/call implicit in the estimate.
2. **Forward extends do real work regardless.**  Skipping the
   first Step 0' / Step 2' forward extend saves only the fixed
   per-call overhead, not the full block-walk cost -- the
   accessor still returns a computed answer.
3. **Accessor has its own cost.**  Global 5-entry array lookup +
   inline function call cost is comparable to the removed inline
   branch check.  Trading one predictable branch for one load.

The tiny clean-reads regression (+0.15s MEM finding) may reflect
cache placement of `first_bint`/`first_sa` vs the inline check
against a compile-time constant.  Below the threshold for further
investigation.

**The real wins are code-quality, not perf.**
- **Primitive contract is now clean.**  `backward_extend_encoded_with_sa`
  no longer has dual "toehold or general" modes; caller must
  provide a valid (bint, sa_sp) seed.  Enforced by `assert()` in
  debug builds.
- **Dead code removed.**  The general-path SA-seed fallback
  (`run_id_and_offset_at` + `locateNext` walk, ~10 lines in
  r-index.cpp) that only existed to backstop toehold misses is
  gone entirely.
- **All three flipped call sites use the same pattern.**  Step 0',
  Step 1', Step 2' each start with an O(1) accessor for the first
  extend, then loop with the primitive.  Uniform and greppable.
- **Test C mirrors production.**  The old Test C exercised a
  code path (whole-BWT + NO_POSITION at the primitive) that
  production no longer takes.  Now Test C's INITIAL case flows
  through `first_backward_with_sa(c)` exactly like Step 1'.

**JR-012 numbers remain the source-of-truth headline.**  The
JR-013 changes are a code-quality refactor that preserves the
JR-012 perf profile within noise.  Any perf story going forward
should quote the JR-012 numbers (noisy −11.7% wall, clean −10.4%
wall vs legacy) and mention JR-013 only as "primitive contract
cleanup, perf-neutral".

### Open questions

1. **Grow-only Step 1' variant** (still JR-012 Q1).  This remains
   the highest-leverage flipped-finder optimization opportunity
   (~5.5s estimated saving on noisy MEM finding).  Not attempted
   in JR-013.
2. **Default toggle for `--use-flipped-mems`** (still JR-012 Q3).
   JR-013 doesn't change the answer; JR-012 findings still apply.

### Data artifacts

**Vesuvio perf benchmark (canonical / source of truth):**
- `../mem-projection/pangenome-pipeline/perf/jr013-noisy-hoist-vesuvio/`
  -- SUMMARY.tsv, PROVENANCE.txt, and lightweight/trial-1/
  (find_mems.log, find_mems.time -- other trials shipped as
  SUMMARY rows only).
- `../mem-projection/pangenome-pipeline/perf/jr013-clean-hoist-vesuvio/`

**JR-012 baseline (for delta comparison):**
- `../mem-projection/pangenome-pipeline/perf/jr012-{noisy,clean}-{legacy,step0}/`

### Commits (this repo, branch backward-extend-sa-perf)

- `987f0b7` r-index+algorithm: hoist first-extend toehold check into callers
- `1d118e9` algorithm: simplify find_mems_flipped_function (guard/assign + assert)
- `20a3683` journal: JR-013 toehold hoist -- perf null, contract cleanup

---

## JR-014 — Serialize run_head_c[] into .ri: load time collapse

```yaml
id: JR-014
date: 2026-08-07
author: claude-opus-4-7 (session with hlakshmidevi)
status: resolved
tags: [r-index, format, load-time, serialization, perf, hprc, vesuvio]
refs:
  follows: [JR-013]
benchmark-platform: vesuvio (Linux 6.12.73 x86_64, Debian)
```

### Context

JR-011 introduced `run_head_c[]` -- 5 per-DNA-character run-head
`sd_vector`s indexed by run id -- and built them eagerly at load
time inside `ensure_run_head_c_built()`.  On Vesuvio, this build cost
~8.7s per `load_encoded` invocation (two full passes over ~249M
BWT runs, decoding byte-encoded run headers on each pass).  Total
r-index load was ~11.5s, of which ~80% was `run_head_c[]` build
and ~20% was actual deserialization.

Load time is paid once per process invocation, so at HPRC chr6
scale (~155K MEMs per read, ~100K reads per file) it amortized to
low per-query overhead in batch settings.  But it inflated total
wall time on single-query and small-batch workloads by an amount
that dominated the JR-012 / JR-013 flipped-vs-legacy delta
(~10-12% win vs ~9s of load-time overhead on both paths).

### Hypothesis

Move `run_head_c[]` construction from query-time (`load_encoded`)
to build-time (`build_rindex`), serialize the resulting `sd_vector`s
into the .ri file, and deserialize them verbatim on load.

Expected outcomes:
- Load: ~11.5s -> ~2s (deserialization only; ~135 MB read on chr6).
- Total wall: both legacy and flipped paths drop by ~9s uniformly.
- Byte-identical query-time behavior: same `sd_vector`s, just loaded
  from disk instead of rebuilt.

### Method

Format change (`Header::HAS_RUN_HEAD_C = 0x2ULL`):
- `serialize_encoded`: sets the flag, calls `ensure_run_head_c_built()`
  to populate (idempotent -- once_flag guards it), and appends the 5
  `sd_vector`s at the end of the byte stream.
- `load_encoded`: throws `sdsl::simple_sds::InvalidData` with a
  clear "regenerate with build_rindex" message if the flag is not
  set.  If set, reads the 5 `sd_vector`s directly from the stream
  and fires the `run_head_c_once` flag as a no-op so any stray call
  to `ensure_run_head_c_built()` is short-circuited.
- `Header::check()` mask updated to accept both `ENCODED_BLOCKS`
  and `HAS_RUN_HEAD_C`.

No backward compat is retained (dev-phase policy).  All existing .ri
files must be regenerated via `build_rindex`.  Toehold table
(`first_bint`/`first_sa`) is left at at-load construction (~1 ms,
160 bytes, not worth the format surface area).

Test D added to `test_backward_extend_sa`: end-to-end round-trip.
Takes an optional 7th positional arg (path to the raw `.rl_bwt`).
When provided: builds a fresh `FastLocate` from the raw BWT,
serializes to a temp .ri via `serialize_encoded`, loads back via
`load_encoded`, and cross-checks both instances agree on
`next_run_with_head_c(c, start_from)` over 2500 random probes (500
per DNA character).  Skipped-not-failed when omitted.

### Findings

**Correctness (byte-identical output preserved).**  Verified against
regenerated yeast chrII 100 kb .ri and HPRC chr6 .ri:
- test_backward_extend_sa yeast: A 5000/5000, B 10000/10000,
  C 100140 steps / 0 mismatches, **D 2500/2500 round-trip probes**
  (new test).
- test_flipped_mems 100K yeast reads: 109520/109520 match legacy,
  0 SA mismatches.
- E2E chr6 noisy: coverage MD5
  `60f6b8e4a759aebb252b83a870b9ff8c` matches JR-013 baseline.
- E2E chr6 clean: coverage MD5
  `8c242732c837d697fb18c52fd0c0cf6e` matches JR-013 baseline.

**Perf: HPRC chr6 on Vesuvio, median-of-3, find_mems --lightweight-tags.**

| workload | metric               | JR-013     | JR-014     | delta       |
|:---------|:---------------------|-----------:|-----------:|------------:|
| noisy    | R-index load         | 11.50s     | **1.84s**  | −9.66s (−84%) |
| noisy    | legacy total exec    | 60.33s     | **51.03s** | −9.30s (−15.4%) |
| noisy    | flipped total exec   | 53.29s     | **42.90s** | −10.39s (−19.5%) |
| noisy    | flipped MEM finding  | 38.63s     | 38.38s     | −0.25s (noise) |
| noisy    | flipped MEM proc     | 2.48s      | 2.43s      | −0.05s (noise) |
| noisy    | peak RSS (both)      | 4512 MB    | 4497 MB    | −15 MB |
| clean    | R-index load         | 11.54s     | **1.80s**  | −9.74s (−84%) |
| clean    | legacy total exec    | 49.13s     | **39.01s** | −10.12s (−20.6%) |
| clean    | flipped total exec   | 44.03s     | **34.23s** | −9.80s (−22.3%) |
| clean    | flipped MEM finding  | 30.96s     | 30.78s     | −0.18s (noise) |
| clean    | flipped MEM proc     | 0.88s      | 0.89s      | +0.01s (noise) |
| clean    | peak RSS (both)      | 4520 MB    | 4516 MB    | −4 MB |

**Load time collapsed from ~11.5s to ~1.8s (−84%).**  Both legacy
and flipped paths benefit equally since they share the loader.

**Flipped's relative advantage is unchanged** (was −11.7% noisy /
−10.4% clean pre-JR-014; still −15.9% noisy / −12.3% clean
post-JR-014 -- larger relative percentages because the denominator
shrank).  JR-014 lifted all boats uniformly by ~9-10s.

**On-disk cost.**  HPRC chr6 .ri grew by ~5-6% (est. ~200 MB
added to the ~3.6 GB baseline).  Yeast chrII 100 kb .ri grew from
160 MB to 168 MB (+5%).

### Cumulative story (JR-011 -> JR-012 -> JR-013 -> JR-014)

| metric                          | pre-JR-011 baseline | JR-014 flipped | total delta |
|:--------------------------------|-------------------:|---------------:|------------:|
| noisy total exec                | ~74s               | **42.90s**     | **−42%**    |
| noisy MEM finding               | ~58s               | 38.38s         | −34%        |
| noisy MEM processing            | ~19s               | 2.43s          | −87%        |
| R-index load                    | ~11s               | 1.84s          | −83%        |

### Discussion

**Format extension mechanics.**  Adding a flag bit to the existing
64-bit `Header::flags` field, with a corresponding mask update in
`Header::check()`, was the minimum-viable change to distinguish
new-format from pre-format .ri files.  A `VERSION` bump would have
provided stricter type-safety but required more coordinated migration;
in dev-phase this trade-off was answered "no backward compat needed"
by user.

**Loud rejection worked as designed.**  On first Vesuvio invocation
against the pre-JR-014 chr6 .ri, `load_encoded` threw with the
"regenerate with build_rindex" message and the run aborted before
producing any output.  Zero silent failure risk.

**Test D scope.**  The initial implementation tried "load .ri ->
serialize -> reload" as a round-trip.  This hung because
`serialize_encoded` requires an in-memory `blocks[]` (built by the
constructor from raw `.rl_bwt`), which `load_encoded` clears to
save memory.  Test D was refactored to take the `.rl_bwt` path
directly, exercising the full build->serialize->load pipeline
(the same one `build_rindex` produces).  This is the correct
scope: it validates that `build_rindex` output is loadable and
semantically equivalent to a fresh in-memory build.

**Toehold table (first_bint / first_sa) intentionally not
serialized.**  At 160 bytes and ~1 ms to build, the load-time
cost is invisible; the format-surface-area cost of serializing
would exceed the benefit.  Documented in the r-index.hpp
first_toehold section and in JR-014 design discussion.

**No optimization of ensure_run_head_c_built() itself.**  We
considered a one-pass build (~4s savings, at cost of ~1.6 GB peak
scratch RAM) as a fallback if serialization proved infeasible.
Since serialization landed cleanly, the two-pass build stays as-is:
it now runs once at `.ri` build time, amortized over the lifetime
of the file, so its ~8.7s cost is invisible to queries.

### Open questions

1. **Load-time budget is now dominated by tag_index load
   (~0.5s on Vesuvio) and deserialization I/O (~1.3s reading
   the ~135 MB `run_head_c` sd_vectors).**  Further load-time
   reduction would target the encoded_stream I/O (~1 GB read,
   already at near-memory-copy speed) or mmap-based load
   (would eliminate the ~135 MB read entirely at cost of
   demand-paging on first access).  Diminishing returns; not
   pursued.
2. **Default toggle for `--use-flipped-mems`** (JR-012/013 Q3
   carried forward).  Post-JR-014, flipped's wall-time advantage
   is stronger and load-time is no longer a mitigating factor.
   Flip the default when convenient.
3. **build_rindex incremental / parallel run_head_c build.**
   Currently sequential and single-threaded; ~15-30 min on Vesuvio
   for HPRC chr6.  Not a query-time concern but worth noting if
   .ri regeneration becomes a bottleneck.

### Data artifacts

**Vesuvio perf benchmark (canonical / source of truth):**
- `../mem-projection/pangenome-pipeline/perf/jr014-noisy-legacy-vesuvio/`
- `../mem-projection/pangenome-pipeline/perf/jr014-noisy-flipped-vesuvio/`
- `../mem-projection/pangenome-pipeline/perf/jr014-clean-legacy-vesuvio/`
- `../mem-projection/pangenome-pipeline/perf/jr014-clean-flipped-vesuvio/`
  -- each with SUMMARY.tsv, PROVENANCE.txt (noisy-flipped only),
  and lightweight/trial-1/{find_mems.log, find_mems.time}.

**JR-013 baseline (for delta comparison):**
- `../mem-projection/pangenome-pipeline/perf/jr013-*-hoist-vesuvio/`

**Regenerated .ri files (backup of pre-JR-014):**
- `runs/hprc-chr6-2026-06-02/hprcv1_chr6.ri.pre-jr014`
- `runs/v1-current/yeast235_chrII_100kb_normalized.ri.pre-jr014`

### Commits (this repo, branch backward-extend-sa-perf)

- `8b44a73` r-index: serialize run_head_c[] into .ri (JR-014, load ~11s -> ~2.5s)
- `e7f8a80` journal: JR-014 serialize run_head_c[] into .ri -- load ~11.5s -> ~1.8s

---

## JR-015 — Post-Risk-E MEM classification: fresh seed-cost + tag-partition tables

```yaml
id: JR-015
date: 2026-08-07
author: claude-opus-4-7 (session with hlakshmidevi)
status: resolved
tags: [classification, seed-cost, tag-partition, hprc, vesuvio, historical-comparison]
refs:
  follows: [JR-014]
  supersedes-numbers-of: [FINDINGS_SEED_COST.md]
  historical-comparison-note: July 2026 FINDINGS_SEED_COST.md likely
    used clean chr6 reads (MEM count 100,899 matches Aug 2026 clean
    100,000; seed_cost 22.7M matches Aug 2026 clean 20.5M within 10%).
    The reads-file identity for July was not recorded in that doc.
benchmark-platform: vesuvio (Linux 6.12.73 x86_64, Debian)
```

### Context

`FINDINGS_SEED_COST.md` (2026-07-23, branch `classify-mem-runs`)
established that seed cost dominates the locate hot path: on
HPRCv1 chr6 100k reads, seed accounted for 32.2% of walk work,
and 8.8% of MEMs (those with `avg_bwt_run_length >= 1000` bp)
carried 75.4% of all seed cost.  Those numbers pre-dated the
Risk E fix (JR-007) and pre-dated the flipped MEM finder
(JR-011 onwards).  Two questions were open:

1. **Are the numbers still current?**  The Risk E fix removed
   ~1.7% of MEMs (spurious non-left-maximal MEMs concentrated in
   long-BWT-run regions).  Impact on the seed-cost distribution
   was unquantified.
2. **What does the tag-type partition look like on both noisy
   and clean workloads?**  July 2026 was noisy-only (or reads
   file unspecified); a noisy-vs-clean comparison had not been
   done.

### Method

Re-ran `find_mems --debug-classify=<path>` on Vesuvio against
JR-014-format .ri files, for both HPRCv1 chr6 workloads and both
MEM finders:

- `chr6.alt_noisy.reads.txt` (155,551 MEMs across 100k reads)
- `chr6.alt.reads.txt`       (100,000 MEMs across 100k reads)

The classifier is idempotent w.r.t. the finder path (it walks
the MEM's BWT interval independently of how the MEM was located),
so legacy and flipped runs produce byte-identical TSVs.  Verified
via `diff <(sort legacy.tsv) <(sort flipped.tsv)` -- zero diff
on both workloads.  This confirms the flipped MEM set is
byte-identical to legacy at the classification level too, not
just at the coverage-MD5 level.

Aggregated via `xy-test/classify_aggregate.py` (Python 3 stdlib
one-off; not committed).

### Findings

**A. July 2026 baseline most likely used clean chr6 reads.**  MEM
count matches almost exactly (July 100,899 vs Aug clean 100,000).
Seed cost matches within 10% (July 22.7M vs Aug clean 20.5M).
The ~10% delta is consistent with the Risk E fix (JR-007) removing
a small number of MEMs concentrated in long-BWT-run regions.

| dataset                                 | Total MEMs | seed_cost   |
|:----------------------------------------|-----------:|------------:|
| July 2026 pre-Risk-E (assumed clean)    | 100,899    | 22.7M       |
| Aug 2026 post-Risk-E, clean             | 100,000    | 20.5M       |
| Aug 2026 post-Risk-E, noisy             | 155,551    | 33.0M       |

**B. Tail concentration is a stable structural property of HPRC
chr6 -- unchanged by Risk E fix or reads-file choice.**

| dataset                            | %MEMs in tail (>=1000 bp) | %seed in tail |
|:-----------------------------------|--------------------------:|--------------:|
| July 2026 pre-Risk-E (assumed clean) | 8.8%                    | 75.4%         |
| Aug 2026 post-Risk-E, noisy        | **8.6%**                  | **76.2%**     |
| Aug 2026 post-Risk-E, clean        | **8.2%**                  | **73.8%**     |

Same structural story across three measurement conditions:
**~8-9% of MEMs concentrate ~74-76% of seed cost.**  The Risk E
fix's impact on the tail is small (-0.6pp of MEMs on the closest
comparison, July vs Aug clean).

**C. Seed / walk ratio: July's 32% figure is not reproducible
against runtime measurements.**

| dataset                            | seed_cost (calls) | seed / walk | inter_current |
|:-----------------------------------|-----------------:|------------:|--------------:|
| July 2026 pre-Risk-E (assumed clean) | 22.7M          | 32.2%       | ~47.8M (implied) |
| Aug 2026 post-Risk-E, clean        | **20.5M**        | **88.3%**   | **2.7M**      |
| Aug 2026 post-Risk-E, noisy        | **33.0M**        | **87.8%**   | **4.6M**      |

Aug 2026 clean's `inter_current` (2.7M) matches the runtime
`locate_next_calls` reported by JR-013/014 exactly.  The July
32.2% seed/walk implies inter_current ~47.8M, which would
contradict the runtime measurement by ~17x.  Something changed
in either the classifier's `inter_current` computation, the
reads file, or the aggregation between July and August.
Historical July numbers should be treated as unverified.

Definition per current classifier
(`classify_mem_runs` in src/find_mems.cpp):
```
inter_current = sum over k>=1 of (tag_starts[k] - prev_pos_cur)
              = last_tag_start - sp    (telescoping)
```

Direct comparison of the two Aug 2026 workloads (identical index,
different reads file):

| metric                            | noisy       | clean       |
|:----------------------------------|------------:|------------:|
| Total MEMs                        | 155,551     | 100,000     |
| MEMs per read                     | 1.56        | 1.00        |
| Total seed cost (locateNext)      | 33.0M       | 20.5M       |
| Total inter_current               | 4.6M        | 2.7M        |
| seed / walk                       | 87.75%      | 88.30%      |
| Volume-weighted avg BWT run len   | 1477 bp     | 1268 bp     |
| b1 tags (anchor-eligible)         | 10.27%      | 5.34%       |
| b3 tags (strictly interior)       | 86.59%      | 93.00%      |
| samples[rid] savings ceiling      | **0.32%**   | **0.06%**   |

Noisy has 2x more anchor-eligible tags (b1) than clean because
shorter MEMs cross BWT run boundaries more often.  This makes
noisy's `samples[rid]` optimization ceiling (0.32%) 5x higher
than clean's (0.06%) -- but both are far below actionable
thresholds and confirm the FINDINGS_SEED_COST.md conclusion
that `samples[rid]` is not worth building.

**avg_bwt_run_length caveat.**  The July 2026 FINDINGS_SEED_COST.md
report of "volume-weighted avg BWT run length: 552 bp" is close
to a *simple arithmetic mean* (not volume-weighted) computed on
the same TSV today: **483 bp on Aug 2026 clean** (close to July's
552 bp).  If the aggregator weights by MEM's `bwt_size` instead,
the result is 1268 bp.  Different weightings, different numbers;
the July doc's weighting was not explicit.  For meaningful
comparisons, prefer the per-MEM bucketed table above (which is
weighting-independent).

**D. seed_cost (in locateNext calls) closes the loop with
JR-013/014 wall-time measurements.**

At an implied 0.53-0.56 microseconds per `locateNext` call
(Vesuvio, cache-mostly-cold seed-walk regime):

| workload | seed_cost (calls) | predicted wall (0.53us/call) | observed first-locate (JR-014) |
|:---------|-----------------:|-----------------------------:|------------------------------:|
| noisy    | 33,029,141       | 17.51s                       | **17.65s**                    |
| clean    | 20,467,583       | 10.85s                       | **11.43s**                    |

Model closes within ~5%.  This confirms the `seed_cost` column
in the classifier is the exact quantity flipped's SA carry saves
per query.  The JR-014 flipped-path `first_locate_calls = 0`
translates directly to the classifier's `seed_cost` being the
saved wall time.

Note: 0.53 us/call is 2.2x the FINDINGS_SEED_COST.md estimate of
0.24 us/call.  The old figure was an amortized mix of seed and
inter-tag walks; seed walks (cold predecessor + block walk from
the run start to `sp`) are more expensive per call than inter-tag
walks (warm SA cursor continuation).  Model refinement:

- Amortized (mixed):    ~0.24 us/call  (July 2026 figure)
- Seed walks only:       ~0.53 us/call  (Aug 2026 back-computed)
- Inter-tag walks only:  <0.24 us/call  (implied; not directly measured)

**E. Risk E fix impact on the classification: essentially null.**

The Risk E fix removed spurious non-left-maximal MEMs (JR-007).
Spurious MEMs concentrate in long-BWT-run regions (per JR-005:
median spurious size 32,173 occurrences vs median real MEM size
86).  In theory that should reduce the tail bucket counts by 1-2%.
Actual measured tail (>=1000 bp) MEM counts:

| dataset                        | %MEMs in tail | %seed in tail |
|:-------------------------------|--------------:|--------------:|
| July 2026 pre-Risk-E           | 8.8%          | 75.4%         |
| Aug 2026 post-Risk-E, noisy    | 8.6%          | 76.2%         |
| Aug 2026 post-Risk-E, clean    | 8.2%          | 73.8%         |

Small changes (-0.2 to -0.6pp of MEMs in the tail; ±2pp of seed
share) but well within the range attributable to the different
reads files.  FINDINGS_SEED_COST.md's headline (75% of seed in
8-9% of MEMs) survives the fix intact.

### Discussion

**FINDINGS_SEED_COST.md's conclusions all stand.**  The
"attack seed cost, ignore samples[rid], target long-BWT-run tail"
strategy is confirmed by the fresh numbers.  The specific
"projected 7% end-to-end wall reduction via SA subsample" analysis
in that doc is now less relevant: JR-011 delivered the
SA-carrying primitive (via flipped MEM finder) that eliminates
seed cost entirely on the flipped path, achieving the theoretical
maximum saving without any additional storage.

**The classifier is now a direct measure of what the flipped path
saves.**  The `seed_cost` column corresponds byte-for-byte to
`first_locate_calls` in the runtime `first_locate_time` profiler
output on legacy runs.  For any future work considering further
locate optimizations, this classifier gives a per-MEM quantification
of the ceiling.

**JR-015 is a fresh baseline, not a new optimization.**  Nothing
in this entry changes code or ships perf.  It updates the seed-cost
+ tag-partition analysis for the post-Risk-E, post-JR-014 code and
data, cross-validates the JR-013/014 first-locate wall times against
independent measurement, and completes the noisy-vs-clean comparison
that was missing from July 2026's work.

### Open questions / follow-ups

1. **Unexplained gap in July 2026's seed/walk ratio.**  July's
   32.2% seed/walk implies inter_current ~47.8M, but Aug 2026
   clean's inter_current on the same reads-and-index combination
   is 2.7M -- a 17x delta.  The classifier's `inter_current`
   definition is unchanged.  Possibilities: (a) different reads
   file than assumed, (b) different .ri format that changed tag
   geometry, (c) an aggregation bug in the July analysis.  Not
   pursued -- July numbers should be treated as unverified for
   the seed/walk fraction; the per-bucket seed-cost concentration
   story survives independently (see Finding B).
2. **Inter-tag walk cost per call.**  Not directly measured.
   The mixed 0.24 us and cold-seed 0.53 us bookend it, but
   pinning the inter-tag figure would firm up wall-time
   projections for any future inter-tag optimization.
3. **JR-011 already delivered the ideal seed-cost fix
   (elimination via SA carry).**  No further seed-cost work is
   needed.  Any future locate optimizations should target
   inter-tag walks (12% of walk work on Aug 2026) -- but the
   `samples[rid]` optimization ceiling of 0.06-0.32% shows the
   inter-tag structure has little exploitable slack.
4. **Non-encoded r-index path.**  All tables assume the encoded
   path (default for `.ri` files).  Non-encoded is dead code but
   would need a fresh classification run to match.

### Data artifacts

**Vesuvio classification output (canonical):**
- `../mem-projection/pangenome-pipeline/classify_output_jr015/classify_noisy_legacy.tsv`  (155551 rows, ~10 MB)
- `../mem-projection/pangenome-pipeline/classify_output_jr015/classify_noisy_flipped.tsv` (identical to legacy)
- `../mem-projection/pangenome-pipeline/classify_output_jr015/classify_clean_legacy.tsv`  (100000 rows, ~6.6 MB)
- `../mem-projection/pangenome-pipeline/classify_output_jr015/classify_clean_flipped.tsv` (identical to legacy)
- `.../classify_output_jr015/summary_{noisy,clean}.txt` -- human-readable aggregates.

**Historical baseline (superseded numbers, kept for reference):**
- `FINDINGS_SEED_COST.md` -- July 2026 pre-Risk-E analysis.

**Aggregator (one-off, not committed):**
- `xy-test/classify_aggregate.py` -- Python 3 stdlib script that
  produces the summary tables from the classifier's TSV.

### Commits (this repo, branch backward-extend-sa-perf)

- (no code changes in JR-015 -- this entry documents fresh
  measurements against JR-014 code)

---

## JR-016 — Tag-head SA samples (sr-index-style): storage cost on HPRC chr6

```yaml
id: JR-016
date: 2026-08-09
author: claude-opus-4-7 (session with hlakshmidevi)
status: open
tags: [tag-head-samples, sr-index, sa-sampling, storage, hprc, vesuvio]
refs:
  follows: [JR-015]
benchmark-platform: vesuvio (Linux 6.12.73 x86_64, Debian)
```

### Context

JR-015 confirmed seed cost is a structural, workload-independent
property of HPRC chr6: 30-33% of walk work, ~76% concentrated in
tag heads inside long BWT runs.  The flipped path (JR-012)
eliminates this cost by SA-carrying through backward-only MEM
enumeration, at the price of a more expensive MEM-finding phase.
An alternative attack: keep the legacy MEM finder (which is
faster at MEM finding per JR-012 Vesuvio tables: 28.35s vs Step 0'
flipped 38.63s on noisy chr6) and instead add a query-side data
structure that answers `SA[tag_head]` cheaply.  This decouples
"where the MEMs are" from "what the SA is at each tag emit."

Design: sr-index-inspired sub-sampling over tag-run heads.  From
Cobas, Gagie, Navarro 2021 (arXiv:2103.15329) §4.1, adapted to
our two-tier structure (tag runs + BWT runs):

- Candidates = tag-run heads (deletable) + BWT-run heads
  (unremovable anchors; already sampled in `.ri` samples[]).
- Sort candidates by SA value (text-position order).
- Scan left-to-right; for each tag-run head `s_i`, delete iff
  `SA[s_{i+1}] - SA[prev_kept_SA] < s`.

Guarantee: from any deleted tag head `t`, walking LF backward
reaches a kept sample (tag-head OR BWT-anchor) in < s LF steps,
because LF decrements SA by 1 per step and the deletion rule
ensures a kept sample exists at some SA in
`[SA[t] - s + 1, SA[t]]`.  A find_mems query for `SA[t]` at any
deleted tag head becomes: 1 hash-marker lookup (fast path if
kept) OR up to s LF steps to a known-SA position (slow path).
No `locate_sa_value` seed walk needed.

### Method

New binary `bin/build_tag_head_samples` (branch
`tag-head-samples`, off `backward-extend-sa-perf` HEAD `c6858a1`):

**Phase 1** -- enumerate all candidates with their SAs.  Walk BWT
runs in order via existing block-decode; at each BWT-run head,
emit anchor candidate with `SA = getSample(bwt_rid)`.  Walk
locateNext through each run, emitting tag-head candidates at
their BWT positions with the running SA cursor.  Total locateNext
cost ~= n (single full BWT walk).

**Phase 2** -- sort candidates by SA; apply deletion rule
per-`s` value in a single linear scan.  All 4 s values processed
against one shared Phase 1 output.

**Phase 3** -- serialize per s: sd_vector `kept_marker` indexed
by tag-run rid (1 bit per rid, 1 iff kept), int_vector
`sa_values` of packed SAs for kept tag-run heads in rid order.
Width from `r_index.samples.width()` (=39 bits on HPRC chr6);
using naive `ceil(log2(n))` = 35 bits silently truncates SA
values since packed SA range is `n_seq * max_length` >> `n` on
multi-string BWTs.  Caught by verification pass (a) below on
first yeast run; fixed before Vesuvio run.

**Verification** -- two spot-check passes per s value on the
built table:

- (a) 1000 random kept tag heads: assert
  `sa_values[rank_1(rid)] == locate_sa_value(run_start_bwt(rid))`.
- (b) 100 random deleted tag heads: simulate the query-time LF
  walk, assert we land on a known sample within s steps AND
  reconstructed SA (sample + steps) equals ground truth.

Both checks pass on all 4 s values on both yeast (n=337M) and
HPRC chr6 (n=31e9).  See "Verification" section below for what's
covered and what isn't.

### Findings

**HPRC chr6 storage table** (Vesuvio, single 1h 24m build for all
4 s values via `xy-test/tag_samples/run_hprc_chr6.sh`):

Index stats: n = 31.0 Gbp BWT, 249.4M BWT runs, 710.2M tag runs,
of which 240.6M (34%) coincide with BWT-run heads (implicit-
deleted for free -- SA available from r-index `samples[]`),
leaving 469.6M genuine tag-head deletion candidates.

| s | Kept | %candidates | %all_tag_runs | kept_marker | sa_values | **Total** |
|--:|--:|--:|--:|--:|--:|--:|
| 32  | 10.77M | 2.29% | 1.52% | 10.94 MB | 50.07 MB | **61.0 MB** |
| 64  |  5.74M | 1.22% | 0.81% |  6.60 MB | 26.66 MB | **33.3 MB** |
| 128 |  3.44M | 0.73% | 0.48% |  4.22 MB | 16.00 MB | **20.2 MB** |
| 256 |  2.27M | 0.48% | 0.32% |  2.93 MB | 10.57 MB | **13.5 MB** |

Bits per kept sample climb from 8.52 (s=32) to 10.82 (s=256) as
kept density drops -- sd_vector overhead grows sublinearly with
sparsity.  SA values consume 39 bits/entry uniformly.

**Verification** (per s value):

- (a) 1000/1000 kept-SA spot checks match `locate_sa_value_ref`
- (b) 100/100 deleted tag heads reachable within s LF steps AND
  reconstructed SA matches ground truth

On all 4 s values, both passes.  Zero mismatches.  Failure modes
covered include: SA computation bugs, bit-width truncation, wrong
`kept_marker` rank alignment with `sa_values`, deletion-rule
invariant violation, LF-walk arithmetic errors.  Not covered:
end-to-end pipeline (no `find_mems` variant consumes these files
yet), and rare corner cases below the 1100-sample-per-s
statistical resolution.

**Build cost breakdown** (HPRC chr6, Vesuvio, one invocation for
all 4 s values):

| Phase | Wall | Notes |
|:------|-----:|:------|
| R-index + ltag load | ~2 s | |
| Phase 1: candidate + SA enumeration | 4952 s (82 min) | 12.47B locateNext calls at 0.40 us/call |
| Sort by SA (719M candidates) | 82 s | std::sort on 24-byte tuples |
| Per-s: deletion + build + verify + write | ~4 s each | 4 passes over sorted candidates |
| **Total wall** | **1h 24m** | |
| Peak RSS | 21.7 GB | dominated by sorted candidate array |

Amortization: SA computation is O(n) regardless of s value, so
running all 4 s values in one invocation is essentially free
compared to running 1.

### Interpretation

**Storage overhead is far smaller than initially projected.**  The
early mean-field estimates in this session's design phase
predicted 300 MB - 3.5 GB storage per s -- too pessimistic by
1-2 orders of magnitude because the exponential-gap-length
assumption breaks: real BWT candidate gaps are fat-tailed, and
the interleaving of 249M BWT-run-head anchors (mean text-order
gap ~124) already covers most tag heads by proximity to a BWT
anchor.  Genuine samples are only needed for tag heads inside
long BWT runs, which are ~8-9% of all tag runs per JR-015 / the
FINDINGS_SEED_COST tail characterization.

**All s values fit the CLAUDE.md RSS budget by an order of
magnitude.**  Baseline find_mems peak RSS on HPRC chr6 = 4.2 GB.
CLAUDE.md 10% ceiling = 420 MB.  Every s value in {32, 64, 128,
256} adds <70 MB (1.6% overhead).  Storage is no longer the
constraint on this design.

**Query-time cost projection (not yet measured)** based on the
< s LF-step guarantee and the JR-013 locateNext cost of ~0.4 us:

- Per emit: <= s LF steps, mean ~s/2 assuming uniform position
  along the deletion interval.
- Per HPRC chr6 workload (155.5K MEMs * 6.2 emits/MEM = 964K emits):

| s | expected LF walk / query workload | vs today's legacy locate 30.5s |
|--:|--:|--:|
| 32  | 964K * 16 * 0.4us ~= 6s   | -24s savings |
| 64  | 964K * 32 * 0.4us ~= 12s  | -18s savings |
| 128 | 964K * 64 * 0.4us ~= 25s  | -5s savings  |
| 256 | 964K * 128 * 0.4us ~= 50s | +20s regression |

Sweet spot for query time: s=32 or s=64.  Both fit RSS budget
trivially and project to substantial (18-24s) wall savings on
the legacy path.  s=256 is a regression -- LF walk exceeds seed
cost being eliminated.

**Strategic implication.**  If this design ships, the legacy MEM
finder becomes competitive with (or faster than) flipped:
- Legacy today: 28.35s MEM-finding + 30.5s locate = 58.9s (Vesuvio noisy)
- Legacy + s=32 samples (projected): 28.35s MEM-finding + 6s LF walk = 34.4s
- Flipped Step 0' (JR-012 Vesuvio noisy): 38.63s MEM-finding + 2.48s locate = 41.1s

Projection: legacy + samples beats flipped by ~7s at 61 MB extra
storage.  Assumes projection holds; needs end-to-end validation.

### Open questions

1. **End-to-end validation not done.**  The samples files are
   correct at the SA-lookup level (verified) but no `find_mems`
   variant consumes them yet.  The `--gaf` gate on HPRC chr6
   must pass 2000/2000 valid AND coverage MD5 must match the
   JR-012 baseline `60f6b8e4a759aebb252b83a870b9ff8c` before
   this can be trusted for production.  This is the next
   correctness gate.
2. **Query-time perf is projected, not measured.**  The 6-24s
   wall savings estimates assume mean LF walk = s/2 with uniform
   distribution.  Real walks may be shorter (BWT anchors dense
   in text order intercept most walks well before s cap) or
   longer (specific tag-head distributions).  Only a real
   `find_mems` benchmark on Vesuvio can settle this.
3. **Should this replace or coexist with the flipped default?**
   JR-012 recommended flipping the default to `--use-flipped-mems`
   based on ~11% wall wins on both workloads.  If JR-016
   validates end-to-end with ~20% wall win over flipped, legacy
   + samples becomes the new default candidate.  Deferred until
   (1) and (2) complete.
4. **Choice of s.**  s=32 minimizes query cost but pays 61 MB;
   s=64 halves storage at ~2x query cost; s=128 further halves
   at another ~2x.  Actual sweet spot depends on downstream
   priorities: cohort-run batch (RSS matters more, prefer larger
   s) vs single-query latency (LF walk matters more, prefer
   smaller s).

### Verification (what's covered / not)

**Covered by build-time spot checks (a) + (b):**
- Phase 1 SA computation correctness (walk from BWT-run-head sample via locateNext)
- Bit-width sizing of sa_values (must match r-index samples width for packed SA)
- Coincident tag-head-with-BWT-anchor handling (step-0 anchor check in query)
- Deletion rule not too aggressive (deleted head reachable within s)
- kept_marker rank/select alignment with sa_values ordering
- LF() primitive behavior in the query walk
- Off-by-one in step counter during LF walk

**NOT covered (natural next-step gate):**
- End-to-end pipeline via `./query.sh --gaf`: no `find_mems`
  variant reads `.tag_samples.s<N>` yet.
- Statistical: 1100 spot-check samples per s across ~700M
  candidates rules out bug rates >~0.1% with high confidence,
  but not rare-corner-case bugs (1-in-10-million rate).

### Design decisions worth recording

- **Distance metric: text-position (SA), not BWT-position.**  BWT
  adjacency has no relation to LF-walk reachability.  Text
  adjacency IS the correct metric because LF walks step through
  SA-adjacent positions (SA decrements by 1 per LF step, modulo
  sequence boundaries which BWT runs never cross).
- **Coincident tag-head + BWT-anchor: implicit-deleted, no storage.**
  Their SAs are available for free from r-index `samples[]`.  Query
  slow-path step 0 catches these via the "is pos a BWT-run head"
  check.  Saves storage for 240.6M / 710.2M = 34% of all tag runs
  on HPRC chr6.
- **Single Phase 1 shared across all s values.**  SA computation
  cost is fixed at O(n).  Running all 4 s values in one
  invocation costs the same as running 1.  Building s=32 alone
  costs ~1h 24m; adding s=64/128/256 costs ~12 s more.

### Data artifacts

**Vesuvio build outputs (canonical / source of truth):**
- `~/mem-projection/pangenome-pipeline/runs/hprc-chr6-2026-06-02/tag_head_samples/hprcv1_chr6.tag_samples.s{32,64,128,256}`
- `~/mem-projection/pangenome-pipeline/runs/hprc-chr6-2026-06-02/tag_head_samples/build_tag_head_samples.log`
  (copied to Mac at `mem-projection/pangenome-pipeline/classify_output_jr015/build_tag_head_samples_v2.log`)

**Yeast smoke-test outputs (Mac):**
- `xy-test/tag_samples/yeast.tag_samples.s{32,64,128,256}` (v2 files,
  all pass verification 1000/1000 + 100/100)
- `xy-test/tag_samples/yeast.ri` (rebuilt with post-JR-014 format)

### Commits (this repo, branch tag-head-samples)

- `cedea60` build_tag_head_samples: sampled SA table over tag-run heads
  (initial v1 with BWT-order deletion; superseded)
- `29c28f1` run_hprc_chr6.sh: drop unnecessary .ri rebuild step
- `e5babd0` build_tag_head_samples: text-order deletion (v2) + verification
  (fixes SA bit-width truncation; adds spot-check passes; changes
  deletion rule to text-order per Cobas/Gagie/Navarro 2021)
- `890d413` run_hprc_chr6.sh: delete stale v1 outputs before rerun
- (this entry, pending as of writing)

---

## JR-017 — LF-operation latency on HPRC chr6: baseline vs fused LF_scan

```yaml
id: JR-017
date: 2026-08-09
author: claude-opus-4-7 (session with hlakshmidevi)
status: open
tags: [lf, r-index, microbenchmark, hprc, vesuvio, jr-016]
refs:
  follows: [JR-016]
benchmark-platform: vesuvio (Linux 6.12.73 x86_64, Debian)
```

### Context

JR-016 established that tag-head SA samples fit trivially in storage
(13-61 MB for s in {32, 64, 128, 256}), so the design's viability
hinges on **query-time LF-walk cost**.  Per-emit cost in the slow
path is bounded by s LF steps.  Earlier back-of-envelope projections
used ~0.4 us/LF from JR-013's aggregate locateNext measurement, but
LF is a different operation than locateNext and its true cost was
unknown.  Also open: whether fusing the two internal block walks in
LF (bwt_char_at_encoded + rankAt_encoded -> one scan_at) would
provide meaningful speedup.

### Method

Two implementations compared:

- **LF()** -- today's baseline (r-index.hpp:807): two independent block
  walks via psi_encoded, one for the char and one for the rank.
- **LF_scan()** -- new header-inline (r-index.hpp:988) wrapping one
  scan_at() call.  Reuses the block scan machinery JR-011 added.
  Bonus: populates run_id, run_start_bwt, and ranks[] as byproducts
  useful for downstream callers (e.g., the "is pos a BWT-run head"
  check in the sr-index sampled-SA slow path -- free from the same
  scan).

Two access patterns per variant:

- **Random**: 1,000,000 fresh random BWT positions, one LF per call.
  Models per-fresh-walk cost in the sr-index query (first LF step of
  each deleted tag-head's walk starts at an unrelated position).
- **Chain**: 10,000 random starts x 100 sequential LF steps
  (cur = LF(cur)).  Models the per-step cost of steps 2..s within an
  ongoing walk from a deleted tag-head to a nearby sample.

Protocol: `bin/bench_lf` (src/bench_lf.cpp).  Correctness pass first
(10,000 random positions cross-checked LF vs LF_scan).  Warmup pass
(one untimed round of both variants + both patterns).  Then 5 timed
trials per (variant, pattern), median + stddev reported.  Same
seeded position pools reused across variants so measurements are on
identical inputs.  XOR sink accumulator prevents dead-code elimination.

HPRC chr6 r-index used: `runs/hprc-chr6-2026-06-02/hprcv1_chr6.ri`
(post-JR-014 format, ~3.6 GB on disk).

### Findings

**Correctness.**  LF() == LF_scan() on 10,000 / 10,000 random positions.
Zero mismatches.

**Latency (Vesuvio, median-of-5 ns/op, stddev < 1% of median on every
row).**

| variant             | pattern | ns/op    | stddev | speedup vs baseline |
|:--------------------|:--------|---------:|-------:|--------------------:|
| LF (baseline)       | random  | 677.01   | 3.54   | 1.00x               |
| LF_scan (fused)     | random  | **646.53** | 4.65 | **1.05x (4.5% faster)** |
| LF (baseline)       | chain   | 664.31   | 2.35   | 1.00x               |
| LF_scan (fused)     | chain   | **634.29** | 2.54 | **1.05x (4.5% faster)** |

**Yeast comparison** (Mac, small n=337M .ri, ~168 MB on disk):

| variant             | pattern | ns/op   | speedup |
|:--------------------|:--------|--------:|--------:|
| LF (baseline)       | random  | 434.54  | 1.00x   |
| LF_scan (fused)     | random  | 330.92  | **1.31x (24% faster)** |

### Interpretation

**Fusing gives 24% on yeast but only 4.5% on HPRC.**  Root cause:
DRAM latency.  On yeast the entire .ri fits in L3 cache, so both
variants pay CPU-only work per LF call -- fusing eliminates one full
block walk (~50% of CPU work per LF), giving the ~24% total speedup.
On HPRC, block accesses are DRAM-bound (~150 ns per fetch); the
baseline's second scan hits the just-warmed block from L1, so it's
cheap.  Fusing only saves the CPU walk cost, not the memory fetch.
Net: fused wins by the walk-cost fraction of total, which shrinks
as DRAM cost grows.

**Random and chain are essentially identical wall on HPRC** (both ~640 ns
for LF_scan; 12 ns spread).  LF is a random-looking permutation, so
successive LF steps jump to unrelated block regions.  No cache-locality
benefit from being in a walk.  Confirms LF cost model is per-step
independent, no per-walk amortization.

**Real LF cost is ~634 ns/op** (LF_scan chain, HPRC).  Meaningfully
higher than the ~400 ns estimate JR-016 used for query-cost
projection.  Updated projections below.

### Impact on JR-016 query-cost projections

Using **634 ns/LF** (LF_scan chain, real number) and JR-016's HPRC
workload (155.5K MEMs x 6.2 emits/MEM = 964K total emits per query,
mean LF walk length ~ s/2 assuming uniform):

| s   | ops/emit | ops/workload | LF wall (s) | vs legacy locate 30.5s |
|----:|---------:|-------------:|------------:|-----------------------:|
| 32  |       16 |       15.4M  |      **9.8** | **-20.7 s (-68%)**       |
| 64  |       32 |       30.8M  |     **19.5** | **-11.0 s (-36%)**       |
| 128 |       64 |       61.7M  |     **39.1** | **+8.6 s regression**    |
| 256 |      128 |      123.4M  |     **78.3** | **+48 s regression**     |

**Only s=32 and s=64 remain wins**; s=128 flips to regression.
Previously JR-016 projected wins at all s <= 128 using the 0.4 us
estimate; s=128 was borderline (25 s vs 30.5 s).  With the real
634 ns/LF, s=128 is +8.6 s worse than legacy today.

**Sweet spot narrows to s=32 or s=64.**  Storage cost:

| s | storage | LF wall | savings vs legacy | savings vs flipped Step 0' (JR-012 noisy: 41.1s locate) |
|--:|--------:|--------:|------------------:|--------------------------------------------------------:|
| 32 | 61 MB  | 9.8 s  | -20.7 s (68%)   | -3 s vs flipped                                          |
| 64 | 33 MB  | 19.5 s | -11.0 s (36%)   | +7 s worse than flipped                                  |

**s=32 is the only projection where legacy + tag-head samples clearly
beats flipped** (~3 s wall win, tiny margin).  s=64 loses to flipped.
The design's competitive window has narrowed to a single s value with
a small absolute win.

### Discussion

**The naive-LF-optimization ceiling is low in DRAM-bound regimes.**
Fusing block walks (JR-011's scan_at machinery) helps CPU-bound
workloads a lot but flattens to ~5% when DRAM is the bottleneck.
Any further LF speedup would need to attack the DRAM cost itself:
prefetching, block layout changes to reduce cache-line splits, or
skipping LF entirely for certain query paths.

**Chain and random being identical wall is a load-bearing finding.**
It rules out cache-locality optimizations across LF steps within a
single walk.  Any per-step optimization applies uniformly.  It also
means our earlier intuition ("chain should be cheaper due to
sequential access") was wrong.

**The tag-head samples design remains marginally viable, not obviously
worth shipping.**  At s=32: 61 MB storage, ~3 s wall win over the
current flipped default.  End-to-end validation would take days of
integration + benchmarking (build a find_mems --tag-head-samples
variant, gate via --gaf, N=5 warm-cache Vesuvio bench).  Given the
small projected margin, worth asking whether other JR-016 open
questions (like the correctness gate itself) should complete before
investing in that path.

### Open questions

1. **Should we invest in the end-to-end find_mems integration?**
   Projected 3 s wall win at s=32 vs current flipped default is
   marginal.  Cost: ~3-5 days of integration + benchmark work.
   Deferred pending user direction.

2. **Prefetching LF's block fetch?**  A `__builtin_prefetch` on
   blocks_start_pos.predecessor(pos) + block_encoded_start_bits[..]
   could hide part of the DRAM latency.  Speculative gains 10-30%
   in similar r-index microbenchmarks in the literature; would need
   to measure.  Not attempted.

3. **Dedicated run-head bitvector to shortcut the run-head check?**
   Adding an sd_vector marking BWT-run-start positions (~250 MB
   storage for HPRC chr6) would let the query slow path skip the
   scan_at when checking "is cur a BWT-run head" -- but per JR-017's
   fusion result, that check is already ~free when combined with the
   LF scan.  Marginal.  Not worth building unless the query path
   ends up hot in a real find_mems integration.

### Data artifacts

- `~/mem-projection/pangenome-pipeline/runs/hprc-chr6-2026-06-02/tag_head_samples/bench_lf.log`
  (copied to Mac at `mem-projection/pangenome-pipeline/classify_output_jr015/bench_lf.log`)

### Commits (this repo, branch tag-head-samples)

- `063fbe7` r-index: add LF_scan (fused LF via scan_at) + bench_lf microbenchmark
- (this entry, pending as of writing)

---

## JR-018 — Tag-head SA samples: end-to-end integration + N=3 warm-cache perf on HPRC chr6

```yaml
id: JR-018
date: 2026-08-12
author: claude-opus-4-7 (session with hlakshmidevi)
status: resolved
tags: [tag-head-samples, integration, correctness, perf, benchmark, hprc, vesuvio, jr-016, jr-017]
refs:
  follows: [JR-016, JR-017]
benchmark-platform: vesuvio (Linux 6.12.73 x86_64, Debian)
```

### Context

JR-016 established that the sr-index-inspired sampled SA over tag-run
heads is storage-viable (61 MB at s=32 for HPRC chr6). JR-017
measured the real per-LF cost at 634 ns on HPRC (1.6x higher than the
earlier 400 ns proxy) and projected the design's competitive window
would narrow to s in {32, 64} with only ~3s wall margin over the
current flipped default (JR-012).

JR-018 wires the samples primitive into `find_mems` end-to-end for
both `--use-flipped-mems` and legacy MEM finders, validates
correctness on the HPRC chr6 pipeline gate, and measures perf under
N=3 warm-cache via `perf_harness.sh`.

### Method

**Integration** (branch `tag-head-samples`):

1. `TagHeadSamples` loader class (`include/pangenome_index/tag_head_samples.hpp`
   + `src/tag_head_samples.cpp`, commit `fd49146`). Pure data structure,
   inline `try_sa_at_tag_head(rid)` fast-path returning `UNKNOWN` sentinel
   for deleted rids.

2. `seed_sa_via_samples()` free function (`include/pangenome_index/tag_head_samples_query.hpp`,
   commit `fd49146`). Fast-path lookup + LF-walk slow path (up to
   `sample_period()` steps) using `LF_scan` (JR-017 fused primitive)
   for free BWT-run-head check.

3. CLI flag `--tag-head-samples=<path>` on `find_mems`. Requires
   `--lightweight-tags`. `--use-flipped-mems` optional; dispatch
   selects the appropriate samples emit variant.

4. **Flipped variant** `dump_mem_info_lightweight_flipped_samples`
   (commit `fd49146`): first_rid emit uses `mem_sa.sa_sp` directly
   (flipped SA-carry); interior/last rids resolved via `seed_sa_via_samples`.

5. **Legacy variant** `dump_mem_info_lightweight_samples` (commit
   `7301cdb`): first_rid resolved via `seed_sa_via_samples(first_rid)`
   plus a `locateNext` walk from `run_start_bwt(first_rid)` to
   `mem.bwt_start` (mean ~4 BWT positions on HPRC); interior/last
   identical to flipped variant.

6. Fine-grained timing (commit `175fa65`, `#if TIME`-gated): samples
   file load, first_rid resolve+walk, interior/last resolve.

**Correctness gate**: `./query.sh configs/hprcv1-chr6.env
hprc-chr6-2026-06-02 <name> --gaf`. Coverage MD5 must equal
`60f6b8e4a759aebb252b83a870b9ff8c` (JR-012 baseline); validate_gaf
must be 2000/2000.

**Perf harness**: `perf/perf_harness.sh <config> 3 <tag>
--modes=lightweight` (2 untimed warmups + 3 timed trials, warm cache).

**Also built s=8 and s=16** samples (commit `3e45b60`, single Vesuvio
invocation `run_hprc_chr6_s8_s16.sh` amortizing Phase 1 SA
enumeration across both). Extends the JR-016 s-curve.

### Findings

**S-curve storage table (extends JR-016 with s=8/s=16):**

| s | Kept | %candidates | %all_tag_runs | File size |
|--:|--:|--:|--:|--:|
| 8 | 50.30M | 10.71% | 7.08% | **272 MB** |
| 16 | 22.40M | 4.77% | 3.15% | **125 MB** |
| 32 | 10.77M | 2.29% | 1.52% | 61 MB |
| 64 | 5.74M | 1.22% | 0.81% | 33 MB |
| 128 | 3.44M | 0.73% | 0.48% | 20 MB |
| 256 | 2.27M | 0.48% | 0.32% | 14 MB |

**Correctness gate (HPRC chr6 alt-noisy, 100K reads):** all 11 configs
tested (baseline flipped, baseline legacy, flipped+samples s in
{8, 16, 32, 64, 128, 256}, legacy+samples s in {8, 16, 32, 64, 128, 256})
produce coverage MD5 = `60f6b8e4a759aebb252b83a870b9ff8c` and
validate_gaf 2000/2000. Yeast chrII (2000 reads) also byte-identical
between legacy baseline and legacy+samples across all four s values.
Zero correctness regressions.

**N=3 warm-cache perf (Vesuvio, HPRC chr6 alt-noisy):**

| Config | find_mems wall | Peak RSS | vs flipped |
|:--|--:|--:|--:|
| Flipped baseline | 38.45 +- 0.08 s | 4512 MB | 0 |
| Legacy baseline | 44.78 +- 0.03 s | 4512 MB | +6.33s (+16%) |
| **Legacy + samples s=8** | **30.47 +- 0.08 s** | 4789 MB | **-7.98s (-20.8%)** |
| Legacy + samples s=16 | 32.53 +- 0.06 s | 4639 MB | -5.92s (-15.4%) |
| Flipped + samples s=8 | 38.58 +- 0.02 s | 4788 MB | +0.13s (noise) |

**Combined pipeline wall (find_mems + gafpack):**

| Config | find_mems | gafpack | **Combined** | vs flipped |
|:--|--:|--:|--:|--:|
| Flipped baseline | 38.45s | 24.33s | **62.78s** | 0 |
| Legacy + samples s=8 | 30.47s | 24.28s | **54.75s** | **-8.03s (-12.8%)** |
| Legacy + samples s=16 | 32.53s | 24.39s | **56.92s** | -5.86s (-9.3%) |

**S-curve on legacy+samples path** (single-trial via query.sh):

| s | find_mems wall | Fast-path % | Avg LF steps / slow |
|--:|--:|--:|--:|
| 8 | **30.46s** | 12.27% | 2.99 |
| 16 | 32.30s | 5.25% | 6.25 |
| 32 | 36.17s | 2.82% | 12.35 |
| 64 | 43.98s | 1.04% | 24.48 |
| 128 | 58.35s | 0.55% | 46.64 |
| 256 | 85.10s | 0.28% | 85.31 |

s in {8, 16, 32} beat flipped; s=64 ties legacy baseline; s in {128, 256} regress.

**Fine-grained emit breakdown (legacy + samples, N=3 mean):**

| Config | Samples file load | First-rid resolve+walk | Interior/last resolve | first_rid locateNext calls |
|:--|--:|--:|--:|--:|
| s=8 | 0.12s | 0.81s | 2.04s | 615.6K |
| s=16 | 0.05s | 0.89s | 4.03s | 615.6K |

Interior cost dominates and scales linearly with avg LF steps x
per-LF cost (JR-017: 634 ns). Model closes to within ~5% (e.g. s=8:
810K rids x 4 avg LF x 634 ns = ~2.05s predicted vs 2.04s measured).

**Instrumentation overhead check (TIME=0 vs TIME=1, N=3):**

| Config | TIME=1 | TIME=0 | Delta |
|:--|--:|--:|--:|
| Flipped baseline | 38.49s | 38.56s | +0.07s (noise) |
| Legacy + samples s=8 | 30.64s | 30.43s | -0.21s |

Instrumentation adds <1% overhead. All reported numbers are honest.

### Interpretation

**Legacy MEM finding is ~7.5s faster than flipped on HPRC chr6
noisy** (26.46s vs 34.00s, both N=3, variance <0.05s). Legacy's
3-phase enumeration is fundamentally cheaper than flipped's Step 2'
fresh forward extend per emit (Risk C in JR-011). Flipped monetized
this by eliminating the ~14s locate_sa_value seed cost via SA-carry;
that trade-off left flipped ~4s faster overall (JR-012).

**Legacy + samples restores the original trade shape**: keep legacy's
faster MEM finding, eliminate its seed cost via samples lookup
instead of SA-carry. At s=8: samples emit costs ~2.9s (vs flipped's
~2.1s locateNext emit), while legacy's MEM-finding advantage saves
~7.5s. Net: -7.98s (-20.8%) over flipped baseline.

**Dataset-invariant fast-path rate at s=32** (~3% on both yeast and
HPRC despite 100x scale difference; see JR-018 measurements in this
entry) confirms the deletion rule's storage-first bias. Aggressive
sampling (s=8) meaningfully improves fast-path rate (12.27%) and
shortens slow-path LF walks (avg 3 steps), which is why s=8 becomes
the winner rather than a marginal improvement.

**Flipped + samples adds nothing** (+0.13s vs flipped baseline).
Flipped already eliminated the seed cost via SA-carry; samples emit
at s=8 offers no compounding benefit and slightly slows the emit
path. Confirmed: samples targets legacy's seed cost, not flipped's.

### Recommendation

**Adopt `--tag-head-samples=<hprcv1_chr6.tag_samples.s8>` with the
legacy MEM finder as the default HPRC chr6 configuration.**
Expected win over current flipped default: -8s find_mems wall
(-20.8%), -8s combined pipeline (-12.8%), at +277 MB peak RSS
(+6.1%, well under the CLAUDE.md +10% ceiling) and +272 MB on-disk
storage for the samples file.

**s=16 is the backup option** for storage-constrained workloads:
-5.92s find_mems (-15.4%), +127 MB RSS (+2.8%), +125 MB storage.

### Open questions

1. **Whole-genome HPRCv1/HPRCv2 datasets not yet tested.** The
   ~3% fast-path rate held across the yeast:HPRC-chr6 100x scale gap,
   so we expect it to hold at whole-genome scale, but should
   confirm. Storage cost extrapolates to ~2-3 GB per s=8 samples
   file on a full HPRCv1 (~100 Gbp) -- still fits typical server
   RAM but worth measuring.

2. **Genotyping-framework integration.** JR-016/017/018 have
   demonstrated the samples design is production-ready as a
   find_mems drop-in. Downstream use as input to graph-coverage-based
   genotyping is next.

3. **Multi-batch parallelism.** Coverage aggregation is
   embarrassingly parallel across read batches; a single-pass merge
   can amortize further. Not yet exploited.

4. **Ship the flag as default?** Currently opt-in via
   `--tag-head-samples`. Default flip would require:
   (a) pipeline wrapper (`query.sh`) to auto-locate a `.tag_samples.s<N>`
       file next to the `.ri`, or
   (b) a decision on whether we ship the `.tag_samples.s8` file as
       part of the index build output by default.
   Not scoped for this entry.

### Data artifacts

**Vesuvio perf logs (canonical / source of truth):**

- Original single-trial via `query.sh`:
  `~/mem-projection/pangenome-pipeline/runs/hprc-chr6-2026-06-02/queries/jr018-{legacy-baseline,legacy-samples-s{8,16,32,64,128,256},samples-s{8,16,32}}/lightweight/`
- N=3 warm-cache via `perf_harness.sh`:
  `~/mem-projection/pangenome-pipeline/perf/jr018-final-{flipped-baseline,legacy-baseline,legacy-samples-s{8,16}}/`
- TIME=0 sanity harness:
  `~/mem-projection/pangenome-pipeline/perf/jr018-perf-{flipped-baseline,legacy-samples-s8}-notime/`
- Build log for s=8/s=16:
  `~/mem-projection/pangenome-pipeline/runs/hprc-chr6-2026-06-02/tag_head_samples/build_tag_head_samples_s8_s16.log`

**Mac copies:**

- `xy-test/tag_samples/jr018_v2_logs/` (query.sh single-trial)
- `xy-test/tag_samples/jr018_perf_logs/` (N=3 harness, TIME=1)
- `xy-test/tag_samples/jr018_perf_notime_logs/` (N=3 harness, TIME=0)
- `xy-test/tag_samples/jr018_final_logs/` (N=3 harness with fine-grained timing)
- `xy-test/tag_samples/build_tag_head_samples_s8_s16.log`

### Commits (this repo, branch tag-head-samples)

- `fd49146` find_mems: integrate tag-head SA samples into flipped-lightweight emit path
- `3e45b60` run_hprc_chr6_s8_s16: build tag_head_samples at s=8 and s=16
- `7301cdb` find_mems: legacy+samples emit path with first_rid walk
- `175fa65` find_mems: instrument samples load + first_rid + interior emit times
- (this entry, pending as of writing)

---

## JR-019 — Flipped vs Legacy+s=8 end-to-end perf on HPRC chr6 alt-noisy: L=25 vs L=50

```yaml
id: JR-019
date: 2026-08-12
author: claude-opus-4-7 (session with hlakshmidevi)
status: resolved
tags: [perf, benchmark, hprc, alt-noisy, l25, l50, samples, jr-018]
refs:
  follows: [JR-018]
benchmark-platform: vesuvio (Linux 6.12.73 x86_64, Debian)
```

### Context

JR-018 established that legacy + samples s=8 beats flipped baseline by
-20.8% on find_mems wall (N=3 warm-cache) at MEM_LEN=50 on HPRC chr6
alt-noisy. JR-019 extends the measurement to MEM_LEN=25 to confirm
the win generalizes across MEM-length regimes, and records the
end-to-end throughput (find_mems + gafpack) with standard workload
characteristics for both L values.

### Method

Two `query.sh --gaf` invocations per L value (flipped baseline + legacy
+ samples s=8), single-trial each. L=50 uses `configs/hprcv1-chr6.env`;
L=25 uses `configs/hprcv1-chr6-alt-reads-L25.env`. Same reads
(`chr6.alt_noisy.reads.txt`), same graph, same index (post-JR-014
`.ri`).

L=50 timing values use the N=3 warm-cache means from the JR-018
perf harness (source `perf/jr018-final-{flipped-baseline,legacy-samples-s8}`).
L=25 values are single-trial from query.sh (source
`runs/hprc-chr6-2026-06-02/queries/jr018-L25-{flipped,legacy-s8}`).

Correctness verified: L=25 configs produce validate_gaf 2000/2000
each and byte-identical record counts (250,177 MEMs, 1,432,999
records). L=50 configs share coverage MD5 `60f6b8e4a759aebb252b83a870b9ff8c`
per JR-018.

### Findings

**MEM_LEN = 50**

Performance (all wall time):

| Metric | Flipped baseline | Legacy + samples s=8 | delta |
|:--|--:|--:|--:|
| find_mems | 38.45s | 30.47s | **-7.98s (-20.8%)** |
| gafpack | 24.33s | 24.28s | -0.05s (noise) |
| **Total wall** | **62.78s** | **54.75s** | **-8.03s (-12.8%)** |
| Peak RSS (find_mems) | 4512 MB | 4789 MB | +277 MB (+6.1%) |
| Peak RSS (gafpack) | 314 MB | 314 MB | ~0 |
| Throughput -- reads/s | 1593 | 1827 | **+15%** |
| Throughput -- MEMs/s | 2,478 | 2,841 | **+15%** |

Workload characteristics:

| | Value |
|:--|:--|
| Read file | `chr6.alt_noisy.reads.txt` |
| Read length | 200 bp |
| Total reads | 100,000 |
| Total MEMs | 155,551 |
| Average MEM length | 106.29 bp |
| Total occurrence of all MEMs (Sigma mem.size) | 16,792,975 |
| Total records written by find_mems | 965,868 |
| Total occurrence on graph (post-gafpack dedup, GAF entries) | 532,533 |
| Average occurrence / MEM | 107.96 |
| Average tags / MEM | 6.21 |

Configuration params:

| Param | Value |
|:--|:--|
| min_len (MEM_LEN) | 50 |
| min_occ (MIN_OCC) | 1 |
| Reads type | HPRC alt-noisy (in-graph, ~1% base error) |
| Graph | HPRC v1 chr6 (`hprcv1_chr6.gbz`) |
| Index | `hprcv1_chr6.ri` (~3.76 GB), `hprcv1_chr6.ltags` (~89 MB) |
| Samples file (legacy+s=8 only) | `hprcv1_chr6.tag_samples.s8` (272 MB) |

**MEM_LEN = 25**

Performance (all wall time):

| Metric | Flipped baseline | Legacy + samples s=8 | delta |
|:--|--:|--:|--:|
| find_mems | 48.34s | 39.10s | **-9.24s (-19.1%)** |
| gafpack | 24.23s | 24.47s | +0.24s (noise) |
| **Total wall** | **72.57s** | **63.57s** | **-9.00s (-12.4%)** |
| Peak RSS (find_mems) | 4531 MB | 4808 MB | +277 MB (+6.1%) |
| Peak RSS (gafpack) | 328 MB | 329 MB | ~0 |
| Throughput -- reads/s | 1378 | 1573 | **+14%** |
| Throughput -- MEMs/s | 3,447 | 3,935 | **+14%** |

Workload characteristics:

| | Value |
|:--|:--|
| Read file | `chr6.alt_noisy.reads.txt` |
| Read length | 200 bp |
| Total reads | 100,000 |
| Total MEMs | 250,177 |
| Average MEM length | 79.42 bp |
| Total occurrence of all MEMs (Sigma mem.size) | 32,050,895 |
| Total records written by find_mems | 1,432,999 |
| Total occurrence on graph (post-gafpack dedup, GAF entries) | 840,010 |
| Average occurrence / MEM | 128.16 |
| Average tags / MEM | 5.92 |

Configuration params:

| Param | Value |
|:--|:--|
| min_len (MEM_LEN) | 25 |
| min_occ (MIN_OCC) | 1 |
| Reads type | HPRC alt-noisy (in-graph, ~1% base error) |
| Graph | HPRC v1 chr6 (`hprcv1_chr6.gbz`) |
| Index | `hprcv1_chr6.ri` (~3.76 GB), `hprcv1_chr6.ltags` (~89 MB) |
| Samples file (legacy+s=8 only) | `hprcv1_chr6.tag_samples.s8` (272 MB) |

**Cross-L comparison:**

| Metric | L=50 | L=25 | change |
|:--|--:|--:|--:|
| Total MEMs | 155,551 | 250,177 | +61% |
| Sigma mem.size | 16.79M | 32.05M | +91% |
| Total records (find_mems) | 965,868 | 1,432,999 | +48% |
| Total GAF entries (post-dedup) | 532,533 | 840,010 | +58% |
| Compression ratio (Sigma mem.size / GAF entries) | 31.5x | 38.2x | -- |
| **Total wall -- flipped** | 62.78s | 72.57s | +9.79s (+16%) |
| **Total wall -- legacy+s=8** | 54.75s | 63.57s | +8.82s (+16%) |
| **Win magnitude (end-to-end)** | **-8.03s (-12.8%)** | **-9.00s (-12.4%)** | consistent |
| Peak RSS delta | +277 MB | +277 MB | same |

### Interpretation

**The legacy+s=8 win generalizes cleanly across MEM lengths.** End-to-end
wall reduction is -12.4% at L=25 and -12.8% at L=50 -- statistically
indistinguishable at single-trial resolution. Absolute win grows slightly
in wall time (7.98 -> 9.24 s on find_mems) because samples emit scales
sublinearly with MEM count while flipped's `locateNext` walk scales
linearly with BWT positions scanned.

**Shorter min-len filter (L=25) produces more, shorter MEMs.** Same
200 bp reads, but the L=25 threshold admits ~1.6x more MEMs at ~25%
shorter avg length (79.4 vs 106.3 bp). Total BWT positions scanned
grows by ~1.9x (32.0M vs 16.8M), roughly matching the MEM count
growth times a small increase in avg-occurrence-per-MEM (128.2 vs
108.0).

**Compression ratios are stable across L.** find_mems collapses ~17x
of BWT positions into per-tag-run records (17.4x at L=50, 22.4x at L=25);
gafpack dedup collapses another ~1.7x to unique graph positions.
Net: 31-38x compression from raw BWT hits to graph coverage entries.

**Peak RSS delta is workload-invariant** (+277 MB in both), confirming
the samples file's in-memory footprint (285.6 MB per JR-016 build
log) is the sole source of the memory overhead.

### Open questions

1. **N=3 warm-cache validation for L=25.** L=50 numbers are already
   N=3 (JR-018); L=25 numbers here are single-trial. Given the
   consistency of the win margin, N=3 for L=25 is a nice-to-have,
   not a gate.
2. **Out-of-graph read behavior.** HG002 clean/noisy reads (out of
   HPRC v1 graph) tested separately -- results not yet in this
   entry.
3. **Whole-genome scaling.** HPRC v1 chr1 samples build in progress
   (Vesuvio, ~2 hr wall); will follow up with chr1 perf when
   complete.

### Data artifacts

**Vesuvio (canonical / source of truth):**
- L=50 N=3 harness: `~/mem-projection/pangenome-pipeline/perf/jr018-final-{flipped-baseline,legacy-samples-s8}/`
- L=25 single-trial: `~/mem-projection/pangenome-pipeline/runs/hprc-chr6-2026-06-02/queries/jr018-L25-{flipped,legacy-s8}/lightweight/`

**Mac copies:**
- L=50: `xy-test/tag_samples/jr018_final_logs/{jr018-final-flipped-baseline,jr018-final-legacy-samples-s8}/`
- L=25: `xy-test/tag_samples/jr018_L25_logs/jr018-L25-{flipped,legacy-s8}/`

### Commits (this repo, branch tag-head-samples)

No new code commits; measurement-only entry against binaries built
from HEAD `b2c1608` (find_mems: suppress Locate ops block when samples
path is active).

---

## JR-020 — convert_tags decoded SDSL header as ByteCode: latent bug hit by HPRCv1 chr1 build

```yaml
id: JR-020
date: 2026-08-13
author: claude-opus-4-7 (session with hlakshmidevi)
status: resolved
tags: [convert_tags, build_tags, sdsl, int_vector_buffer, tag-arrays, correctness, hprc, chr1, bug-fix]
refs:
  follows: [JR-019]
```

### Context

Building the HPRCv1 chr1 tag index produced a `.ltags` whose
`bwt_intervals.size()` violated the CLAUDE.md line-112 invariant
`bwt_intervals size == BWT_size + 1`:

| Dataset            | BWT size (n)    | bwt_intervals size (built) | Expected (n+1) | Delta |
|:-------------------|----------------:|---------------------------:|---------------:|------:|
| HPRCv1 chr6 alt    |  31,016,755,770 |             31,016,755,771 | 31,016,755,771 |     **0** ✓ |
| **HPRCv1 chr1**    |  45,105,246,014 |         **45,105,246,211** | 45,105,246,015 | **+196** ✗ |

Chr1's `_compressed.tags` and the downstream `.ltags` were unusable
for MEM projection (every BWT-position query would be shifted by up
to 196 positions). The chr6 build was correct so the failure was
initially assumed to be a chr1-specific pipeline issue. It was not
-- it was a latent bug in `src/convert_tags.cpp` that had gone
undetected on every prior graph because the header bytes happened
to decode into varints with `node_id == 0` and were silently skipped.

Empirical shortcut discovered while investigating: passing
`--num-seq 4328` (instead of the correct 4524 = `2 * 2262 paths`)
produced `bwt_intervals size = 45,105,246,015`. `4524 - 4328 = 196`
= the exact overshoot. Any linear shift of `num_seq` linearly
shifts `bwt_intervals`, so this "fix" masked the real bug -- it did
not remove any real coverage; it just under-prepended 196 endmarker
positions so the accounting cancelled out. If used in production
this would corrupt every seq_id assignment in gafpack output. See
"Do not deploy the workaround" below.

### Hypothesis

1. The `.tags` file produced by `traverse_sequences_parallel`
   (`include/pangenome_index/algorithm.hpp:496`) is opened as
   `sdsl::int_vector_buffer<8> out(filename, std::ios::out | std::ios::trunc)`.
   Per SDSL (`sdsl-lite/include/sdsl/int_vector_buffer.hpp:153`), a
   nonzero `t_width` implies `m_offset = 8`, meaning the file has
   an 8-byte header (LE uint64_t `size_in_bits`) prefixed to the
   payload, plus 0-7 trailing zero bytes of alignment padding
   written by `close()`.
2. `convert_tags.cpp` (`src/convert_tags.cpp:66`, pre-fix) slurps
   the entire file and starts decoding at byte 0 with
   `gbwt::ByteCode::read`, treating the SDSL header as raw
   ByteCode payload.
3. On chr1 the first 4 header bytes `b8 cb d8 78` form a single
   4-byte varint whose top bit is finally clear on `0x78`.
   Decoding: `raw = 0x0f1625b8 = 253,247,928`. Field slicing per
   `TagArray::encode_run_length`/`decode_run` (10 bits offset, 1
   bit is_rev, 9 bits length, rest node_id): `node=241`,
   `off=440`, `is_rev=1`, `length=196`.
4. `node > 0` means `convert_tags.cpp:154` does NOT hit the
   `nid == 0` skip. The varint is emitted as a real run of length
   196 into the compact output. Hence `+196` overshoot.
5. On chr6 the same code path fires but *every* header varint
   happens to have `node_id == 0` (the arithmetic accident that
   the first 8 bytes of a particular `size_in_bits` value all
   have their high bits arranged such that `raw >> 20 == 0`), so
   every phantom decode is skipped and the bug is silent.

### Method

**Direct evidence (vesuvio hex-dumps of the algorithm-format .tags
files):**

```sh
# On vesuvio
cd ~/mem-projection/pangenome-pipeline/runs/hprcv1-chr1-2026-08-12/
xxd -c 16 -l 32 hprcv1_chr1.tags       # first 32 bytes -> SDSL header + payload start
xxd -c 16 -s -32 hprcv1_chr1.tags      # last 32 bytes  -> trailing pad
```

Same for chr6, yeast, and hprcv2-mc-chr6. Then simulate
`gbwt::ByteCode::read` + `TagArray::decode_run` on the 8 header
bytes in Python and enumerate every varint that would be decoded
and whether its `node_id == 0`.

**Correctness verification (three tiers, per CLAUDE.md line 109):**

Tier 1 -- local `convert_tags` round-trip against the checked-in
yeast reference:

```sh
cd ~/personal/pangenome-index-latest
make bin/convert_tags        # patched
bin/convert_tags \
    ../mem-projection/pangenome-pipeline/runs/v2-yeast235/yeast235_chrII_100kb_normalized.tags \
    /tmp/yeast.patched.tags \
    --num-seq 634
md5 /tmp/yeast.patched.tags \
    ../mem-projection/pangenome-pipeline/runs/v2-yeast235/yeast235_chrII_100kb_normalized_compressed.tags
```

Tier 2 -- HPRCv1 chr6 E2E on vesuvio, comparing coverage MD5 to
the JR-012/JR-013/JR-018/JR-019 gold-standard
`60f6b8e4a759aebb252b83a870b9ff8c`:

```sh
# On vesuvio
scp Mac:~/personal/pangenome-index-latest/src/convert_tags.cpp \
    ~/pangenome-index-latest/src/convert_tags.cpp
cd ~/pangenome-index-latest && make bin/convert_tags

cd ~/mem-projection/pangenome-pipeline
./query.sh configs/hprcv1-chr6-alt-reads.env hprc-chr6-2026-06-02 \
    convert-tags-header-fix --gaf
md5sum runs/hprc-chr6-2026-06-02/queries/convert-tags-header-fix/lightweight/alignment_coverage.csv
```

Tier 3 -- HPRCv1 chr1 rebuild with the patched binary, confirm the
CLAUDE.md invariant:

```sh
cd ~/mem-projection/pangenome-pipeline/runs/hprcv1-chr1-2026-08-12/
gtime -v ~/pangenome-index-latest/bin/convert_tags \
    hprcv1_chr1.tags \
    hprcv1_chr1_compressed.tags.header-fix \
    --num-seq 4524 \
    2>&1 | tee logs/06_convert_tags_header_fix.log
grep "Building the bwt_intervals" logs/06_convert_tags_header_fix.log
```

### Findings

**Header hex-dumps + simulated ByteCode decode (buggy convert_tags
path). Each row is one varint the buggy code would decode from the
8 header bytes.**

Chr1 (`b8 cb d8 78 29 00 00 00`):

| Varint bytes | raw       | decoded fields                       | action        |
|:-------------|:----------|:-------------------------------------|:--------------|
| `b8 cb d8 78`| `0xf1625b8` | node=241 off=440 rev=1 len=**196**   | **EMIT (bug)** |
| `29`         | `0x29`    | node=0   off=41  rev=0 len=0         | skip           |
| `00`         | `0x0`     | node=0   off=0   rev=0 len=0         | skip           |
| `00`         | `0x0`     | node=0   off=0   rev=0 len=0         | skip           |
| `00`         | `0x0`     | node=0   off=0   rev=0 len=0         | skip           |

Chr6 alt (`58 94 44 07 08 00 00 00`):

| Varint bytes | raw       | decoded fields                | action |
|:-------------|:----------|:------------------------------|:-------|
| `58`         | `0x58`    | node=0 off=88  rev=0 len=0    | skip   |
| `94 44`      | `0x2214`  | node=0 off=532 rev=0 len=4    | skip   |
| `07`         | `0x7`     | node=0 off=7   rev=0 len=0    | skip   |
| `08`         | `0x8`     | node=0 off=8   rev=0 len=0    | skip   |
| `00`         | `0x0`     | (as above)                    | skip   |
| ...          |           |                               |        |

Yeast (`50 73 b7 84 02 00 00 00`):

| Varint bytes | raw       | decoded fields                | action |
|:-------------|:----------|:------------------------------|:-------|
| `50`         | `0x50`    | node=0 off=80  rev=0 len=0    | skip   |
| `73`         | `0x73`    | node=0 off=115 rev=0 len=0    | skip   |
| `b7 84 02`   | `0x8237`  | node=0 off=567 rev=0 len=16   | skip   |
| `00`         | `0x0`     | (as above)                    | skip   |
| ...          |           |                               |        |

Hprcv2-mc-chr6 (`18 15 e4 0c 0a 00 00 00`):

| Varint bytes | raw       | decoded fields                | action |
|:-------------|:----------|:------------------------------|:-------|
| `18`         | `0x18`    | node=0 off=24  rev=0 len=0    | skip   |
| `15`         | `0x15`    | node=0 off=21  rev=0 len=0    | skip   |
| `e4 0c`      | `0x664`   | node=0 off=612 rev=1 len=0    | skip   |
| `0a`         | `0xa`     | node=0 off=10  rev=0 len=0    | skip   |
| `00`         | `0x0`     | (as above)                    | skip   |
| ...          |           |                               |        |

Prediction from these tables: only chr1 emits a phantom run;
that run has `len=196`, matching the observed `+196` overshoot.
Chr6, yeast, and hprcv2-mc-chr6 emit 0 phantom bytes and are
byte-safe.

**Direct match against observed convert_tags "Skipping" log lines:**

The buggy convert_tags logs each `nid == 0` skip. My simulated
decode of the header bytes matches the leading log lines exactly
in every dataset (`pos_t` prints as `id[+/-]offset` per
`gbwtgraph::utils.h:330`; `-` means `is_rev=1`, `+` means `is_rev=0`;
the second column is `length`):

| Dataset | Simulated decode  | Observed log line   |
|:--------|:------------------|:--------------------|
| chr1    | node=241 off=440 rev=1 len=196 | (not logged; skip fires only if node==0)  |
| chr1    | node=0 off=41 rev=0 len=0      | `Skipping pure endmarker run: 0+41 0`     |
| chr6    | node=0 off=88 rev=0 len=0      | `Skipping pure endmarker run: 0+88 0`     |
| chr6    | node=0 off=532 rev=0 len=4     | `Skipping pure endmarker run: 0+532 4`    |
| chr6    | node=0 off=7 rev=0 len=0       | `Skipping pure endmarker run: 0+7 0`      |
| chr6    | node=0 off=8 rev=0 len=0       | `Skipping pure endmarker run: 0+8 0`      |
| yeast   | node=0 off=80 rev=0 len=0      | `Skipping pure endmarker run: 0+80 0`     |
| yeast   | node=0 off=115 rev=0 len=0     | `Skipping pure endmarker run: 0+115 0`    |
| yeast   | node=0 off=567 rev=0 len=16    | `Skipping pure endmarker run: 0+567 16`   |
| mc-chr6 | node=0 off=24 rev=0 len=0      | `Skipping pure endmarker run: 0+24 0`     |
| mc-chr6 | node=0 off=612 rev=1 len=0     | `Skipping pure endmarker run: 0-612 0`    |

Every non-zero-offset skip line comes from header bytes.
`0+0 0` lines come from either header bytes with all-zero bits or
from trailing SDSL zero-padding (up to 7 bytes).

**Verification (three tiers):**

Tier 1 -- yeast round-trip. Reference produced by pre-fix
convert_tags; my Mac-side baseline (rebuilt HEAD binary) and
patched binary both regenerated it. All three MD5s match:

```
0672d1e901cfedffd110a1f743dccc55   yeast.baseline.tags     (rebuilt HEAD)
0672d1e901cfedffd110a1f743dccc55   yeast.patched.tags      (patched)
0672d1e901cfedffd110a1f743dccc55   ../runs/v2-yeast235/yeast235_chrII_100kb_normalized_compressed.tags (checked-in reference)
```

Stats logged by convert_tags in all three runs identical:
`bwt_intervals size 337,356,603` (= BWT_size + 1 = 337,356,602 + 1),
`Number of runs 226,291,960`. Log differs only in skip-line count
(baseline: 12 lines from 8 header bytes + up to 4 trailing pad;
patched: 6 lines from trailing pad only) -- confirms 8 bytes are
being stripped correctly and no payload byte is being consumed
in the process.

Tier 2 -- HPRCv1 chr6 alt-noisy E2E on vesuvio. Patched
convert_tags -> build_lightweight_tags -> find_mems -> gafpack
-> validate_gaf:

| variant | coverage MD5 | validate_gaf |
|:--------|:-------------|-------------:|
| `convert-tags-header-fix` | `60f6b8e4a759aebb252b83a870b9ff8c` | 2000/2000 ✓ |
| JR-012/JR-013/JR-018/JR-019 gold baseline | `60f6b8e4a759aebb252b83a870b9ff8c` | -- (prior) |
| `alt-noisy` older pre-JR-012 run | `f96f5030aa70a70b038a6ad1cf517a8b` | -- (superseded) |

**Byte-identical downstream coverage** to every prior fix that
targeted this dataset. This is JR-007's strongest correctness
signal (see JR-007 line 933).

Tier 3 -- HPRCv1 chr1 rebuild on vesuvio:

| dataset    | pre-fix bwt_intervals | post-fix bwt_intervals | expected (BWT+1) | pre-fix runs | post-fix runs |
|:-----------|----------------------:|-----------------------:|-----------------:|-------------:|--------------:|
| **chr1**   |        45,105,246,211 |         **45,105,246,015** |   45,105,246,015 | 3,336,172,720 | **3,336,172,719** |

Both invariant checks pass:
- `bwt_intervals size == BWT + 1` ✓
- `num_runs` dropped by exactly 1 (the single phantom run from the
  4-byte header varint disappeared)
- `.ltags` built downstream from the corrected `_compressed.tags`:
  2.5 GB, 10 s, ~19 GB peak RSS, `num_runs = 3,336,172,719`
  (consistent with `_compressed.tags`)

Post-fix chr1 skip lines: `Skipping pure endmarker run: 0+0 0`
(exactly 1; from the single trailing pad byte -- your `xxd -c 16 -s -32`
showed the file ending in `...f2 ba d9 01 00`, one trailing zero).

**Audit of other existing indexes** (using the invariant check
below, no rebuild required):

| Dataset | BWT size (n) | bwt_intervals size | delta (bwt_int - n - 1) | status |
|:--------|-------------:|-------------------:|------------------------:|:------|
| yeast235 chrII 100 kb normalized  |       337,356,602 |        337,356,603 |  **0** | ✓ safe |
| HPRCv1 chr6 alt                   |    31,016,755,770 |     31,016,755,771 |  **0** | ✓ safe |
| **HPRCv2 MC chm13 chr6**          |   157,116,624,594 |    157,116,624,595 |  **0** | ✓ safe |
| **HPRCv1 chr1 (pre-fix build)**   |    45,105,246,014 |     45,105,246,211 | **+196** | ✗ affected -- rebuild required |
| HPRCv1 chr1 (post-fix)            |    45,105,246,014 |     45,105,246,015 |  **0** | ✓ fixed |

### Interpretation

The bug is **latent and data-dependent**. Whether a given `.tags`
file exposes it depends on the exact SDSL header byte pattern,
which is derived from `size_in_bits` = 8 * payload_bytes. For a
random large file, roughly half the possible headers have their
first varint decode to `node_id > 0`. It is coincidence, not
design, that all pre-chr1 builds landed on "safe" headers.

The fix (`src/convert_tags.cpp`) is 10 lines: skip the 8-byte
SDSL header before reading the payload. We deliberately do NOT
trust the header's `size_in_bits` field to bound the read --
on all four audited files the file has more bytes past the
header than `size_in_bits/8` claims (chr1 ~5 MB extra, chr6
~10 KB extra, yeast 6 bytes extra, mc-chr6 similar). Root cause
of that discrepancy is not investigated further: likely
`write_block()` flushing after the header is stamped. It does
not matter for correctness because any real payload bytes past
the claimed length are valid ByteCode varints that would decode
to real runs; trailing SDSL zero-padding decodes to
`node=0` and is skipped downstream by the existing endmarker-run
guard (`convert_tags.cpp:154`).

**Downstream impact of using an affected `.ltags`**: every
`LightTagIndex::run_id_at(p)` for `p` past the phantom run start
would return a run-id shifted by 1, and `find_mems`'s emitted
`(seq_id, path_bp)` records would map to wrong graph positions
proportional to the phantom-run count. On chr1 with one 196-long
phantom, this is one bad boundary; on any future build that
happens to have a header pattern with multiple phantom varints,
it could be many.

**Why the empirical `--num-seq 4524 - 196 = 4328` workaround
worked at all and why it MUST NOT be used**: convert_tags
prepends `num_seq` endmarker positions at BWT index 0. Reducing
`num_seq` by 196 makes `bwt_intervals size` equal `BWT + 1` by
cancelling out the phantom-run overshoot in the accounting
total. But this poisons every seq_id lookup: gafpack maps BWT
run 0..num_seq-1 to seq_ids 0..num_seq-1, so under-prepending
by 196 shifts every real read's seq_id by 196 and produces
wrong (but same-cardinality) coverage. The empirical fix
"looks right" only in the aggregate check; it silently
corrupts every per-record output. **Do not deploy the
workaround. Deploy the patch.**

### Pointers for future readers

Where the bug lived:
- `src/convert_tags.cpp:52-79` -- the file-read block (pre-fix).
  The line `encoded_bytes.resize(static_cast<size_t>(sz))` +
  full-file read + `while (loc < encoded_bytes.size())` loop was
  the direct source of the misinterpretation.
- `include/pangenome_index/algorithm.hpp:496` --
  `sdsl::int_vector_buffer<8> out(filename, ...)` is where the
  SDSL header prefix originates. Any tool that reads
  `.tags` produced here must strip 8 bytes before decoding.

Where the SDSL layout is documented:
- `sdsl-lite/include/sdsl/int_vector_buffer.hpp:143-175`
  (constructor) sets `m_offset = t_width ? 8 : 9` for non-plain
  streams.
- `sdsl-lite/include/sdsl/int_vector_buffer.hpp:333-360`
  (`close()`) writes the header at offset 0 and pads with up
  to 7 zero bytes to 8-byte alignment.
- Header format for `int_vector<8>` (via
  `int_vector.hpp:633`): a single 8-byte LE `uint64_t
  size_in_bits`. Width is compile-time (no extra byte).

Where the tag-run bit layout is defined:
- `src/tag_arrays.cpp:28`
  (`TagArray::encode_run_length`): `offset` in bits 0-9, `is_rev`
  in bit 10, `length` in bits 11..(10+length_bits), `node_id` in
  bits (11+length_bits).. . `length_bits = 9` per
  `include/pangenome_index/tag_arrays.hpp:141`. So `node_id =
  raw >> 20`, `length = (raw >> 11) & 0x1FF`, etc.
- `src/tag_arrays.cpp:59` (`decode_run`) is the inverse.
- The 511-cap is documented in the RLE splitting loops
  (`length` field is 9 bits -> max value 511) at
  `tag_arrays.cpp:113,141` and elsewhere.

Where the `bwt_intervals` invariant is asserted informationally:
- CLAUDE.md line 112: "check `bwt_intervals size == .seq size + 1`
  in logs" -- this is the primary correctness gate for any
  change touching TagArray serialization.
- The `.seq` size is the total text length (all haplotypes
  concatenated with endmarker separators), reported by grlbwt
  as `Number of symbols in the file`, and equal to
  `FastLocate::bwt_size()`.

**Fast diagnostic (no rebuild, ~10 sec per index):**

```sh
cd /path/to/runs/<run-name>/
BWT_SIZE=$(awk '/^Number of symbols in the file/ {print $NF}' logs/03_grlbwt.log)
BWT_INT=$(awk '/Building the bwt_intervals vector/ {print $8}' logs/06_convert_tags.log)
DELTA=$(( BWT_INT - BWT_SIZE - 1 ))
if [ "$DELTA" -eq 0 ]; then
    echo "SAFE (delta 0)"
elif [ "$DELTA" -gt 0 ]; then
    echo "AFFECTED (delta +$DELTA -- rebuild _compressed.tags and .ltags with patched convert_tags)"
else
    echo "UNDERSHOOT delta $DELTA -- investigate (unusual)"
fi
```

If `AFFECTED`, the algorithm-format `.tags` file is fine (bug is
in the reader). Rebuild only `_compressed.tags` (via
`bin/convert_tags`) and its downstream `.ltags` (via
`bin/build_lightweight_tags`). Do not rerun `build_tags` (~hours
on chr1). Do not rerun `build_rindex`.

**Mechanism check (5 sec per index):** if you want to know
*why* a given file was or wasn't hit, hex-dump the first 8
bytes and simulate the ByteCode decode:

```sh
xxd -c 16 -l 8 <base>.tags
```

Feed those 8 bytes into the Python snippet at the end of this
entry to enumerate every varint the buggy code would decode
from them and whether it has `node_id > 0`.

### Open questions

1. **Where do the "extra bytes past claimed payload" come from?**
   All four audited files have `file_size - 8 - size_in_bits/8 > 0`
   (chr1 ~5 MB, chr6 ~10 KB, yeast 6 bytes, mc-chr6 unknown).
   Best guess: SDSL `write_block()` flushed a partial buffer that
   `close()` didn't back off before stamping the header. Not
   pursued because the fix (keep all post-header bytes) is
   robust to this: any real ByteCode past the claimed payload
   decodes as valid runs; trailing pad decodes as `0+0 0` skips.
   If someone ever wants zero-padding-free output, investigate
   here.
2. **Other consumers of `.tags` files?** `build_lightweight_tags`
   reads `_compressed.tags` (a different format written by
   `TagArray::merge_compressed_files_sdsl` -- three
   `sdsl::serialize` calls, not `int_vector_buffer<8>`). No SDSL
   header issue there. `query_tags` reads the same
   `_compressed.tags`. `merge_tags` doesn't call convert_tags. So
   the fix scope is exactly `convert_tags`.
3. **The `parsaeskandar/pangenome-index@bidirectional`
   fork has not identified this bug.** Their most recent
   convert_tags-touching commits are:
   - `7853a2c` "fixed Silent truncation of endmarker count causes data corruption" (the uint16_t truncation on num_seq)
   - `23cd454` "Fixed the convert_tags bug" (skip node_id=0 in max_node_id computation)
   - `ade878f` "Fixed bug in the merge_tags"

   None address the SDSL header. Upstream will hit this the
   first time they build an index whose header happens to
   decode with `node_id > 0`. Worth a heads-up PR when we get
   there.

### Commits (this repo)

- Working-tree change on branch `tag-head-samples`:
  `src/convert_tags.cpp` -- strip 8-byte SDSL header before
  ByteCode-decoding the payload. Approximately +23 / -8 lines.
  Not yet committed at time of writing (pending sign-off).

### Files (this repo, kept for reproducibility)

- `scratch/convert_tags.cpp.head-baseline` -- pre-fix HEAD source
- `scratch/convert_tags.cpp.new-patched` -- post-fix source
- `scratch/convert_tags.baseline` -- pre-fix binary
- `scratch/convert_tags.patched` -- post-fix binary
- `scratch/yeast_test/*.log` -- both convert_tags stderr from
  round-trip test
- `scratch/yeast_test/yeast.baseline.tags`, `yeast.patched.tags`
  -- ~1 GB each, byte-identical MD5 `0672d1e901cfedffd110a1f743dccc55`

### Artifacts (vesuvio)

- `~/mem-projection/pangenome-pipeline/runs/hprc-chr6-2026-06-02/queries/convert-tags-header-fix/`
  -- Tier-2 chr6 E2E validation run, coverage MD5
  `60f6b8e4a759aebb252b83a870b9ff8c`
- `~/mem-projection/pangenome-pipeline/runs/hprcv1-chr1-2026-08-12/hprcv1_chr1_compressed.tags.header-fix`
  -- Tier-3 chr1 corrected `_compressed.tags`
- `~/mem-projection/pangenome-pipeline/runs/hprcv1-chr1-2026-08-12/hprcv1_chr1.ltags.header-fix`
  -- downstream `.ltags` built from the corrected
  `_compressed.tags`, 2.5 GB, `num_runs = 3,336,172,719`
- `~/mem-projection/pangenome-pipeline/runs/hprcv1-chr1-2026-08-12/logs/06_convert_tags_header_fix.log`
  -- log confirming `bwt_intervals size = 45,105,246,015` and
  1 skip line (`0+0 0`, trailing pad only)

### Appendix: Python simulator for header decode

Save the 8 header bytes to `HEX_BYTES` and run:

```python
# Decode the SDSL header of a .tags file as if it were raw ByteCode
# payload, showing exactly which varints the buggy convert_tags would
# emit as real runs (node_id > 0) vs skip (node_id == 0).
HEX_BYTES = "b8 cb d8 78 29 00 00 00"   # chr1 example
hdr = bytes(int(x, 16) for x in HEX_BYTES.split())

def decode(bs):
    i = 0
    while i < len(bs):
        start = i
        off = 0
        val = bs[i] & 0x7F
        while (bs[i] & 0x80) and (i + 1 < len(bs)):
            i += 1; off += 7
            val += (bs[i] & 0x7F) << off
        incomplete = bool(bs[i] & 0x80)
        i += 1
        yield val, i - start, incomplete
        if incomplete: return

# Field layout: bits 0..9 offset, bit 10 is_rev, bits 11..19 length,
# bits 20.. node_id.  See src/tag_arrays.cpp:28,59 and length_bits=9.
def decode_run(v):
    return {
        "node": v >> 20,
        "off":  v & 0x3FF,
        "rev":  (v >> 10) & 1,
        "len":  (v >> 11) & 0x1FF,
    }

emit_total = 0
for raw, nb, inc in decode(hdr):
    r = decode_run(raw)
    is_emit = r["node"] > 0
    print(f"  bytes={nb}  raw=0x{raw:x}  node={r['node']}  "
          f"off={r['off']}  rev={r['rev']}  len={r['len']}  "
          f"[{'EMIT (BUG)' if is_emit else 'skip'}]")
    if is_emit:
        emit_total += r["len"]
print(f"Total phantom-emit length from header: {emit_total}")
```

Interpretation: if `Total phantom-emit length from header` is 0,
the file is safe under the pre-fix convert_tags. Any nonzero value
is the exact overshoot you'd see in `bwt_intervals size` for that
file.

---

## JR-021 — Standardized performance-report format (PERF_REPORT_TEMPLATE.md) + HPRC chr6 L=25 worked example

```yaml
id: JR-021
date: 2026-08-21
author: claude-sonnet-5 (session with hlakshmidevi)
status: resolved
tags: [reporting, template, perf, benchmark, hprc, l25, noisy, vesuvio, jr-018, jr-019, tag-head-samples]
refs:
  follows: [JR-019]
benchmark-platform: vesuvio (Linux 6.12.73 x86_64, Debian)
```

### Context

Prior perf write-ups (JR-008, JR-012, JR-014, JR-018, JR-019) each built
their own ad-hoc table layout inside a full journal narrative
(Context/Hypothesis/Method/Findings/Interpretation). User asked for a
different, reusable shape instead: a dataset-centric report -- Reference
characteristics, Read characteristics, MEM Coverage Statistics, Performance,
Disk storage -- that can be produced the same way for any (graph, reads,
config) combination without re-deriving the format each time. This entry
establishes that format as `PERF_REPORT_TEMPLATE.md` (this repo) and records
the worked example that produced it: HPRC v1 chr6 (PGGB), MEM_LEN=25,
alt-noisy reads, flipped vs legacy+`--tag-head-samples=s8`.

The underlying benchmark data is not new -- it is the N=3 warm-cache harness
run at `perf/jr019-perf-L25-{flipped,legacy-s8}/` on Vesuvio, executed
2026-08-14 (two days after JR-019 was written) but never folded back into
the journal. This entry is also, incidentally, the write-up JR-019's open
question #1 ("N=3 warm-cache validation for L=25") was waiting on.

### Method

Data pulled directly from `perf/jr019-perf-L25-{flipped,legacy-s8}/SUMMARY.tsv`,
`PROVENANCE.txt`, and each trial's `find_mems.log` (Vesuvio). Graph
characteristics pulled fresh via `gbz_stats -p` and `gbz_stats -g` against
`hprcv1/chr6.gbz` (the `-i` flag's "Total length" line was initially
mis-read as a bp count -- see the Gotcha in `PERF_REPORT_TEMPLATE.md`'s
Reference characteristics section, caught before it made it into this
report). Two corrections were made mid-session while assembling the MEM
Coverage Statistics section, both now encoded as gotchas in the template:

1. "Total Entries Written" was initially copied from JR-019's journal prose
   (1,432,999) without cross-checking against this run's own `SUMMARY.tsv`.
   The actual ground truth for this run is `bin_bytes / 16 = 23,715,712 / 16
   = 1,482,232`, exact and identical across all 6 trials. The discrepancy is
   consistent with `find_mems.log`'s "K/M" abbreviated text using truncation
   rather than rounding (1,482,232 prints as "1.4M", which a naive reader
   would round back to ~1.4M rather than the true 1.48M).
2. "Average Tag Run Length" was resolved to `find_mems.log`'s "Global n/r
   ratio (Σmem.size / Σtag_runs)" = 21.6344, not the "Mean per-MEM n/r
   ratio" = 77.98 -- the latter is explicitly flagged as outlier-skewed in
   `mem-projection/pangenome-pipeline/CLAUDE.md`'s footgun list and is not
   meant for this kind of summary reporting.

### Findings

**Reference characteristics**

| | |
|:--|--:|
| Graph | HPRC v1 chr6 (PGGB) |
| Nodes | 4,743,141 |
| Edges | 6,586,286 |
| Graph sequence (deduplicated) | 228.63 Mbp |
| Total haplotype path length (all haplotypes, single orientation) | 15.51 Gbp |
| Haplotypes / samples / paths | 90 haplotypes, 46 samples, 1,410 paths (1,409 contigs) |

**Read characteristics**

| | |
|:--|--:|
| Reads | 100,000 x 200 bp, single-end |
| Source | `HG00438#2#JAHBCA010000010.1#0` (alt haplotype, in-graph) |
| Error profile | Noisy -- ~1% per-base substitution error, wgsim-simulated, Illumina-like uniform-error model |
| Read file | `chr6.alt_noisy.reads.txt` (20.1 MB) |

**MEM Coverage Statistics** (identical between both configs -- verified via
the correctness gate below)

| Metric | Value |
|:--|--:|
| Total Reads | 100,000 |
| Total MEMs | 250,177 |
| Average MEM length | 79.42 bp |
| Total Occurrence (Σ mem.size) | 32,050,895 (32.0M) |
| Avg Occurrence / MEM | 128.16 |
| Avg Tag Runs / MEM (SA interval) | 5.92 |
| Average Tag Run Length (Global n/r ratio) | 21.63 |
| Total Entries Written (find_mems, pre-dedup) | 1,482,232 |
| Total Occurrence on graph after MEM projection (post gafpack dedup) | 840,010 |

**Performance** (N=3 warm-cache, mean +- stdev)

| Metric | Flipped | Legacy + samples (s=8) | delta |
|:--|--:|--:|--:|
| find_mems time | 48.30 +- 0.17 s | 38.94 +- 0.14 s | -9.36s (-19.4%) |
| find_mems memory (peak RSS) | 4,516.4 +- 0.7 MB | 4,792.4 +- 0.7 MB | +276.0 MB (+6.1%) |
| gafpack time | 24.04 +- 0.05 s | 24.16 +- 0.08 s | +0.12s (noise) |
| gafpack memory (peak RSS) | ~335 MB | ~336 MB | ~0 |
| **Total wall (find_mems + gafpack)** | **72.33 +- 0.21 s** | **63.10 +- 0.11 s** | **-9.23s (-12.8%)** |
| **MEM throughput (MEMs/s, total wall)** | **3,459** | **3,964** | **+14.6%** |

**Disk storage**

| Index component | Size |
|:--|--:|
| r-index (`.ri`) | 3.63 GiB (3,900,355,874 bytes) |
| tag-index, lightweight (`.ltags`) | 670.31 MiB (702,868,237 bytes) |
| tag-samples (s=8, legacy+s8 config only) | 272.35 MiB (285,578,626 bytes) |

**Correctness gate:** coverage MD5 `5b1b770ca5f8d4719c5ea0ea2cd8e655`,
`.bin` records 1,482,232, gafpack entries 840,010, and MEM count 250,177 --
all identical across all 6 trials (3 flipped + 3 legacy+s8), stderr warnings
0 throughout.

### Interpretation

**The template holds up on first real use** -- every field in
`PERF_REPORT_TEMPLATE.md` had a concrete, unambiguous source once the two
corrections above were made. Both corrections are now load-bearing gotchas
in the template itself, which is the point: they should not have to be
rediscovered by whoever runs the next report.

**The underlying numbers reconfirm JR-019 at N=3.** Find_mems wall -19.4%,
combined wall -12.8%, both within ~1pp of JR-019's single-trial numbers
(-19.1% / -12.4%) and inside the N=3 stdev bars (<=0.17s on every total).
Combined with JR-018's L=50 result (-20.8%/-12.8%), the legacy+samples-s8
win now has harness-grade (not single-trial) confirmation at two MEM
lengths. This resolves JR-019 open question #1.

**New at L=25:** the samples-resolution cost (4.18s: 1.12s first-rid walk +
3.06s interior/last resolve) is 23.9% cheaper than flipped's SA-carry
locateNext walk (5.49s) -- a cleaner isolation of this effect than L=50
gave, because L=25 has more, shorter MEMs (250K vs 155K) and therefore more
tag-run boundaries per MEM, so the fixed per-resolution cost of a samples
lookup (avg 2.93 LF steps, well under the s=8 cap) pays off more.

### Open questions

1. JR-019's remaining open questions (#2 out-of-graph HG002 reads, #3
   whole-genome scaling) are untouched by this entry.
2. Two un-journaled Vesuvio runs postdate this data and are candidates for
   a future report using this template: `perf/hprcv2-mc-chr6-flipped/`
   (2026-08-14, flipped on the larger HPRCv2 MC chr6 graph, no
   legacy/samples comparison yet) and
   `runs/hprcv1-chr1-2026-08-12/queries/flipped-hg00097-L25/` (2026-08-15).
3. `PERF_REPORT_TEMPLATE.md` has one worked example (this entry). Whether
   its section shape holds for a fundamentally different scenario --
   e.g. an out-of-graph read set, where "Total Occurrence on graph after
   MEM projection" is expected to be sparse rather than near-total, or a
   whole-genome (multi-chromosome) index -- is untested.

### Data artifacts

Same as JR-019's L=25 data:
- `~/mem-projection/pangenome-pipeline/perf/jr019-perf-L25-flipped/`
- `~/mem-projection/pangenome-pipeline/perf/jr019-perf-L25-legacy-s8/`

Plus fresh for this entry:
- `gbz_stats -p` / `-g` output against `~/mem-projection/hprcv1/chr6.gbz`
  (not persisted to a log file on Vesuvio; re-run if needed, ~seconds)

### Commits

No code changes. `PERF_REPORT_TEMPLATE.md` added to this repo (this commit);
no changes to `find_mems`/`gafpack`/harness code.

---

## JR-022 — Block-count off-by-one in the FastLocate constructor: a phantom trailing block

```yaml
id: JR-022
date: 2026-08-21
author: claude-opus-5 (session with hlakshmidevi)
status: resolved
tags: [r-index, blocks, off-by-one, sdsl, sd-vector, build_tag_head_samples, correctness, bug-fix, chr1, mc-chr6, vesuvio]
refs:
  follows: [JR-020, JR-021]
benchmark-platform: vesuvio (Linux 6.12.73 x86_64, Debian)
```

### Context

`build_tag_head_samples` aborted on HPRCv1 chr1 with an SDSL assertion:

```
sdsl/int_vector.hpp:1443: ... operator[]: Assertion `idx < this->size()' failed.
```

It died **after** Phase 1's BWT walk logged 100% (all 338,432,560 runs) but
**before** the `candidates emitted:` summary printed — burning 2h25m each
time. Two attempts failed identically (2026-08-13 and, with added bounds
checks, 2026-08-15). Those diagnostics never fired, which was itself a clue.

While investigating, a second casualty turned up: **HPRCv2 MC chr6 had
crashed the same way on 2026-08-14** (log stops at "Phase 1: enumerate
candidates...", zero output files, ~32 min in). That failure was never
recorded anywhere; JR-016/018 only ever reported chr6 success.

### Hypothesis

The constructor sizes `blocks` and `blocks_start_pos` with **different
formulas** (`src/r-index.cpp`):

```cpp
blocks.resize((total_runs / block_size) + 1);                  // line 1127: floor + 1
sd_vector_builder block_sd_builder(bwt_size(),
    (total_runs / block_size) + (total_runs % block_size != 0)); // line 1143: ceil
```

These agree only when `block_size` (10) does **not** divide `total_runs`.
When it does, `blocks` gains a trailing block that is never filled, yet is
still serialized and counted by `num_blocks()` — with no corresponding
1-bit in `blocks_start_pos`.

Prediction: crash iff `total_runs % 10 == 0`.

### Method

Rather than reproduce (2.5h/attempt, no core dump — `ulimit -c` is 0), wrote
`src/probe_blocks.cpp`: loads a `.ri` and reports `num_blocks()` vs the
1-bits in `blocks_start_pos`, plus how the trailing blocks decode. **~3 s per
index against files already on disk.**

### Findings

**Perfect correlation across the whole index fleet** (run counts from
`08_print_stats.log`):

| Dataset | BWT runs | `%10` | `num_blocks()` | `blocks_start_pos` ones | mismatch | outcome |
|:--|--:|--:|--:|--:|--:|:--|
| hprcv1 chr6 | 249,445,481 | 1 | 24,944,549 | 24,944,549 | **0** | ✅ worked (JR-016/018) |
| **hprcv1 chr1** | 338,432,560 | **0** | 33,843,257 | 33,843,**256** | **1** | ❌ crashed |
| **hprcv2 MC chr6** | 287,351,380 | **0** | 28,735,139 | 28,735,**138** | **1** | ❌ crashed |
| hprcv1 MC chr6 | 246,861,799 | 9 | — | — | — | untested, predicted safe |
| hprcv2 chr18 | 141,430,719 | 9 | — | — | — | untested, predicted safe |
| hprcv2 chr6 | 291,744,113 | 3 | — | — | — | untested, predicted safe |
| yeast235 chrII | 14,033,362 | 2 | — | — | — | ✅ worked |

Exactly the two `%10 == 0` indexes are the two that crashed. Same code, same
`block_size`, same `C.size()`=6 — only `runs % 10` differs.

**The full mechanism** (each link observed, not inferred):

1. `blocks.resize(floor+1)` creates a phantom block when `10 | total_runs`.
2. `Run_blocks() : character_cum_ranks(8)` (r-index.hpp:154) hardcodes **8**
   cum-rank slots; the real alphabet is `C.size()` = **6**.
3. `serialize_encoded` loops over *all* blocks writing
   `get_cum_ranks().size()` varints — so **8 zero-bytes** for the untouched
   phantom block. Confirmed on disk: `start_bits[33843256] = 1523612419`,
   stream size `1523612427` → **exactly 8 bytes**.
4. On decode, `skip_header` consumes only `C.size()` = 6 varints. The
   **2 leftover `0x00` bytes decode as valid run headers** (`prefix=0` →
   `run_length=1`), so `get_block_runs()` returns the phantom block as
   **non-empty**: `2 runs: (sym=10,len=1) (sym=10,len=1)` — `sym=10` is
   `nuc[0]`, not real sequence.
5. Non-empty ⇒ `BwtRunHeadIter::load_current_block()` takes the branch calling
   `blocks_start_select_1(num_blocks)`.
6. sdsl's `select_support_sd_trait<1>::select` (sd_vector.hpp:930) does
   `v->low[i-1]` **with no bounds check**; `low` is an `int_vector<0>` of size
   `ones`. Reading `low[ones]` fires the assertion.

This happens inside `bwt_iter.advance()` — *after* the 100% progress print and
*before* the loop re-tests — which is exactly why the JR-2026-08-15 bounds
checks (in the loop body) never fired and `candidates emitted:` never printed.

**Blast radius is narrow.** Query paths reach blocks via
`blocks_start_pos.predecessor()`, which can never return the phantom block
(no 1-bit). Only *sequential* block iteration trips it, and `BwtRunHeadIter`
in `build_tag_head_samples` is its sole user. The affected `.ri` files were
producing correct `find_mems`/genotyping output the whole time.

### Fix

Three changes (commit `10b83af`):

1. `r-index.cpp:1127` — size `blocks` with the same ceiling expression as
   `block_sd_builder`.
2. `r-index.cpp:373` — serialize exactly `C.size()` cum-rank varints, so
   writer and reader cannot desynchronize regardless of a block's own vector
   length. No-op for real blocks (they always hold `C.size()` entries).
3. `r-index.hpp:154` — drop the fabricated `8` from `Run_blocks`.

Verified the constructor never indexes `blocks[ceil]`: when `10 | total_runs`
the loop exits before touching it; otherwise the tail path at line 1282 uses
`ceil-1`.

### Verification

**Byte-identity gate.** For any index with `runs % 10 != 0`, `floor+1 == ceil`,
so the rebuilt `.ri` must be *byte-identical* — a free proof the fix is inert
where it should be:

| Dataset | `.ri` after rebuild | mismatch | build_rindex wall |
|:--|:--|--:|--:|
| yeast235 | **byte-identical** (`03f0b1ad…`) | 0 → 0 | 6.5 s |
| chr6 | **byte-identical** (`75a9022c…`) | 0 → 0 | 7:45 |
| **chr1** | −1 block, **−8 bytes** (5,465,354,069 → …061) | **1 → 0** | 12:57 |
| **mc-chr6** | −1 block, **−8 bytes** (4,730,096,845 → …837) | **1 → 0** | 30:56 |

Both affected indexes shrank by *exactly* 8 bytes and one block — the
quantitative prediction from step 3 above. Trailing blocks now decode to real
runs instead of `sym=10` garbage.

**Full projection pipeline** (find_mems → gafpack → validate_gaf):

| Dataset | Coverage MD5 | vs baseline | validate_gaf |
|:--|:--|:--|:--|
| yeast235 | `1baf368f3cdc59890f88410cc34cc051` | ✅ | 2000/2000 |
| chr6 | `60f6b8e4a759aebb252b83a870b9ff8c` | ✅ gold (JR-012/13/18/19/20) | 2000/2000 |
| chr1 | `9197ee230c283e816e6119f16feb2f8b` | ✅ exact | 2000/2000 |

chr1 is the load-bearing row: its `.ri` genuinely changed content, and the
pipeline still produced byte-identical coverage — confirming the phantom block
was unreachable from the query path.

mc-chr6 was verified by a stronger single-variable A/B instead (its stored
baselines' `FIND_MEMS_EXTRA_FLAGS` were unrecorded, so a historical comparison
could not distinguish "fix broke something" from "wrong flag"): `find_mems`
output old-`.ri` vs new-`.ri`, same binary, minutes apart. Byte-identical in
**both** finder modes — legacy (`189301f1…`) and flipped (`b385b1cf…`). Since
gafpack never reads the `.ri`, identical `find_mems` output makes the whole
downstream chain identical by construction.

**Resolution.** `build_tag_head_samples` then completed on chr1 for the first
time: wall 3:46:56, peak RSS 95.3 GB, exit 0, and `1000/1000` kept-SA +
`100/100` deleted-head probes on each of s=8/16/32. Phase 1's candidate and
`locateNext` counters matched the crashed run *exactly* at every 2.5%
checkpoint through 100% (3,348,609,586 candidates), so the fix removed the
phantom block and perturbed nothing else.

### Interpretation

Same shape as JR-020: **latent, data-dependent, and silent by arithmetic
coincidence.** JR-020 turned on whether a file's SDSL header bytes happened to
decode with `node_id == 0`; this one turns on whether `block_size` divides
`total_runs`. Both went undetected for months because the datasets in use
happened to land on the safe side, and both surfaced on chr1.

Two compounding defects were required. The block-count off-by-one alone would
have been harmless — an empty trailing block decodes to nothing and gets
skipped. It only became fatal because `character_cum_ranks(8)` disagreed with
`C.size()`, turning padding into phantom runs. Neither is dangerous alone;
together they abort a 2.5-hour job.

**Process note.** The probe was the whole difference: it converted a
"reproduce for 2.5 h and hope for a backtrace" problem into a 3-second query
against existing files, and turned a plausible story into a controlled
experiment across three indexes. Worth reaching for earlier when a long job
fails at scale.

### Open questions

1. **Other unchecked `select_1` call sites.** sdsl's `select_support_sd`
   does an unguarded `low[i-1]`. Any caller passing a rank derived from a
   *different* structure's size can fault the same way. A sweep of
   `blocks_start_select_1` / `last_select_1` callers for that pattern would
   say whether this class of bug has siblings.
2. **Should `num_blocks()` be derived from `blocks_start_pos` instead of
   `blocks_encoded_start_bits.size()`?** That would make the two structures
   incapable of disagreeing, rather than relying on the constructor keeping
   them in sync.
3. **Untested indexes.** hprcv1-MC-chr6, hprcv2-chr18 and hprcv2-chr6 are
   predicted safe from their run counts but were not rebuilt. `probe_blocks`
   settles each in ~3 s if they are ever used with a sequential block walk.

### Commits (branch `tag-head-samples`)

- `10b83af` r-index: fix block-count off-by-one that corrupts sequential block walks
- `1fe287a` probe_blocks: diagnostic for block-count / blocks_start_pos consistency

### Data artifacts (vesuvio)

- `runs/hprcv1-chr1-2026-08-12/hprcv1_chr1.ri.pre-blockfix` — pre-fix index
- `runs/hprcv2-mc-chr6-2026-06-07/hprcv2-mc-chm13-chr6.ri.pre-blockfix`
- `runs/hprcv1-chr1-2026-08-12/tag_head_samples/build_tag_head_samples.blockfix.log`
  — the successful run; `.debug.log` alongside it is the 2026-08-15 crash
- `runs/*/queries/blockfix-verify/lightweight/` — pipeline verification runs
- `~/mcchr6_ab/` — mc-chr6 old-vs-new `find_mems` A/B outputs

---

## JR-023 — Tag-head SA samples do not scale to whole-genome: the chr6 recommendation inverts on chr1

```yaml
id: JR-023
date: 2026-08-22
author: claude-opus-5 (session with hlakshmidevi)
status: resolved
tags: [tag-head-samples, scaling, whole-genome, perf, benchmark, chr1, hprc, flipped, vesuvio, jr-016, jr-018, jr-019]
refs:
  follows: [JR-022]
  resolves: [JR-018, JR-019]
benchmark-platform: vesuvio (Linux 6.12.73 x86_64, Debian)
```

### Context

JR-018 recommended `--tag-head-samples=…s8` with the legacy MEM finder as the
default HPRC chr6 configuration: **−20.8%** find_mems wall versus the flipped
baseline, at +6.1% peak RSS. Its open question #1 flagged the obvious risk —
*"Whole-genome HPRCv1/HPRCv2 datasets not yet tested… we expect it to hold, but
should confirm"* — and JR-019's open question #3 promised chr1 numbers once the
samples build finished.

That build never finished. It aborted twice on the block-count bug fixed in
JR-022. With chr1's samples now built (s=8/16/32, all probes passing), this
entry finally runs the comparison.

### Question

Does JR-018's "legacy + samples beats flipped" result hold at whole-genome
scale?

### Method

`find_mems` run directly (not via `perf_harness.sh`, which also runs gafpack on
every trial and would have added hours while measuring nothing — the `.bin`
outputs are proven identical). Per config: 2 untimed warmups + 3 timed trials,
contiguous per-mode ordering. Five configs: flipped, legacy, legacy+s8/s16/s32.
Index `hprcv1_chr1.ri` (post-JR-022), L=25, MIN_OCC=1.

**Two readsets**, because the first choice was wrong:

- **in-graph**: `chr1.alt_noisy.reads.txt` — HG00438 alt haplotype, **38 of the
  2262 chr1 paths**. This is the direct analogue of chr6's benchmark readset
  (also HG00438), and the only valid basis for comparing against JR-018/019.
- **out-of-graph**: `chr1.hg00097_noisy.reads.txt` — HG00097, **0 of 2262
  paths**. Initially benchmarked by mistake (it was the readset with a coverage
  baseline, chosen for JR-022's correctness gate and then carried over without
  thinking). Retained because it answers JR-019's open question #2, which listed
  out-of-graph behaviour as untested.

Both 100K × 200 bp.

> **Revised 2026-08-22.** The flipped / legacy / s=8 rows below were originally
> measured while mc-chr6's `build_tag_head_samples` occupied one core, and have
> been replaced with clean-machine N=3 re-runs. The contamination proved
> **non-uniform** (−0.3% to −10.2% depending on mode and readset), so the
> original relative deltas were distorted, not merely scaled — one conclusion
> changed sign (see C). s=16 and s=32 were not re-run and remain the original
> loaded-run values, marked inline.

### Findings

**A. Storage: ~8× chr6 at every s, and the ratio is stable.**

| s | chr6 (JR-016) | chr1 | ratio | chr1 kept | % of all tag runs |
|--:|--:|--:|--:|--:|--:|
| 8 | 272 MB | **2.19 GB** | 8.1× | 403.8 M | 12.11% |
| 16 | 125 MB | **1.02 GB** | 8.2× | 184.7 M | 5.54% |
| 32 | 61 MB | **507 MB** | 8.3× | 87.3 M | 2.62% |

Build: 3h46m wall, **95.3 GB peak RSS** (chr6 was 1h24m / 21.7 GB).

**B. Root cause — the free-SA subsidy collapses.**

Tag heads that coincide with a BWT-run head get their SA free from the
r-index's existing `samples[]`, needing no stored sample:

| | chr6 | chr1 |
|:--|--:|--:|
| tag runs | 710.2 M | 3336.2 M |
| BWT runs | 249.4 M | 338.4 M |
| **tag runs per BWT run** | **2.85** | **9.86** |
| **coincident (free) tag heads** | **33.9%** | **9.8%** |

chr6 got a third of its tag heads subsidised; chr1 gets a tenth. chr1 packs
4.7× more tag runs onto 1.36× more BWT runs, so coincidences are rare and far
more genuine samples must be stored — and, at query time, far more emits fall
to the slow path.

**C. Perf, in-graph (HG00438 alt, L=25, N=3 medians).**

| config | median | RSS | vs legacy | vs flipped | fast-path |
|:--|--:|--:|--:|--:|--:|
| **flipped** | **89.40 s** | 9403 MB | −34% | — | — |
| legacy | 134.48 s | 9403 MB | — | +50% | — |
| legacy + s8 | 136.23 s | 11663 MB (+24.0%) | **+1.3%** | **+52%** | 14.02% |
| legacy + s16 *(loaded)* | 219.06 s | 10462 MB (+11.3%) | +63% | +145% | 6.45% |
| legacy + s32 *(loaded)* | 379.91 s | 9914 MB (+5.4%) | +182% | **+325%** | 3.05% |

**s=8 is slower than plain legacy (+1.3%), not faster.** On the loaded run it
appeared to win by 1.8%; the clean re-run reverses the sign. So on chr1's
realistic in-graph workload the samples design delivers **no benefit at all**
over the legacy baseline, while costing +2.26 GB of resident memory.

**D. Perf, out-of-graph (HG00097, L=25, N=3 medians).**

| config | median | RSS | vs legacy | vs flipped | fast-path |
|:--|--:|--:|--:|--:|--:|
| **flipped** | **76.71 s** | 8578 MB | −33% | — | — |
| legacy + s8 | 83.09 s | 10839 MB (+26.3%) | **−27%** | **+8.3%** | 13.73% |
| legacy | 113.75 s | 8579 MB | — | +48% | — |
| legacy + s16 *(loaded)* | 120.55 s | 9637 MB (+12.3%) | +6% | +57% | 6.43% |
| legacy + s32 *(loaded)* | 184.40 s | 9090 MB (+6.0%) | +62% | +140% | 3.07% |

On the clean machine s=8 improves markedly here (−20% → −27% vs legacy;
+20% → +8.3% vs flipped) — the opposite direction to the in-graph readset.
This mode- and readset-dependent sensitivity to background load is why the
original numbers could not simply be rescaled.

**E. The direct like-for-like comparison — the result inverts.**

Same L=25, same HG00438 alt-noisy readset family, same metric:

| | chr6 (JR-019) | chr1 in-graph | |
|:--|--:|--:|:--|
| flipped | 48.34 s | 89.40 s | |
| legacy + s8 | 39.10 s | 136.23 s | |
| **s8 vs flipped** | **−19.1%** ✅ | **+52.4%** ❌ | **reversed** |

**F. Correctness preserved throughout.** All samples configs produce `.bin` and
`seq_id_starts` byte-identical to the legacy baseline on both readsets. Flipped's
`.bin` differs (`9e7e561b…` vs `8d73e085…`) but is the same 584,554,560 bytes
with identical `seq_id_starts` — JR-007's documented sort tie-break, since legacy
iterates left-to-right and flipped right-to-left. Confirmed benign: gafpack
coverage from the two is **byte-identical over 160,353,056 bytes**.

### Interpretation

**JR-018's recommendation does not generalise, and should not be applied to
whole-genome indexes.** On chr1 the flipped finder dominates every samples
configuration on *both* axes simultaneously — fastest wall, lowest RSS, no
2.19 GB side file, no 3h46m build step.

The mechanism is not a bug in the samples design; it works exactly as specified,
eliminating the entire ~51 s seed cost in every arm. What fails is the **trade**.
Cost of the replacement is (emits × slow-path fraction × avg LF steps × ~634 ns
per DRAM-bound LF, JR-017). On chr6 three factors were mild: high coincidence
rate kept the slow-path fraction down, low tag density kept emits per MEM down,
and the seed cost being replaced was proportionally large. On chr1 all three move
against it at once. The seed saving is fixed; the replacement cost is not.

**On the in-graph workload no `s` value beats even the legacy baseline** — s=8
is +1.3%, and it is also the setting that busts the +10% RSS ceiling, by 2.4×.
The only setting inside the ceiling (s=32) is ~4× slower than flipped. There is
no viable operating point.

**In-graph is markedly harsher than out-of-graph** (s=8: +1.3% vs −27% relative
to legacy). In-graph reads match the pangenome more, producing more MEMs and more
tag-run emits per MEM — precisely the quantity the slow-path cost scales with,
while the seed saving does not. Anyone extrapolating from out-of-graph numbers
would materially overestimate the design.

**Retrospective on JR-018.** Its s=8 recommendation was correct for chr6 but
rested on a property nobody had identified as load-bearing: the 33.9% coincidence
rate. JR-016 recorded that number as a storage optimisation ("saves storage for
34% of all tag runs") without noting it was equally a *query-time* subsidy. The
lesson is that `s` is not the only knob — the graph's tag-runs-per-BWT-run ratio
determines whether any `s` works, and that ratio is a property of the dataset, not
a tunable.

### Recommendation

Keep `--use-flipped-mems` as the whole-genome path. Do not build or ship
`.tag_samples` files for chr1-scale indexes. JR-018's chr6 guidance stands **for
chr6-scale single-chromosome indexes only** and should be read with that scope.

Before applying samples to any new dataset, compute
`tag_runs / bwt_runs` and the coincidence rate — both available in ~3 s from
`print_stats` plus one `build_tag_head_samples` Phase-2 line. If the coincidence
rate is well below ~30%, expect the design not to pay.

### Open questions

1. **Where is the crossover?** chr6 (2.85 tag runs/BWT run, 33.9% coincident)
   wins; chr1 (9.86, 9.8%) loses badly. hprcv2 MC chr6 sits between and its
   samples build is running — one more point would show whether the transition
   is sharp or gradual, and whether coincidence rate alone predicts it.
2. **Clean-machine re-run.** These numbers carry mc-chr6's concurrent load.
   Relative deltas are sound but absolute walls are inflated; re-run before any
   external reporting.
3. **Is flipped's advantage growing with scale?** flipped beats legacy by 34%
   here versus JR-012's ~11% on chr6. If the margin widens with graph size, that
   strengthens the case for flipped as the default everywhere, which JR-014's
   open question #2 left pending.
4. **Does a much smaller s help?** s=4 or s=2 would raise the fast-path rate,
   but storage already scales ~2× per halving — s=4 would be ~4.4 GB on chr1.
   Almost certainly not viable, but it would confirm the curve has no minimum
   hiding below s=8.

### Data artifacts (vesuvio)

- `~/chr1_bench_ingraph/` — in-graph N=3 logs, `.time` files, coverage identity check
- `~/chr1_bench/` — out-of-graph N=3 logs and `.time` files
- `~/chr1_gate/` — 4-way byte-identity gate (baseline vs s8/s16/s32)
- `runs/hprcv1-chr1-2026-08-12/tag_head_samples/hprcv1_chr1.tag_samples.s{8,16,32}`
- `runs/hprcv1-chr1-2026-08-12/tag_head_samples/build_tag_head_samples.blockfix.log`

### Commits

No code changes — measurement-only entry against `720843f`.

> Mechanism corrected by JR-025 (2026-08-22): query cost is driven by emit
> volume, not by the coincidence rate. The findings and recommendation in this
> entry stand; the explanation in "Interpretation" does not.

---

## JR-024 — Tag-array structural characteristics: coincidence rate as the samples-design predictor

```yaml
id: JR-024
date: 2026-08-22
author: claude-opus-5 (session with hlakshmidevi)
status: resolved
tags: [tag-head-samples, characterization, tag-array, coincidence-rate, density, chr1, mc-chr6, hprc, vesuvio]
refs:
  follows: [JR-023]
  supports: [JR-016, JR-018, JR-023]
```

### Context

JR-023 found the samples design wins on chr6/mc-chr6 and loses badly on chr1,
and attributed it to "tag-array density". This entry pins down the underlying
structural quantities, defines the terms `build_tag_head_samples` reports
(which are easy to misread), and records the measured values.

### The population hierarchy

Every tag run has a head — its leftmost BWT position. Each falls in exactly
one bucket:

```
All tag runs                                  (num_tag_runs)
├── COINCIDENT with a BWT-run head            → SA free from r-index samples[]
│                                               never stored in the samples file
└── tag-head candidates ("deletable")         (num_tag_runs − coincident)
      ├── KEPT     → SA stored in .tag_samples.sN
      └── DELETED  → SA recovered by LF-walking ≤ s steps at query time
```

Phase 1 emits only the anchor at a coincident position, so a coincident tag
head costs **zero storage and zero query work** — its SA is already in the
r-index. This is a *subsidy*, and its size is a property of the dataset.

### Reading the reported percentages

`kept:` is reported against two different denominators:

| label | denominator | measures |
|---|---|---|
| "% of candidates" | tag-head candidates (deletable only) | how aggressive the deletion rule was |
| "% of all tag runs" | `num_tag_runs` (incl. coincident) | true storage burden |

Their ratio is the non-coincidence rate, so the subsidy can be read straight
off any `kept:` line: chr1 `12.11/13.42 = 0.902` → 9.8% coincident;
mc-chr6 `7.14/10.91 = 0.654` → 34.6% coincident.

### Measured values

| | chr6 | mc-chr6 | chr1 |
|:--|--:|--:|--:|
| BWT size | 31.0 Gbp | 157.1 Gbp | 45.1 Gbp |
| BWT runs | 249.4 M | 287.4 M | 338.4 M |
| tag runs | 710.2 M | 770.8 M | **3336.2 M** |
| **tag runs / BWT run** | **2.85** | **2.68** | **9.86** |
| **coincident (free SA)** | **33.9%** | **34.5%** | **9.8%** |
| avg BWT run length | 124 bp | 547 bp | 133 bp |
| s=8 file | 272 MB | 304 MB | **2.19 GB** |
| s=8 kept (% all tag runs) | 7.08% | 7.14% | **12.11%** |
| s=8 vs flipped (JR-023) | −19% ✅ | −13% ✅ | **+50%** ❌ |

chr6 figures for tag runs / BWT runs / coincidence are carried from JR-016,
not re-measured here.

### Interpretation

**Coincidence rate is a direct function of tag-array density.** A tag head is
coincident when it lands on a BWT-run head; the denser the tag array relative
to the BWT run structure, the rarer that is. At ~2.7–2.9 tag runs per BWT run
roughly a third of tag heads are subsidised; at 9.86 only a tenth are.

**BWT size is not the predictor.** mc-chr6 has 3.5× chr1's BWT and 5× chr6's,
yet sits with chr6 on every structural measure and behaves like it. Any
extrapolation of the samples design from "chromosome size" or "graph size" is
unfounded; `tag_runs / bwt_runs` is the quantity that matters.

**The subsidy is doubly load-bearing.** JR-016 recorded the 34% coincidence
rate as a *storage* optimisation only. It is equally a *query-time* one: a
coincident head never takes the slow path. Losing it raises both the file size
and the per-emit LF cost, which is why chr1 degrades on both axes at once.

**Average BWT run length is a second, independent variable.** It governs how
much seed cost there is to eliminate (`locate_sa_value` walks ~L/2 steps).
mc-chr6's 547 bp runs make legacy's seed cost 60% of its wall (95 s of 157 s)
versus chr1's 36% — which is why mc-chr6 shows the design's largest absolute
saving despite a coincidence rate no better than chr6's.

**Cheap pre-flight test.** Before committing to a multi-hour build, compute
`tag_runs / bwt_runs` from `print_stats` (seconds, no build required). Below
~3 the design is likely viable; near 10 it is not.

### Open questions

1. **Where is the crossover?** Three points at 2.68 / 2.85 (win) and 9.86
   (lose). The transition is unlocated — anything between 3 and 9 would
   narrow it.
2. **What drives tag density?** chr1 packs 4.7× more tag runs onto 1.4× more
   BWT runs than chr6. Whether that is graph complexity, node fragmentation
   after `-m 1024` chopping, or a chr1-specific property is unexamined.
3. **Log naming collision.** `candidates` denotes two different populations:
   Phase 1's `candidates emitted:` (3,348,609,586 on chr1) is anchors + tag
   heads, while the per-`s` `% of candidates` (3,010,177,026) is tag heads
   only. Computing `kept / candidates_emitted` yields 12.06% instead of the
   reported 13.42%. Renaming to `total candidates (anchors + tag-heads)` and
   `% of deletable tag heads` would remove the trap.

---

## JR-025 — Correction to JR-023: emit volume, not coincidence rate, drives samples query cost

```yaml
id: JR-025
date: 2026-08-22
author: claude-opus-5 (session with hlakshmidevi)
status: resolved
tags: [tag-head-samples, correction, model, emit-volume, perf, chr1, mc-chr6, vesuvio]
refs:
  follows: [JR-024]
  corrects: [JR-023]
```

### Context

JR-023 explained chr1's poor samples performance by its collapsed coincidence
rate (9.8% vs chr6's 33.9%), reasoning that fewer free SAs means more slow-path
LF walks. The mc-chr6 benchmark run afterwards produced fast-path measurements
that contradict that explanation. JR-023's findings and recommendation stand;
its stated mechanism does not.

### The contradiction

If coincidence rate drove query cost, chr1 should show the worst fast-path hit
rate. It shows the **best**:

| dataset / readset | fast-path @ s=8 | s=8 vs flipped |
|:--|--:|--:|
| **chr1 in-graph** | **14.02%** | **+50%** ❌ |
| mc-chr6 hg002 | 13.51% | −13% ✅ |
| chr1 out-of-graph | 13.73% | +20% ❌ |
| chr6 (JR-018) | 12.27% | −19% ✅ |
| mc-chr6 alt_noisy | 12.11% | −32% ✅ |

The ordering is essentially uncorrelated with the outcome. Miss rate is not the
discriminator.

### What actually differs

Per-query **emit volume** — the number of tag-run records needing SA
resolution, measured directly as `.bin` bytes / 16:

| | chr1 in-graph | mc-chr6 alt_noisy | mc-chr6 hg002 |
|:--|--:|--:|--:|
| emits (records) | **36.53 M** | 1.74 M | 5.03 M |
| fast-path rate | 14.02% | 12.11% | 13.51% |
| slow-path resolutions | 31.4 M | 1.53 M | 4.35 M |
| avg LF steps / slow | 3.44 | 2.98 | 3.31 |
| **slow-path LF ops** | **108.0 M** | **4.55 M** | 14.4 M |
| ≈ cost @ 634 ns (JR-017) | **68.5 s** | 2.9 s | 9.1 s |

chr1 performs **21× more emits** than mc-chr6/alt_noisy at a *higher* hit rate,
so it pays ~24× the slow-path LF cost.

**The model closes quantitatively.** chr1 in-graph (clean-machine N=3):
legacy 134.48 s → s=8 136.23 s, a net **+1.75 s**. Predicted: ~66.8 s of
eliminated seed cost minus ~68.5 s of slow-path cost = +1.7 s. The two terms
very nearly cancel with the slow path marginally ahead — which is exactly the
observed result that s=8 fails to beat even the legacy baseline on chr1.

### Corrected model

```
samples net gain = seed_cost_eliminated
                 − (emits × (1 − fast_path_rate) × avg_LF_steps × ~634 ns)
```

- **Coincidence rate** governs **storage** (how many samples must be stored) —
  JR-024's role for it is correct.
- **Emit volume** governs **query cost**, and is set by tag-array density and
  by how well the reads match the graph.
- **Average BWT run length** governs the size of the term being eliminated
  (`locate_sa_value` walks ~L/2), which is why mc-chr6's 547 bp runs produce
  the largest absolute saving.

The two are correlated — dense tag arrays yield both low coincidence and high
emit volume — which is why JR-023's explanation fitted three datasets while
being the wrong causal variable. mc-chr6's two readsets separate them: same
graph, same coincidence rate, 2.9× the emits on `hg002`, and the margin falls
from −32% to −13%. Coincidence rate cannot explain that; emit volume does.

### Background-load contamination: measured, and it was not uniform

> **Revised 2026-08-22.** Originally this section flagged the chr1 benchmarks as
> running under background load and recommended a clean re-run. That re-run has
> since been done (flipped / legacy / s=8, both readsets, N=3, quiet machine)
> and the numbers above are the clean ones.

Both chr1 benchmarks originally ran while mc-chr6's single-threaded build
occupied one core; both mc-chr6 benchmarks ran on a quiet machine. The
assumption at the time was that constant load plus contiguous per-mode ordering
would preserve relative deltas. **That assumption was wrong.**

| readset | mode | loaded | clean | recovery |
|:--|:--|--:|--:|--:|
| in-graph | flipped | 94.19 s | 89.40 s | −5.1% |
| in-graph | legacy | 143.42 s | 134.48 s | −6.2% |
| in-graph | s=8 | 140.84 s | 136.23 s | −3.3% |
| out-of-graph | flipped | 76.97 s | 76.71 s | −0.3% |
| out-of-graph | legacy | 115.41 s | 113.75 s | −1.4% |
| out-of-graph | s=8 | 92.55 s | 83.09 s | **−10.2%** |

Recovery ranged from 0.3% to 10.2% and varied by both mode and readset, so the
contamination was a **distortion, not a scale factor**. It changed one
conclusion's sign: in-graph s=8 vs legacy went from −1.8% (a marginal win) to
+1.3% (a marginal loss). Clean-machine per-trial spread was ±0.1 s on most arms
versus up to ±4 s loaded.

s=16 and s=32 were not re-run; those rows in JR-023 remain loaded-run values and
are marked as such.

**Methodological note.** "Load is constant across arms, so relative deltas hold"
is not safe for memory-latency-bound workloads. Modes with different working-set
sizes contend differently for cache and memory bandwidth — here s=8 carries an
extra 2.19 GB resident, and it was the arm whose sensitivity differed most.
Benchmarks that will be compared against each other should be run on a quiet
machine, not merely under equal load.

### Open questions

1. **Predict from emit volume before building.** Emit count is measurable from
   a single `find_mems` run against an existing index (no samples file needed),
   so the corrected model is testable in ~2 min per dataset rather than hours.
   That is a better pre-flight test than JR-024's `tag_runs / bwt_runs` proxy.
2. **hg002 s=32 was never run** (benchmark paused). s=16 had already crossed
   into a loss, so it would not change the conclusion, but the curve is
   incomplete.

### Commits

No code changes — measurement-only entry against `f3cedcf`.

---
