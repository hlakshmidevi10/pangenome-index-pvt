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
  ../mem-projection/yeast-235/yeast-235-chrI/S288C_chrII_N100K_R1_200_reads.txt \
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
- (this entry, pending as of writing)

---
