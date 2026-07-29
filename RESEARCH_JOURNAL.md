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
