# Design: Zero-Seed-Cost MEM Finder via Flipped Iteration + R-Index SA Carrying

**Status:** Stages 0-0.85 complete (BWT primitive prerequisite implemented + verified). Stages 1-5 pending. See §5 for the implementation plan and completed prerequisites.

**Branch:** `classify-mem-runs`

**Companion artifacts:**
- `FINDINGS_SEED_COST.md` — measurement data motivating this design
- `VALIDATION_GUIDE.md` — validation protocol (100% valid on HPRC chr6 is the E2E bar)
- `CLAUDE.md` — engineering principles (correctness > performance > minimality > elegance)
- `xy-test/mem-prototype/compare_mem_algorithms.py` — Python prototype (1,014/1,014 tests establish algorithm correctness)
- `src/test_r_index_sa.cpp` — Tier-1 unit test pattern (validated `leftmost_c_in_interval` on yeast, 105K+ tests, 0 failures)

**Do not** validate against `xy-test/xy.ri` — it is a known-bad fixture. Use the canonical validation datasets from `VALIDATION_GUIDE.md`:
- Yeast (fast iteration): `../mem-projection/pangenome-pipeline/runs/v2-yeast235/yeast235_chrII_100kb_normalized.ri` with reads at `../mem-projection/yeast-235/yeast-235-chrI/S288C_chrII_N100K_R1_200_reads.txt` (config `configs/yeast235-chrII-normalized.env`).
- HPRC chr6 (production scale): `../mem-projection/pangenome-pipeline/runs/hprc-chr6-2026-06-02/hprcv1_chr6.ri` with reads at `../mem-projection/hprcv1/chr6.alt.reads.txt` (config `configs/hprcv1-chr6-alt-reads.env`).

Any subsetting (e.g., a "500 reads" quick loop) must be a deterministic prefix of the canonical reads file — pass `[max_reads=N]` to the test binaries rather than committing a separate on-disk subset.

---

## 1. Objective

Eliminate the seed cost in `dump_mem_info_lightweight` (`src/find_mems.cpp:604`).

**Measurement basis:** Instrumented `find_mems` classifier (branch `classify-mem-runs`, commit `40c7d5c`) on HPRCv1 chr6, 100,000 reads, 100,899 MEMs. Result: `seed_cost` is 32.2% of total walk cost. Seed cost concentrates in long BWT runs: 8.8% of MEMs (those with `avg_bwt_run_length >= 1000`) contribute 75.4% of all seed cost.

**Structural constraint:** MONI-style SA-carrying during backward search does not directly apply to our FMD indexer because `forward_extend_encoded` (`src/r-index.cpp:787-793`) is implemented as swap → backward-extend-on-RC → swap. The LF identity `SA[LF(k)] = SA[k] - 1` used to propagate SA values is valid only for `backward_extend_encoded` on the forward BWT, and the current three-phase MEM finder (`include/pangenome_index/algorithm.hpp:655-738`) uses `forward_extend_encoded` in Step 2. Any SA value carried through Step 1 is destroyed by Step 2.

**Approach:** Restructure MEM enumeration so it uses *only* backward extends, making the r-index SA-carrying technique applicable end-to-end. Deliver `SA[sp]` for each MEM's interval "for free" as part of the discovery walk.

**Projected impact:** If seed cost drops to zero, walk cost decreases by 32%. Applying the same downstream ratios used in `FINDINGS_SEED_COST.md` (walk ≈ 93% of MEM-processing time; MEM-processing ≈ 34% of `find_mems`): `0.32 × 0.93 × 0.34 ≈ 10.1%` end-to-end `find_mems` wall reduction.

---

## 2. Reference: current three-phase algorithm

Located at `include/pangenome_index/algorithm.hpp:655-738`, function `find_mems_function`. Outer loop `find_all_mems` at line 741.

For each right-anchor `x` (iterated left-to-right, `x_0 = 0`):

- **Step 1** (backward extend seed): from empty interval, backward-extend `P[x + min_len - 1], P[x + min_len - 2], ..., P[x]`. If size drops below `min_occ` at position `j`, return `j + 1` (no MEM emitted, next anchor advances past the failure).
- **Step 2** (forward extend to right-maximality): starting from Step 1's interval, forward-extend `P[x + min_len], P[x + min_len + 1], ...` until size drops or `j = n`. Let `e = j`. Emit MEM `[x, e - 1]` with size = last successful interval size, `bwt_start = sp`.
- **Step 3** (backward extend for next-x): from empty interval, backward-extend `P[e], P[e-1], ..., P[x+1]` until size drops at some `j`. Return `j + 1` as next anchor.

**Locate path** (`src/find_mems.cpp:604`): given emitted MEM's `[sp, ep]`, call `locate_sa_value(sp)` to obtain `SA[sp]`. This is the **seed cost** — walks `locateNext` from the BWT run's start up to `sp`, worst case `avg_bwt_run_length / 2` steps per MEM.

**Bug flag:** if Step 2 extends to `j = n`, Step 3 accesses `P[n]` (out of bounds). Applies to the current C++ implementation; noted here for the record. Python prototype guards against this.

---

## 3. Proposed algorithm: flipped iteration

Outer loop iterates right-anchors **right-to-left**, `x_0 = n - 1`. Uses only backward extends for MEM discovery.

### 3.1 Algorithm specification

At right-anchor `x`:

**Step 1' (backward extend, produces MEM).** From empty interval, backward-extend `P[x], P[x-1], ..., P[0]` until:
- Size drops below `min_occ` at some position `j`. Set `i = j + 1` (leftmost successful position). Interval `[sp, ep]` corresponds to substring `P[i..x]`.
- Or reach `j = -1` (extended all the way to the pattern's left boundary). Set `i = 0`.

If `x - i + 1 >= min_len`, emit MEM `[i, x + 1)` with the current interval `[sp, ep]` and size = `ep - sp + 1`.

**Step 2' (forward extend, identifies next-x).** From empty interval, forward-extend `P[i-1], P[i], P[i+1], ...` until:
- Size drops below `min_occ` at some position `j`. Return `j - 1` as the next anchor.
- Or reach `j = n` (extended all the way to the pattern's right boundary). Return `-1` (terminate outer loop).

Special case: if `i == 0` after Step 1', the MEM covers the pattern's left boundary. No `P[i-1]` exists. Return `-1` (terminate).

### 3.2 Correctness argument

**Property 1 (Step 1' produces left-maximal MEMs at anchor x).** The backward extend stops at the leftmost `i` such that `P[i..x]` has `>= min_occ` occurrences. By construction:
- If `i > 0`, then `P[i-1..x]` has `< min_occ` occurrences → left-maximal in the strict global sense.
- If `i = 0`, then `P[i..x]` reaches the pattern's left boundary → left-maximal trivially.

**Property 2 (Step 2' identifies the next MEM's right anchor).** After emitting MEM `P[i..x]`, we know:
- `P[i-1..x]` has `< min_occ` occurrences.
- Consequence: any MEM whose interval starts at position `<= i - 1` must end before `x`.

Step 2' forward-extends from `P[i-1]`. The extend succeeds through positions `i-1, i, i+1, ...` as long as size stays `>= min_occ`. Suppose it fails at position `j`. Then:
- `P[i-1..j-1]` has `>= min_occ` occurrences (the interval before adding `P[j]`).
- `P[i-1..j]` has `< min_occ` occurrences.
- Any MEM with left endpoint `<= i - 1` cannot extend right past `j - 1` (adding any character at position `>= j` would take a superstring of `P[i-1..j]` which has `< min_occ` occurrences).

Set next anchor `x' = j - 1`. At the next iteration, Step 1' will backward-extend from `x'` and find the true left endpoint of the MEM ending at `x'`. That left endpoint is guaranteed to be `<= i - 1` (because `P[i-1..j-1]` already succeeds), so the next MEM extends at least as far left as `i - 1`.

**Property 3 (no MEMs are missed).** Between anchors `x` and `x' = j - 1` (exclusive on the left, inclusive on the right), every possible right-anchor position corresponds to a right-endpoint that is either:
- Inside the just-emitted MEM `P[i..x]` — covered.
- Or equal to `j - 1` — will be processed at the next iteration.

Positions `x - 1, x - 2, ..., j` are strictly inside `P[i..x]` (they fall in the interior of the emitted MEM), so any MEM with those right-endpoints would be a substring of `P[i..x]` — not left-maximal, therefore not a MEM in its own right.

Positions `j - 1, j - 2, ..., i` are covered by the next MEM (Step 1' at `x' = j - 1` will discover a MEM whose interval spans them).

Positions `< i` are handled by subsequent iterations after that.

**Empirical validation.** `xy-test/mem-prototype/compare_mem_algorithms.py` implements both algorithms plus a brute-force MEM oracle (enumerate all `(a, b)` and check left+right maximality directly). Tested on:
- 14 handcrafted cases (boundaries, repeats, low-complexity, min_occ edge cases): 14/14 pass.
- 500 random cases, ACGT alphabet, text 20-100, pattern 3-15: 500/500 pass.
- 500 stress cases, mixed alphabets (AT / ACGT / ACGTN), text 30-300, pattern 4-30: 500/500 pass.

**Total: 1,014 / 1,014 test cases produce identical MEM sets across current, flipped, and brute-force algorithms.**

### 3.3 Comparison table with current algorithm

| Aspect | Current three-phase | Proposed flipped |
|---|---|---|
| Right-anchor iteration | Left-to-right | Right-to-left |
| Step 1 (seed) | Backward extend `min_len` chars | Backward extend until failure |
| Step 2 (MEM extension) | Forward extend to right-maximality | (folded into Step 1') |
| Step 3 (next-x) | Backward extend from MEM's right end | Forward extend from `P[i-1]` |
| MEM emission | After Step 2 | After Step 1' |
| Uses `forward_extend_encoded`? | Yes (Step 2, Step 3 uses backward) | Only Step 2' |
| SA-carry compatibility of MEM-discovery walk | Broken by Step 2 | Compatible (Step 1' is pure backward) |
| Extension work per MEM | ~ (x - i) + (e - x) + (e - x) chars | ~ (x - i) + (j - i + 1) chars |

Extension work is comparable in magnitude; the flipped algorithm's Step 2' work is roughly proportional to the discovered MEM length plus one extra character to the left.

---

## 4. SA-carrying mechanism during Step 1'

Adapt Section 3 of the r-index paper (Cobas, Gagie, Navarro 2021, arXiv:2103.15329) to our codebase's samples-at-run-*starts* convention.

### 4.1 Codebase invariants

From `include/pangenome_index/r-index.hpp:338-345`:

```cpp
// (sequence id, sequence offset) samples at the start of each run.
sdsl::int_vector<0> samples;

// Mark the text positions at the end of a run.
sdsl::sd_vector<> last;

// If last[i] = 1, last_to_run[last_rank(i)] is the identifier of the run.
sdsl::int_vector<0> last_to_run;
```

Key primitive semantics (verified in prior investigation):
- `samples[r]` = SA value at the start (leftmost BWT position) of run `r`.
- `locateNext(SA[k])` returns `SA[k+1]`, i.e., BWT-order-forward step. Used in `locate_sa_value` (`src/find_mems.cpp:178-191`) to walk from a run's start to any BWT position within it.
- No `locatePrev` primitive exists.

**Contrast with the r-index paper.** The paper's `Samples[p]` = `SA[j] - 1` where `j` is the *last* BWT position of the p-th run (run *end* samples, one-off from text position). Our `samples[r]` = SA value at the *first* BWT position of run `r` (run *start* samples). This is a mirror-image convention.

**Consequence:** the paper's mechanism gives `SA[ep]` (rightmost interval endpoint) during backward search. Adapting to our convention gives `SA[sp]` (leftmost interval endpoint) — which is exactly what we need for tag emission (since `locateNext` walks BWT-forward from `sp` through the interval).

### 4.2 SA carry update rule

**State carried through Step 1':** an SA value `sa_sp` that always equals `SA[sp]` for the current interval `[sp, ep]`.

**Initialization (before backward extend starts):** `sp` and `ep` are undefined (empty match). No SA value to carry yet. This is the "first character" case.

**First backward extend by `c`** (equivalent to the very first Step 1' iteration): interval becomes `[C[c] + 1, C[c] + count(c)]` covering all occurrences of character `c` in the BWT. `sp` = the leftmost such position. The run containing `sp` starts at `sp` itself (since `sp` is the first `c` in BWT order — it must be at a run boundary, either the start of a run of `c` or immediately after a run of a smaller character). Look up `samples[run_id_at(sp)]` = `SA[sp]`. Initialize `sa_sp = samples[run_id_at(sp)]`.

**Subsequent backward extends by `c`** (state: current interval `[sp, ep]`, current `sa_sp = SA[sp]`, extending by `c`):

Compute new interval `[sp', ep'] = backward_extend([sp, ep], c)` via existing `backward_extend_encoded` primitive. If empty, return failure (Step 1' terminates).

**Determine `SA[sp']`:**

Two cases based on whether the old `sp` position has BWT character equal to `c`:

**Case A (easy):** `BWT[sp] == c`.

The LF map sends `sp` to some position in the new interval. Specifically, if `sp` is the leftmost position in `[sp, ep]` with `BWT == c`, then `LF(sp) = sp'`. In this case, `SA[sp'] = SA[sp] - 1 = sa_sp - 1`.

But `sp` may not be the leftmost `c` in `[sp, ep]`; the leftmost `c` could be at a later position. We need to check whether `sp` itself has BWT == c.

Refined Case A: `BWT[sp] == c` AND `sp` is the leftmost c-position in `[sp, ep]` (i.e., `rank_c(BWT, sp - 1) < rank_c(BWT, sp)`, meaning `sp` is a c-position). Then `sp' = LF(sp)`, `sa_sp' = sa_sp - 1`.

**Case B (hard):** `BWT[sp] != c`, or `sp` has `BWT == c` but is not leftmost.

Need to find the leftmost c-position in `[sp, ep]`. Call it `j`. Then `sp' = LF(j)`, and `SA[j]` must be looked up.

Since `j` is the leftmost `c` in `[sp, ep]` and the c-positions form runs in the BWT, `j` is at a **run boundary** — specifically the start of a run of `c` (unless the entire interval is at run starts). By the samples-at-run-starts convention, `samples[run_id_at(j)]` gives `SA[j]` directly. Then `sa_sp' = SA[j] - 1`.

**Computing `j` (the leftmost c-position in `[sp, ep]`):**

Use `FastLocate::leftmost_c_in_interval(c, sp, ep)` (added in commit `1c9ff21`, see §4.3). Returns `SIZE_MAX` if no `c` exists in the interval — but in Case B we've already narrowed to a non-empty interval that must contain `c`, so a `SIZE_MAX` return here indicates a bug (assert in debug builds).

### 4.3 Required BWT primitives

Availability audited during Stage 0 (see §5.2). Status:

| Primitive | Available? | Where |
|---|---|---|
| `BWT[sp]` — character access at BWT position | ✅ | `FastLocate::bwt_char_at_encoded` (`src/r-index.cpp:519`; made public in commit `1eb0139`) |
| `rank_c(BWT, sp - 1)` — count of `c` in BWT prefix | ✅ | `FastLocate::rankAt_encoded` (`src/r-index.cpp:576`); also cached inside `backward_extend_encoded` at lines 757-758 |
| Leftmost `c` in `[sp, ep]` — replaces raw `select_c` | ✅ | `FastLocate::leftmost_c_in_interval` (added in commit `1c9ff21`) |
| `run_id_at(pos)` and run-start BWT position | ✅ | `FastLocate::run_id_and_offset_at` (`src/r-index.cpp:1212`) |
| `samples[run_id]` — SA at run start | ✅ | Existing `samples` field |

**On "leftmost c" rather than raw `select_c`:** the pangenome-index r-index does not store a general `select_c` over the BWT. Instead, commit `1c9ff21` added per-DNA-character run-start `sd_vector`s (one per char in ACGTN, `~2 MB total` per HPRCv1 chromosome), together with a query `leftmost_c_in_interval(c, sp, ep)` that returns the leftmost BWT position in `[sp, ep]` with `BWT[p] == c`. This is exactly the operation Case B needs, and it exploits the RLE structure (the answer is either `sp` itself when `BWT[sp] == c`, or the start of the leftmost c-run whose start falls in the interval). See `FastLocate::leftmost_c_in_interval` in `src/r-index.cpp:1752` and validation in `src/test_r_index_sa.cpp` (yeast: 105K+ tests, 0 failures).

**Overlap with `backward_extend_encoded`:** the existing `backward_extend_encoded` already computes `rank_at_cached_encoded(k)` and `rank_at_cached_encoded(k + s)` (lines 757-758). Case A's check "`BWT[sp] == c` and `sp` is leftmost c" reduces to `rank_cache_k[c] < rank_cache_ks[c]` AND `bwt_char_at_encoded(sp) == c`. The rank comparison is already computed; only the `bwt_char_at_encoded` call is added. So the marginal per-step cost of the SA-carry logic is one `bwt_char_at_encoded` call in the common case (Case A) plus one `leftmost_c_in_interval` call in the rare case (Case B).

### 4.4 Cost per backward extend step

- **Case A:** one rank query (to check `BWT[sp] == c` and confirm leftmost). Plus the rank/select ops already done by `backward_extend_encoded`. Marginal added cost: O(log(σ + n/r)) for one rank.

- **Case B:** one rank + one sample lookup. Same asymptotic cost.

- **Neither case walks `locateNext`.** This is the whole point — no per-step `locateNext` cost.

**Total Step 1' SA-carry cost:** O(|Step 1' walk|) rank queries. Compare to current `locate_sa_value(sp)` cost: O(`sp - run_start(sp)`) `locateNext` calls, which averages `avg_bwt_run_length / 2` steps.

For the HPRC 100k dataset: seed cost was 22.7M `locateNext` calls; Step 1' walk length averages `~x - i`, which is bounded by pattern length (say 150 bp for typical reads). So SA-carry adds up to `150 × 100,899 ≈ 15M` rank queries total. But most of these ranks are already done by `backward_extend_encoded` internally, so the *marginal* cost is much less.

**Net win estimate:** replace 22.7M `locateNext` calls (at ~0.24 µs each ≈ 5.4s) with ~15M marginal rank ops. Rank ops on `rank_at_cached` (`src/r-index.cpp:1515`) are cache-friendly and fast — likely << 0.24 µs each. Expected total gain in wall time: close to the full 5.4s of seed work, i.e., the projected 32% walk-cost reduction is largely realizable.

---

## 5. C++ implementation plan

### 5.0 Prerequisites completed

**Stage 0 — audit BWT primitives.** Confirmed all except one available (see §4.3 table).

**Stage 0.5-0.85 — resolve missing primitive.** Instead of a general `select_c` (would require a full BWT wavelet tree, ~25 MB per HPRCv1 chromosome), added per-DNA-character run-start `sd_vector`s (~2 MB) with query `leftmost_c_in_interval`. Rationale: exploits our existing RLE structure; the r-index already stores run boundaries via `blocks_start_pos`, so a mirror-image structure per character is a natural fit.

- Committed as `1c9ff21`: `r-index: add per-DNA-character run-start bitvectors for select_c` (primitive + `leftmost_c_in_interval` in `src/r-index.cpp:1752`, header field in `include/pangenome_index/r-index.hpp:356-365`).
- Committed as `1eb0139`: `test: add test_r_index_sa for leftmost_c_in_interval validation` (Tier-1 unit tests + selective exposure of `bwt_char_at_encoded` for tests to reach ground truth).

**Validation of Stage 0.85:**
- Unit: `bin/test_r_index_sa` on yeast chrII (`../mem-projection/pangenome-pipeline/runs/v2-yeast235/yeast235_chrII_100kb_normalized.ri`, n=337M, 14M runs) — 105K+ tests including structural point queries, edge cases, 100K random `(c, sp, ep)` triples. **Zero failures.**
- E2E: `./query.sh configs/hprcv1-chr6-alt-reads.env hprc-chr6-2026-06-02 local-validate --gaf` on HPRC chr6 100K alt reads — **100% valid (2000/2000), 0 invalid, 10.23M GAF entries.** find_mems wall 48s vs baseline 48s (0% delta), RSS 4.599 GB vs baseline 4.6 GB (0% delta). Expected zero regression since the new primitive has no callers yet on the hot path.

### 5.1 Files to modify (remaining stages)

- `include/pangenome_index/algorithm.hpp` — add new `find_all_mems_flipped` and `find_mems_flipped_function` alongside the existing ones. Do not replace; run side-by-side for validation. Per CLAUDE.md "Don't break the default path" — legacy behavior stays byte-identical.

- `include/pangenome_index/r-index.hpp` — add a new struct `bi_interval_with_sa` (or equivalent) carrying `bi_interval + sa_sp`, used only in the flipped code path. **Do not** extend the existing `bi_interval` — 5 existing call sites (r-index.cpp:791, algorithm.hpp:633/669/725, forward_extend_encoded) would each pay a memory-layout cost for a field they never use. Decision made during Stage 0 audit.

- `src/r-index.cpp` — add `backward_extend_encoded_with_sa` as a **new function** alongside the existing `backward_extend_encoded` (do not modify the existing; it's called from 5 sites, all of which want the original behavior). Signature draft:
  ```cpp
  struct bi_interval_with_sa { bi_interval bint; size_t sa_sp; };
  bi_interval_with_sa backward_extend_encoded_with_sa(
      const bi_interval& bint, size_t sa_sp, size_t c) const;
  ```
  Return an "empty" result (bint.size == 0) if the extend fails. The initial toehold (very first backward extend from empty match) must be handled specially — see §4.2.

- `src/find_mems.cpp` — add `dump_mem_info_lightweight_flipped` that receives the emitted MEM's `[sp, ep]` AND its `SA[sp]` (from the flipped MEM finder), skipping the `locate_sa_value(sp)` call entirely. Rest of the tag emission walk is identical.

### 5.2 Rollout stages

Verification bars per CLAUDE.md + VALIDATION_GUIDE.md — every stage that touches output-affecting code paths must pass `--gaf` validation on HPRC chr6 (100% valid / 0 invalid) plus GAF line-count parity with baseline (10,230,764 entries).

**Stage 1: SA-carry primitive standalone.**
- Implement `backward_extend_encoded_with_sa` in `src/r-index.cpp` alongside existing `backward_extend_encoded`.
- Extend `bin/test_r_index_sa` (or add a new binary `bin/test_backward_extend_sa`, following the same pattern) that:
  - Picks random substrings from the yeast index (do NOT use `xy-test/xy.ri` — known-bad fixture).
  - For each substring, runs backward-extend character-by-character both ways: (a) with SA carry via new function, (b) without, then calls existing `locate_sa_value(sp)` at each intermediate state.
  - Assert `sa_sp` from new path equals `locate_sa_value(sp)` at every step.
  - Also count Case A vs Case B frequencies. Log the ratio (informs Risk B in §5.3).
- Bar: 10K+ trials on yeast index, **zero mismatches**.
- E2E: not applicable (new function has no callers yet). Skip.

**Stage 2: Flipped MEM finder.**
- Implement `find_all_mems_flipped` in `include/pangenome_index/algorithm.hpp` using `backward_extend_encoded_with_sa` for Step 1' and existing `forward_extend_encoded` for Step 2'. Step 2' doesn't need SA carry (its purpose is only to identify next-x, not to emit MEMs).
- Write a test that runs both `find_all_mems` and `find_all_mems_flipped` on the canonical yeast reads (`../mem-projection/yeast-235/yeast-235-chrI/S288C_chrII_N100K_R1_200_reads.txt`, subsetted via `max_reads` for iteration speed), diffs MEM sets.
- Bar: bitwise-identical `(start, end, bwt_start, size)` tuples across every yeast read. Also verify `sa_sp` returned by flipped equals `locate_sa_value(bwt_start)` computed via the legacy path — this is the correctness check that ties the algorithm-level result back to the primitive-level guarantee.
- E2E: not applicable yet (find_mems still uses the legacy path).

**Stage 3: Zero-seed-cost dump path.**
- Implement `dump_mem_info_lightweight_flipped` in `src/find_mems.cpp` that takes `(mem, sa_sp)` and emits tags via `locateNext` from `sa_sp`. **No `locate_sa_value` call.**
- Wire into `find_mems` main behind a new flag `--use-flipped-mems` (default off). Default path remains byte-identical per CLAUDE.md.
- Bar (unit): TSV output with `--use-flipped-mems` on yeast 500 reads is bitwise-identical to legacy output (after sorting, since order may differ).
- Bar (E2E): `./query.sh configs/hprcv1-chr6-alt-reads.env hprc-chr6-2026-06-02 stage3-flipped --gaf` with the flipped path enabled. **Must show 100% valid, 0 invalid.** GAF line count must match baseline (10,230,764 entries) exactly — line-count parity catches silently-dropped MEMs that `validate_gaf` would miss.

**Stage 4: Benchmark.**
- Run flipped vs legacy on HPRC chr6 100K reads via `./query.sh`. Compare `timing_summary.txt`.
- Target: find_mems wall reduction ≈ 10% (per §4.4 projection).
- Bar: no RSS regression >10%; --gaf still 100% valid as regression check.
- If wall doesn't drop, or drops less than 3%: revisit Case B frequency (Risk B in §5.3). If Case B dominates due to short intervals inside long runs, mitigation options open.

**Stage 5: Ship.**
- If Stage 4 confirms projection: make flipped path the default. Keep legacy behind `--use-legacy-mems` for regression testing.
- Update `FINDINGS_SEED_COST.md` with measured results (replace "projected" with "measured" numbers).
- Local commits ready for eventual push (user has instructed no automatic pushes).

### 5.3 Risks

**Risk A: BWT primitives not available.** ~~If `select_c` or leftmost-c-in-interval turns out to require a new BWT scan or a significantly more expensive derivation, Case B of the SA carry rule becomes expensive.~~ **RESOLVED in Stage 0.85.** Added `leftmost_c_in_interval` via per-DNA-character run-start `sd_vector`s. Storage overhead ~2 MB per HPRCv1 chromosome (measured by E2E: RSS unchanged from baseline). Case B query is O(1) fast-path (BWT[sp] == c) + O(log r) fallback (successor on the per-c bitvector).

**Risk B: Case B is more common than expected.** The projected win assumes Case A dominates (fast LF-decrement path). If Case B fires on most backward extends, per-step cost climbs and win diminishes. Mitigation: **Stage 1's test must include Case A / Case B frequency counters** and log the ratio. If Case B > 50% of steps on realistic reads, revisit design (options: (a) larger per-c bitvector coverage, (b) fall back to `locate_sa_value` when Case B is deep in a long run so amortized savings still win).

**Risk C: Step 2' doubles pattern-extension work.** The current algorithm does Step 1 (min_len chars) + Step 2 (MEM extension) + Step 3 (up to MEM length). The flipped does Step 1' (MEM length) + Step 2' (MEM length + 1). Total extension work roughly comparable but the split is different. Actual wall impact depends on `backward_extend_encoded` vs `forward_extend_encoded` cost ratio (forward_extend is swap → backward_extend → swap, so ~2x per char). Mitigation: measure in Stage 4; if Step 2' overhead cancels seed cost savings, evaluate whether Step 2' can be trimmed (e.g., early termination heuristics).

**Risk D: bugs during backward extend's early termination.** If Step 1' terminates because interval size drops below `min_occ` at position `j`, the SA carry state was updated for the *pre-termination* interval, not the failed one. Need to be careful about which state is returned and how it interacts with the "no MEM emitted" outcome. Mitigation: the Python prototype at `xy-test/mem-prototype/compare_mem_algorithms.py` handles this correctly (state is the last-successful interval's SA, not the failed one). C++ implementation should mirror the Python state machine exactly. Stage 1's unit test explicitly covers this by verifying `sa_sp` matches `locate_sa_value` at every step, including the last successful one before failure.

**Risk E: the current algorithm's Step 3 out-of-bounds bug** (accessing `P[len]` when Step 2 extends to end). Doesn't affect the flipped path directly (Step 2' handles the boundary condition), but should be fixed independently to avoid latent UB in the current path. Mitigation: fix as a separate commit before adding new code paths. Tracked as an independent todo item.

**Risk F: Validation surface is downstream, not local.** find_mems output feeds gafpack → validate_gaf → cosigt. A bug in the flipped path could produce MEMs that pass local unit tests but corrupt downstream variant calls. Mitigation: Stage 3 gate requires `./query.sh ... --gaf` to show 100% valid AND GAF line-count parity. Anything less is a hard stop. See VALIDATION_GUIDE.md.

---

## 6. Objective, restated

**Delivering `SA[sp]` for each MEM's interval as a byproduct of the MEM-discovery walk, at the cost of restructuring MEM enumeration to iterate right-to-left with pure-backward-extend Step 1'.** Removes 32% of the walk cost in `dump_mem_info_lightweight`, projected 10% end-to-end `find_mems` wall reduction on HPRCv1 chr6.

Correctness of the flipped enumeration is empirically verified on 1,014 diverse test cases via the Python prototype at `xy-test/mem-prototype/compare_mem_algorithms.py`.

SA-carrying mechanism is a direct adaptation of Cobas/Gagie/Navarro 2021 §3, mirrored to our samples-at-run-*starts* convention (paper uses run-*ends*). Required BWT primitives are all either already available or standard.

---

## 7. Appendix: prototype-verified behavior on edge cases

From `xy-test/mem-prototype/compare_mem_algorithms.py`:

| Case | Behavior |
|---|---|
| Pattern extends to text end (`ACGT` in `ACGTACGT`) | Single MEM `[0,4)`, both algorithms match ground truth. |
| All-same pattern (`AAAA` in `AAAA...`) | Single MEM covering the pattern, both match. |
| Overlapping MEMs (`ACGACGT` in `ACGACG...ACGT...`) | Two MEMs `[0,6)` and `[3,7)`, both match. |
| Low-complexity (`AAAAAA` with mixed text) | Multiple overlapping MEMs enumerated correctly. |
| Pattern doesn't match at all | Empty MEM list, both match. |
| min_occ higher than any occurrence | Empty MEM list, both match. |
| MEM at exact pattern boundary (start or end) | Correctly emitted with proper boundary. |
| Random small (500 tests, ACGT, text 20-100, pattern 3-15) | 500/500 identical. |
| Random stress (500 tests, mixed alphabets AT/ACGT/ACGTN, text 30-300, pattern 4-30) | 500/500 identical. |

No divergences observed. High confidence in algorithmic correctness.
