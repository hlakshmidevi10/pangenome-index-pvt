# Design: Zero-Seed-Cost MEM Finder via Flipped Iteration + R-Index SA Carrying

**Status:** Design document. Algorithm correctness verified empirically (1,014 tests). SA-carrying mechanism specified but not yet implemented. Ready for C++ implementation.

**Branch:** `classify-mem-runs`

**Companion artifacts:**
- `FINDINGS_SEED_COST.md` — measurement data motivating this design
- `xy-test/mem-prototype/compare_mem_algorithms.py` — working Python prototype

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

Standard formulation: `j = C[c] + rank_c(BWT, sp - 1) + 1` gives the position of the `(rank_c(BWT, sp - 1) + 1)`-th c-occurrence in the BWT. Actually that's the LF formula for the smallest c-occurrence to the right of `sp - 1`. Equivalently: `j = select_c(rank_c(BWT, sp - 1) + 1)`. Requires a select-by-character primitive on the BWT.

### 4.3 Required BWT primitives

The mechanism needs (per backward-extend step):

1. **`BWT[sp]`** — access character at position `sp`. Should be cheap (rank/select on the BWT bitvector or equivalent). Available via `rankAt_encoded` if we're careful about API.

2. **`rank_c(BWT, sp - 1)`** — number of c-occurrences in `BWT[0..sp-1]`. This is the standard rank primitive; used by `backward_extend_encoded` internally.

3. **`select_c(BWT, k)`** — position of the k-th c-occurrence in the BWT. Needed only in Case B. May not be directly available; if not, `j = C[c] + rank_c(BWT, sp - 1) + 1` gives it via the same rank + `C[]` array we already have.

4. **`run_id_at(pos)`** — run id containing BWT position `pos`. Available (used by `run_id_and_offset_at` at `src/r-index.cpp:1212`).

5. **`samples[run_id]`** — SA value at the run's start. Available (already used).

**Primitive availability audit before implementation:** need to verify (1), (2), (3) are directly available or trivially derivable from existing calls. `backward_extend_encoded` already computes some of these internally; may be able to extend it to expose the values we need instead of computing them twice.

### 4.4 Cost per backward extend step

- **Case A:** one rank query (to check `BWT[sp] == c` and confirm leftmost). Plus the rank/select ops already done by `backward_extend_encoded`. Marginal added cost: O(log(σ + n/r)) for one rank.

- **Case B:** one rank + one sample lookup. Same asymptotic cost.

- **Neither case walks `locateNext`.** This is the whole point — no per-step `locateNext` cost.

**Total Step 1' SA-carry cost:** O(|Step 1' walk|) rank queries. Compare to current `locate_sa_value(sp)` cost: O(`sp - run_start(sp)`) `locateNext` calls, which averages `avg_bwt_run_length / 2` steps.

For the HPRC 100k dataset: seed cost was 22.7M `locateNext` calls; Step 1' walk length averages `~x - i`, which is bounded by pattern length (say 150 bp for typical reads). So SA-carry adds up to `150 × 100,899 ≈ 15M` rank queries total. But most of these ranks are already done by `backward_extend_encoded` internally, so the *marginal* cost is much less.

**Net win estimate:** replace 22.7M `locateNext` calls (at ~0.24 µs each ≈ 5.4s) with ~15M marginal rank ops. Rank ops on `rank_at_cached` (`src/r-index.cpp:1515`) are cache-friendly and fast — likely << 0.24 µs each. Expected total gain in wall time: close to the full 5.4s of seed work, i.e., the projected 32% walk-cost reduction is largely realizable.

---

## 5. C++ implementation plan

### 5.1 Files to modify

- `include/pangenome_index/algorithm.hpp` — add new `find_all_mems_flipped` and `find_mems_flipped_function` alongside the existing ones (do not replace; run side-by-side for validation).

- `include/pangenome_index/r-index.hpp` — extend `bi_interval` to carry an `sa_sp` field (an optional SA value tracked through backward extends). Or add a new struct `bi_interval_with_sa` used only in the flipped code path.

- `src/r-index.cpp` — add a variant `backward_extend_encoded_with_sa` that computes the standard interval extension AND the SA carry rule. Do not modify the existing `backward_extend_encoded` (used by other code paths).

- `src/find_mems.cpp` — add a lightweight dump path `dump_mem_info_lightweight_flipped` that receives the emitted MEM's `[sp, ep]` AND its `SA[sp]` (from the flipped MEM finder), skipping the `locate_sa_value(sp)` call entirely. Rest of the tag emission walk is identical.

### 5.2 Rollout stages

**Stage 1: SA-carry primitive standalone.**
- Implement `backward_extend_encoded_with_sa` in isolation.
- Write a C++ test (`xy-test/sa_carry_test.cpp` or `tests/`) that: pick random substrings of the test index, run backward extend both ways (with and without SA carry), verify `sa_sp` equals `locate_sa_value(sp)` at each step.
- Run on `xy.ri` and yeast index. Expect zero mismatches across ~10K trials.

**Stage 2: Flipped MEM finder.**
- Implement `find_all_mems_flipped` using `backward_extend_encoded_with_sa` for Step 1' and existing `forward_extend_encoded` for Step 2' (Step 2' doesn't need SA carry).
- Write a test that runs both `find_all_mems` and `find_all_mems_flipped` on the same reads, diffs the resulting MEM sets. Expect identical MEM sets across all yeast reads.

**Stage 3: Zero-seed-cost dump path.**
- Implement `dump_mem_info_lightweight_flipped` that takes `(mem, sa_sp)` and emits tags via `locateNext` from `sa_sp`. No `locate_sa_value` call.
- Wire it into `find_mems` main as an alternative code path behind a flag (e.g., `--use-flipped-mems`).
- Verify TSV output is identical to the current path (bitwise, sorted).

**Stage 4: Benchmark.**
- Run both paths on yeast 500 reads and HPRCv1 chr6 (both 1k sample and full 100k).
- Compare wall time via `--time-per-read` or the existing `SUMMARY.tsv` timing.
- Expected: `find_mems` phase drops by ~10% end-to-end on HPRCv1 chr6 full.

**Stage 5: Ship.**
- If Stage 4 confirms projection, make flipped path the default. Keep current path behind `--use-legacy-mems` for regression testing.
- Update `FINDINGS_SEED_COST.md` with measured results.
- Commit + push.

### 5.3 Risks

**Risk A: BWT primitives not available.** If `select_c` or leftmost-c-in-interval turns out to require a new BWT scan or a significantly more expensive derivation, Case B of the SA carry rule becomes expensive. Mitigation: audit primitive availability at Stage 0 (before any code); if missing, consider precomputing what's needed or restructuring the invariant to avoid Case B (e.g., always track SA at `ep` instead of `sp`, then walk `locateNext` from a stored `SA[sp]` at Step 1' completion — reverts partial benefit).

**Risk B: Case B is more common than expected.** The projected win assumes Case A dominates (fast LF-decrement path). If Case B fires on most backward extends, per-step cost climbs and win diminishes. Mitigation: instrument the primitive with Case A / Case B counters during Stage 1 testing. If Case B > 50% of steps, revisit design.

**Risk C: Step 2' doubles pattern-extension work.** The current algorithm does Step 1 (min_len chars) + Step 2 (MEM extension) + Step 3 (up to MEM length). The flipped does Step 1' (MEM length) + Step 2' (MEM length + 1). Total extension work roughly comparable but the split is different. Actual wall impact depends on `backward_extend_encoded` vs `forward_extend_encoded` cost ratio. Mitigation: measure in Stage 4; if Step 2' overhead cancels seed cost savings, reconsider.

**Risk D: bugs during backward extend's early termination.** If Step 1' terminates because interval size drops below `min_occ` at position `j`, the SA carry state was updated for the *pre-termination* interval, not the failed one. Need to be careful about which state is returned and how it interacts with the "no MEM emitted" outcome. Mitigation: prototype the exact state machine in Python before C++.

**Risk E: the current algorithm's Step 3 out-of-bounds bug** (accessing `P[len]` when Step 2 extends to end). Doesn't affect the flipped path directly (Step 2' handles the boundary condition), but should be fixed independently to avoid latent UB in the current path. Mitigation: fix as a separate commit before adding new code paths.

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
