#!/usr/bin/env python3
"""
Prototype comparing the current three-phase MEM algorithm against the
proposed flipped-iteration algorithm.

Context:
  - The seed_cost in dump_mem_info_lightweight is 32% of walk on HPRCv1 chr6.
  - Proposed optimization: iterate right-anchors right-to-left, backward-extend
    from x leftward (which admits r-index-style SA carrying because it's pure
    backward extend), then forward-extend from i+1 up to x to identify the next
    right-anchor. This eliminates the seed cost entirely.
  - Question this prototype answers: does the flipped algorithm enumerate the
    same MEMs as the current three-phase algorithm?

We do NOT model SA carrying here -- that's a separate concern once algorithm
correctness is settled. We use a naive occurrence oracle (Python substring
search) to abstract away the r-index primitives.

Reference for algorithm semantics:
  include/pangenome_index/algorithm.hpp:655-738 (find_mems_function)
  include/pangenome_index/algorithm.hpp:741-751 (find_all_mems)
"""

import sys
import random
from dataclasses import dataclass
from typing import List, Set, Tuple


@dataclass(frozen=True)
class MEM:
    """A MEM emission. Matches C++ MEM struct fields we care about.

    NOTE: current C++ emits (start=x, end=e, bwt_start, size). end is
    exclusive; the actual MEM string is P[start:end]. For comparison we
    use the string content plus size.
    """
    start: int
    end: int  # exclusive
    size: int  # number of occurrences

    def string(self, pattern: str) -> str:
        return pattern[self.start:self.end]


class OccOracle:
    """Naive occurrence counter: number of times P[a:b] appears in text.

    Substitute for the r-index count() primitive. O(|text|) per query,
    but that's fine for small prototype tests.
    """
    def __init__(self, text: str):
        self.text = text

    def count(self, s: str) -> int:
        if not s:
            return len(self.text) + 1  # empty string matches "everywhere"
        # naive substring count with overlap
        c = 0
        i = 0
        while True:
            j = self.text.find(s, i)
            if j < 0:
                break
            c += 1
            i = j + 1
        return c


# ---------------------------------------------------------------------------
# Current three-phase algorithm, faithful transcription of algorithm.hpp:655-738
# ---------------------------------------------------------------------------

def find_mems_current_at_x(pattern: str, min_len: int, min_occ: int, x: int,
                            oracle: OccOracle, mems_out: List[MEM]) -> int:
    """Faithful port of find_mems_function.

    Returns the next x for the outer loop.
    Emits at most one MEM into mems_out.
    """
    n = len(pattern)
    if n - x < min_len:
        return n

    # Step 1: backward-extend P[x + min_len - 1] down to P[x].
    # We track the interval via its string content because the oracle is naive.
    # The C++ code initializes bint = full-BWT and backward-extends. Here we
    # track the substring being matched.
    start = x + min_len - 1
    j = start
    current_str = ""
    while j >= x:
        current_str = pattern[j] + current_str
        size = oracle.count(current_str)
        if size < min_occ or size <= 0:
            return j + 1  # backward-extend failed at j
        if j == 0:
            break
        j -= 1

    # After Step 1, current_str = P[x:x+min_len], size >= min_occ.
    # Step 2: forward-extend P[x + min_len], P[x + min_len + 1], ...
    e = x + min_len
    last_good_str = current_str
    last_good_size = oracle.count(current_str)
    while e < n:
        trial = current_str + pattern[e]
        size = oracle.count(trial)
        if size < min_occ or size <= 0:
            break
        current_str = trial
        last_good_str = current_str
        last_good_size = size
        e += 1
    # After Step 2, MEM emitted is P[x:e] with occurrence count last_good_size.
    mems_out.append(MEM(start=x, end=e, size=last_good_size))

    # Step 3: fresh backward-extend P[e], P[e-1], ..., P[x+1] from full-BWT.
    # NOTE: this accesses P[e] which is out of bounds if e == n. The C++ code
    # has the same issue -- documented as a potential bug in PROOF_FLIPPED_MEM.md.
    # For the prototype, we handle it by returning n (terminate outer loop) in
    # that case, which is the "safe" behavior most likely intended.
    if e >= n:
        return n

    back_str = ""
    j = e
    while j > x:
        back_str = pattern[j] + back_str
        size = oracle.count(back_str)
        if size < min_occ or size <= 0:
            return j + 1
        j -= 1
    return j + 1


def find_all_mems_current(pattern: str, min_len: int, min_occ: int,
                           oracle: OccOracle) -> List[MEM]:
    mems: List[MEM] = []
    x = 0
    n = len(pattern)
    while x < n:
        x = find_mems_current_at_x(pattern, min_len, min_occ, x, oracle, mems)
    return mems


# ---------------------------------------------------------------------------
# Flipped algorithm (proposed).
#
# Iterate right-anchors right-to-left. For each x:
#   Step 1': backward-extend P[x], P[x-1], ... until size < min_occ or reach P[0].
#            Let i = leftmost successful position. MEM interval = P[i..x].
#   Step 2': fresh forward-extend P[i+1], P[i+2], ..., P[j] until size < min_occ
#            or j reaches x.
#     If fail at some j' <= x: return j' - 1 as next x.
#     If reached x without failing: return i - 1 as next x.
# ---------------------------------------------------------------------------

def find_mems_flipped_at_x(pattern: str, min_len: int, min_occ: int, x: int,
                            oracle: OccOracle, mems_out: List[MEM]) -> int:
    """Flipped algorithm at right-anchor x. Returns next x (may be -1).

    User-corrected Step 2':
      After Step 1' gives us MEM P[i..x], we know P[i-1..x] does NOT match.
      Start fresh, forward-extend P[i-1], P[i], P[i+1], ... until it fails at
      some position j. Return j - 1 as next x.
      (If i == 0, no P[i-1] to start from; return -1.)
    """
    n = len(pattern)
    if x + 1 < min_len:
        return -1

    # Step 1': backward-extend from x leftward.
    current_str = ""
    i = x + 1  # will hold leftmost successful position
    j = x
    while j >= 0:
        trial = pattern[j] + current_str
        size = oracle.count(trial)
        if size < min_occ or size <= 0:
            break
        current_str = trial
        i = j
        j -= 1
    # After loop, current_str = P[i..x], size >= min_occ.
    # If we couldn't extend at all (i > x), no MEM possible at this anchor.
    if i > x:
        return x - 1
    # MEM length check
    mem_len = x - i + 1
    if mem_len < min_len:
        # Length insufficient; still need to advance. Use same next-x mechanism.
        pass  # fall through to Step 2' logic, but don't emit
    else:
        mem_size = oracle.count(current_str)
        mems_out.append(MEM(start=i, end=x + 1, size=mem_size))

    # Step 2': fresh forward-extend starting at position i-1.
    # If i == 0, there's no character before the MEM to seed the forward extend.
    # In that case, the MEM already covers the pattern's left boundary, and
    # we've enumerated all possible MEMs -- terminate.
    if i == 0:
        return -1

    fwd_str = ""
    j = i - 1
    failed_at = -1
    while j < n:
        trial = fwd_str + pattern[j]
        size = oracle.count(trial)
        if size < min_occ or size <= 0:
            failed_at = j
            break
        fwd_str = trial
        j += 1

    if failed_at >= 0:
        # Failed at position failed_at; the successful extend ended at failed_at - 1.
        # That's the new x.
        return failed_at - 1
    else:
        # Reached end of pattern without failing. All positions from i-1 onward
        # extend successfully. Terminate.
        return -1


def find_all_mems_flipped(pattern: str, min_len: int, min_occ: int,
                           oracle: OccOracle) -> List[MEM]:
    mems: List[MEM] = []
    x = len(pattern) - 1
    while x >= 0:
        x = find_mems_flipped_at_x(pattern, min_len, min_occ, x, oracle, mems)
    return mems


# ---------------------------------------------------------------------------
# Ground truth: enumerate all MEMs by brute force.
#
# A MEM is a substring P[a..b] (inclusive) such that:
#   1. length b - a + 1 >= min_len
#   2. count(P[a..b]) >= min_occ
#   3. left-maximal: a == 0 or count(P[a-1..b]) < min_occ
#   4. right-maximal: b == n-1 or count(P[a..b+1]) < min_occ
# ---------------------------------------------------------------------------

def enumerate_ground_truth_mems(pattern: str, min_len: int, min_occ: int,
                                 oracle: OccOracle) -> Set[Tuple[int, int, int]]:
    """Return set of (start, end_exclusive, size) tuples for all MEMs."""
    n = len(pattern)
    mems = set()
    for a in range(n):
        for b in range(a + min_len - 1, n):
            s = pattern[a:b + 1]
            cnt = oracle.count(s)
            if cnt < min_occ:
                break  # extending further right only makes count smaller
            # left-maximal check
            if a > 0:
                left_ext = pattern[a - 1:b + 1]
                if oracle.count(left_ext) >= min_occ:
                    continue
            # right-maximal check
            if b < n - 1:
                right_ext = pattern[a:b + 2]
                if oracle.count(right_ext) >= min_occ:
                    continue
            mems.add((a, b + 1, cnt))
    return mems


# ---------------------------------------------------------------------------
# Comparison utilities
# ---------------------------------------------------------------------------

def mems_to_set(mems: List[MEM]) -> Set[Tuple[int, int, int]]:
    return {(m.start, m.end, m.size) for m in mems}


def compare(pattern: str, min_len: int, min_occ: int, text: str,
            verbose: bool = False) -> Tuple[bool, bool, bool]:
    """Return (current_ok, flipped_ok, current_equals_flipped)."""
    oracle = OccOracle(text)
    current = mems_to_set(find_all_mems_current(pattern, min_len, min_occ, oracle))
    flipped = mems_to_set(find_all_mems_flipped(pattern, min_len, min_occ, oracle))
    ground = enumerate_ground_truth_mems(pattern, min_len, min_occ, oracle)

    current_ok = current == ground
    flipped_ok = flipped == ground
    same = current == flipped

    if verbose or not (current_ok and flipped_ok and same):
        print(f"pattern={pattern!r} min_len={min_len} min_occ={min_occ}")
        print(f"  ground truth: {sorted(ground)}")
        print(f"  current:      {sorted(current)}   ok={current_ok}")
        print(f"  flipped:      {sorted(flipped)}   ok={flipped_ok}")
        print(f"  same:         {same}")
        missing_from_current = ground - current
        extra_in_current = current - ground
        missing_from_flipped = ground - flipped
        extra_in_flipped = flipped - ground
        if missing_from_current:
            print(f"  current missing: {sorted(missing_from_current)}")
        if extra_in_current:
            print(f"  current extra:   {sorted(extra_in_current)}")
        if missing_from_flipped:
            print(f"  flipped missing: {sorted(missing_from_flipped)}")
        if extra_in_flipped:
            print(f"  flipped extra:   {sorted(extra_in_flipped)}")

    return current_ok, flipped_ok, same


# ---------------------------------------------------------------------------
# Test cases
# ---------------------------------------------------------------------------

def handcrafted_tests():
    cases = [
        # (name, pattern, text, min_len, min_occ)
        ("simple",        "ACGT",     "ACGTACGT",         2, 2),
        ("all_same",      "AAAA",     "AAAAAAAAAA",       2, 3),
        ("no_match",      "XYZW",     "ACGTACGT",         2, 1),
        ("full_pattern",  "ACGT",     "ACGTACGT" * 3,     2, 2),
        ("boundary_L",    "ACGT",     "ACG" + "T" * 20,   2, 3),
        ("boundary_R",    "ACGT",     "A" * 20 + "CGT",   2, 3),
        ("nested_reps",   "ACGACG",   "ACGACGACG" * 2,    3, 2),
        ("low_complex",   "AAAAAA",   "AAAA" + "T" * 5 + "AAAAAA", 2, 2),
        ("min_len_edge",  "AC",       "ACACAC",           2, 2),
        ("min_occ_high",  "ACGT",     "ACGT",             2, 5),  # unreachable
        ("mem_at_start",  "TACGTAA",  "TACGT" + "X" * 3 + "AACC", 3, 1),
        ("multi_mems",    "AACGT",    "AAAA" + "CGT" * 4 + "AACGT" * 2, 2, 2),
        # Explicit case exercising the P[len] step-3 access
        ("extend_to_end", "ACGT",     "ACGT" * 5,         2, 3),
        # Case where flipped might diverge: multiple MEMs with overlapping structure
        ("interleaved",   "ACGACGT",  "ACGACG" * 3 + "ACGT" * 3, 3, 2),
    ]

    print("=== Handcrafted tests ===")
    n_pass = 0
    n_fail = 0
    n_current_bad = 0
    n_flipped_bad = 0
    n_diff = 0
    for name, pattern, text, min_len, min_occ in cases:
        print(f"\n--- {name} ---")
        cur_ok, flp_ok, same = compare(pattern, min_len, min_occ, text, verbose=True)
        if cur_ok and flp_ok and same:
            n_pass += 1
        else:
            n_fail += 1
            if not cur_ok:
                n_current_bad += 1
            if not flp_ok:
                n_flipped_bad += 1
            if not same:
                n_diff += 1
    print(f"\n=== Handcrafted summary: {n_pass}/{n_pass+n_fail} pass ===")
    print(f"  current wrong: {n_current_bad}, flipped wrong: {n_flipped_bad}, diverged: {n_diff}")
    return n_fail


def random_tests(n_tests: int, seed: int = 42):
    print(f"\n=== Random tests (n={n_tests}, seed={seed}) ===")
    random.seed(seed)
    alphabet = "ACGT"
    n_pass = 0
    n_fail = 0
    n_current_bad = 0
    n_flipped_bad = 0
    n_diff = 0
    failures = []
    for i in range(n_tests):
        text_len = random.randint(20, 100)
        pattern_len = random.randint(3, 15)
        min_len = random.randint(2, min(4, pattern_len))
        min_occ = random.randint(1, 4)
        text = "".join(random.choice(alphabet) for _ in range(text_len))
        pattern = "".join(random.choice(alphabet) for _ in range(pattern_len))
        cur_ok, flp_ok, same = compare(pattern, min_len, min_occ, text, verbose=False)
        if cur_ok and flp_ok and same:
            n_pass += 1
        else:
            n_fail += 1
            if not cur_ok:
                n_current_bad += 1
            if not flp_ok:
                n_flipped_bad += 1
            if not same:
                n_diff += 1
            failures.append((pattern, text, min_len, min_occ, cur_ok, flp_ok, same))

    print(f"=== Random summary: {n_pass}/{n_tests} pass ===")
    print(f"  current wrong: {n_current_bad}, flipped wrong: {n_flipped_bad}, diverged: {n_diff}")

    # Print first few failures in detail
    for pattern, text, min_len, min_occ, cur_ok, flp_ok, same in failures[:5]:
        print(f"\n--- FAILURE (first 5) ---")
        compare(pattern, min_len, min_occ, text, verbose=True)

    return n_fail


def stress_tests(n_tests: int, seed: int = 137):
    """Larger, more adversarial random tests."""
    print(f"\n=== Stress tests (n={n_tests}, seed={seed}) ===")
    random.seed(seed)
    n_pass = 0
    n_fail = 0
    n_current_bad = 0
    n_flipped_bad = 0
    n_diff = 0
    failures = []
    for i in range(n_tests):
        # Mix of alphabets: small (2 chars, high repetition), medium (DNA), large
        alph_choice = random.random()
        if alph_choice < 0.25:
            alphabet = "AT"  # binary-like, forces many repeats
        elif alph_choice < 0.75:
            alphabet = "ACGT"
        else:
            alphabet = "ACGTN"
        # Wider size range
        text_len = random.randint(30, 300)
        pattern_len = random.randint(4, 30)
        min_len = random.randint(2, min(6, pattern_len))
        min_occ = random.randint(1, 6)
        text = "".join(random.choice(alphabet) for _ in range(text_len))
        pattern = "".join(random.choice(alphabet) for _ in range(pattern_len))
        cur_ok, flp_ok, same = compare(pattern, min_len, min_occ, text, verbose=False)
        if cur_ok and flp_ok and same:
            n_pass += 1
        else:
            n_fail += 1
            if not cur_ok:
                n_current_bad += 1
            if not flp_ok:
                n_flipped_bad += 1
            if not same:
                n_diff += 1
            failures.append((pattern, text, min_len, min_occ, cur_ok, flp_ok, same))

    print(f"=== Stress summary: {n_pass}/{n_tests} pass ===")
    print(f"  current wrong: {n_current_bad}, flipped wrong: {n_flipped_bad}, diverged: {n_diff}")

    for pattern, text, min_len, min_occ, cur_ok, flp_ok, same in failures[:5]:
        print(f"\n--- STRESS FAILURE ---")
        compare(pattern, min_len, min_occ, text, verbose=True)

    return n_fail


if __name__ == "__main__":
    fail_hand = handcrafted_tests()
    fail_rand = random_tests(500)
    fail_stress = stress_tests(500)
    total = fail_hand + fail_rand + fail_stress
    print(f"\n=== TOTAL: handcrafted={fail_hand}, random={fail_rand}, stress={fail_stress} ===")
    sys.exit(1 if total else 0)
