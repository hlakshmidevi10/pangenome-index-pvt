# Project TODO

Lightweight backlog for cross-session follow-ups that aren't yet a JR entry
or in-flight work. Keep entries short; promote to a RESEARCH_JOURNAL.md
entry once actually worked.

- [ ] **Fix CLAUDE.md doc gap.** `CLAUDE.md`'s "Where to read first" section
  (bottom of file) points to `VALIDATION_GUIDE.md`, `ALGORITHM_CHANGES_SUMMARY.md`,
  and `mem_projection_1.pdf` as required reading. None of these three exist
  in the repo or anywhere in git history (confirmed via `git log --all`) —
  they were referenced but never created/committed. Either write them or
  remove/update the references.
- [ ] **Write `VALIDATION_GUIDE.md`.** Referenced as "the full validation
  procedure" from the top of `CLAUDE.md` and cited repeatedly through
  RESEARCH_JOURNAL.md (e.g. JR-007's "coverage-file MD5 as the preferred
  check" recommendation), but never actually written. Should consolidate:
  the `--gaf` / `validate_gaf` 100%-valid gate, GAF line-count parity check,
  coverage-file MD5 comparison (JR-007's stronger signal), and the
  `bwt_intervals size == .seq size + 1` invariant check — currently these
  live scattered across `CLAUDE.md`'s "Verification before claiming done"
  table and `mem-projection/pangenome-pipeline/CLAUDE.md`'s "Correctness
  criterion" / "Quick triage" sections.
