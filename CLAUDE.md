# pangenome-index-latest — agent guide

> **Validation is mandatory.** Before claiming any change is complete, run the validation pipeline. See `VALIDATION_GUIDE.md` for the full procedure — the short version is:
> ```bash
> make bin/find_mems  # rebuild
> cd ../mem-projection/pangenome-pipeline
> ./query.sh configs/hprcv1-chr6-alt-reads.env hprc-chr6-2026-06-02 <test-name> --gaf
> # Must show: Valid entries: 2000 (100.00%)
> ```

> **Read the research journal first.** `RESEARCH_JOURNAL.md` is the append-only
> record of hypotheses, experiments, measurements, and open questions on this
> codebase. Skim the index at the top before starting new work — you may find
> your question already answered, or the reasoning behind a design decision
> already documented. When you do non-trivial work, add an entry. Rules for
> new entries are at the top of that file.

## Engineering principles (read first)

This is **high-throughput, low-latency software**. It runs the inner loop of MEM search and tag-array projection on whole-genome pangenomes. Indexes are GB-scale, query workloads are millions of MEMs, and the entire toolkit feeds downstream genotyping pipelines (cosigt) where every ms of `find_mems` wall and every MB of peak RSS shows up in production cohort runs. **A single invalid MEM projection corrupts variant calls downstream — there is no acceptable error rate.**

**Every code change should be evaluated against four constraints, in order:**

1. **Correctness.** Validation must pass 100%. Any `Invalid > 0` in `validate_gaf` output is a bug. See `VALIDATION_GUIDE.md` for the full validation procedure and common pitfalls.
2. **Performance.** Wall time, peak RSS, cache behavior, file I/O. If a change introduces a per-MEM allocation, an extra hash on the hot path, or doubles a working-set size, justify it. Prefer SDSL-backed bitvectors over `std::vector<bool>`, prefer `int_vector<W>` over `vector<uint64_t>` when bit-width is bounded, prefer streaming over slurp-then-process. Measure before you optimize — use the perf harness in `mem-projection/pangenome-pipeline/perf/`.
3. **Minimality.** No new fields "in case we need them later." No commented-out code blocks. No copy-paste of an existing function with one line changed — refactor or templatize. Reject code bloat aggressively; small surface area is easier to audit for correctness and performance both.
4. **Elegance.** Readable, idiomatic C++17. Names match the paper's terminology where applicable (`bwt_intervals`, `tag_runs`, `pos_t`, etc.). One concept per class. Prefer composition over inheritance. RAII for resource ownership.

**When these constraints conflict, correctness always wins. Then performance wins on the hot path; elegance wins on the cold path.** The hot path is: anything inside `find_mems` per-MEM processing (~93% of time is in `FastLocate::locateNext`), anything inside a `TagArray`/`LightTagIndex` query, the inner loops of `algorithm.hpp` build phases. The cold path is: CLI parsing, file headers, error reporting, debug-stats output.

**Document the why, not the what.** Comments explain non-obvious tradeoffs ("this allocation is amortized across all reads, so a vector is fine here") not what the line obviously does. The 511-cap explanation in `tag_arrays.cpp` is a good model: it answers a question a future reader would actually ask.

## What this directory is

C++17 toolkit (namespace `panindexer`) implementing tag-array-backed pangenome indexing per Eskandar, Paten, Sirén 2025 (*Lossless Pangenome Indexing Using Tag Arrays*). Builds:

- An r-index (run-length BWT + samples) for pattern counting and locate.
- A tag array mapping each BWT position to a graph position `(node_id, offset, strand)`.
- A lightweight tag-array variant (`LightTagIndex`, `.ltags` files) for tag-run-only queries used by find_mems → gafpack → cosigt.

Primary consumers: `find_mems` (MEM search + per-MEM graph projection) and `query_tags`/`query_sampled_tags`. Downstream: gafpack consumes find_mems output to produce per-node coverage vectors.

## Layout

```
include/pangenome_index/   public headers (all classes namespace panindexer)
  tag_arrays.hpp           TagArray: full (node, offset, strand) per run
  light_tag_index.hpp      LightTagIndex: bwt_intervals only; .ltags format
  sampled_tag_array.hpp    Wavelet-tree variant (drops offset)
  r-index.hpp              FastLocate: count, locate, locateNext
  algorithm.hpp            build-time MEM finder, BFS tag extension
  bplus_tree.hpp           BplusTree<Run>: in-memory RLE store, build-only

src/                       one .cpp per binary in PROGRAMS (see makefile:64)
                           plus library .cpp files in LIBOBJS (makefile:61)
deps/grlBWT/               vendored grlBWT (submodule)
bin/                       built binaries (gitignored)
docs/                      tool reference; user-facing docs
tests/                     C++ unit tests + reference outputs
```

## Build

```bash
make                       # all binaries into bin/ ; depends on sdsl-lite at ../sdsl-lite
make bin/find_mems         # build a specific tool
make clean
```

macOS: needs `libomp` + `gnu-time` + zstd via Homebrew (see makefile lines 25–57). The Makefile auto-detects Apple Clang and adapts `PARALLEL_FLAGS`.

## Hot paths (where performance matters most)

- `find_mems` per-MEM processing loop: `src/find_mems.cpp:dump_mem_info_unique_runs` (full tags) and `:dump_mem_info_lightweight` (lightweight). The locate-ops chain dominates (~93% of MEM-processing time).
- `TagArray::query_compressed_decoded_runs` / `LightTagIndex::run_id_at`: called once per MEM.
- `FastLocate::locateNext`: called millions of times per query workload. **This is THE hot function.**
- `algorithm.hpp::parallel_kmers_to_bplustree`, `extend_kmers_bfs_parallel`, `traverse_sequences_parallel`: build-time, but multi-hour on whole-genome pangenomes.

**Performance baselines** (see `VALIDATION_GUIDE.md` for full tables):
- HPRC chr6 (100K reads): find_mems ~48s, ~4.6GB RSS
- Yeast-235 chrII: lightweight is ~5% faster and ~55% less RSS than full-tag mode

Before changing anything in these paths, run the validation query and compare timing:

```sh
cd ../mem-projection/pangenome-pipeline
./query.sh configs/hprcv1-chr6-alt-reads.env hprc-chr6-2026-06-02 <your-tag> --gaf
cat runs/hprc-chr6-2026-06-02/queries/<your-tag>/lightweight/logs/timing_summary.txt
```

Any change that regresses end-to-end wall by >5% or peak RSS by >10% needs justification.

## Cold paths (where elegance / clarity wins)

CLI argv parsing, usage messages, debug-stats printing, file format headers, build-time logging. Spend lines on clarity here, not micro-optimization.

## Doing tasks

- **Plan in TodoWrite** for any task >3 steps. Don't batch completions.
- **Measure before optimizing.** The N=5 perf harness in `mem-projection/pangenome-pipeline/perf/` is the source of truth. Single-trial smoke tests lie; warm-cache N≥3 with `gtime -v` is the contract.
- **One commit per logical change.** Don't bundle "add feature X" with "rename Y" with "fix unrelated typo." The recent `lightweight-tags` branch landed in 3 commits (instrumentation → index + tool → find_mems integration) for exactly this reason.
- **Don't break the default path.** New behavior goes behind flags (`--lightweight-tags`, `--debug-stats`, etc.). Default behavior should be byte-identical to the previous build unless the change is an explicit bugfix.
- **Reuse existing iterators.** `TagArray::for_each_run_compact` exists for a reason. Don't reach into private members.

## Verification before claiming done

**See `VALIDATION_GUIDE.md` for detailed procedures.** Quick reference:

| Change | Verification |
|---|---|
| Anything touching find_mems output | `./query.sh ... --gaf` must show 100% valid. Also compare GAF line counts to catch silently dropped MEMs. |
| Anything in TagArray serialization | `convert_tags` round-trip, check `bwt_intervals size == .seq size + 1` in logs |
| Anything in the build path | Full pipeline rebuild from `.gbz`, compare against existing runs |
| Anything in a hot path | Compare timing against baseline; reject if wall regresses >5% or RSS >10% |
| Anything user-facing (CLI, file format) | Update `docs/tools/<tool>.md` |

## Do not

- Add per-MEM heap allocations on the hot path. Reuse buffers. A `vector` allocation inside the MEM loop will tank performance.
- Introduce new dependencies without confirming with the maintainer. The toolchain is already complex (sdsl-lite, gbwt, gbwtgraph, grlBWT, OpenMP, zstd).
- Add code that's "obviously useful" but doesn't have a current caller. YAGNI applies; dead code is a maintenance tax.
- Skip validation because "it's a small change." Off-by-one errors and node ID mismatches are silent until they corrupt downstream results.
- Store `is_rev` (strand) in the v2 record format. Earlier drafts did this incorrectly — GAF strand is derived from `(seq_id & 1) XOR (GFA step orientation)` in gafpack.
- Edit downstream pipeline configs/docs from here. Open an issue in `mem-projection/pangenome-pipeline` instead.

## Where to read first

1. **`VALIDATION_GUIDE.md`** (this repo): How to validate your changes, common pitfalls, debugging procedures. **Read this before making changes.**
2. **`mem_projection_1.pdf`**: Paper draft with the formal definitions (`Tag[i] = (v, o, b)`, tag-run structure, normalization invariant). The terminology used in code matches this.
3. **`ALGORITHM_CHANGES_SUMMARY.md`**: High-level diagrams of the v2 record format and the unique-tag-runs algorithm.
4. **`mem-projection/pangenome-pipeline/CLAUDE.md`**: Pipeline-level documentation, correctness criteria, known footguns.
