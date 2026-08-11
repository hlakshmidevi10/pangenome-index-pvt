#!/usr/bin/env bash
# Build tag_head_samples at s=8 and s=16 for HPRC chr6, to close the
# s-curve investigation started in JR-016/JR-017/JR-018.
#
# Standalone from run_hprc_chr6.sh (which built s=32/64/128/256):
#   - Preserves the existing .tag_samples.s{32,64,128,256} files.
#   - Writes a distinct log file so JR-016 build log is preserved.
#   - Runs a single invocation for BOTH s values (Phase 1 SA enumeration
#     is per-invocation, ~82 min; amortizing across s=8 and s=16 saves
#     ~80 min vs two separate builds).
#
# Prerequisites:
#   - Branch tag-head-samples pulled, bin/build_tag_head_samples built.
#   - .ri is post-JR-014 (HAS_RUN_HEAD_C flag).
#
# Wall estimate (per JR-016 timing table):
#   Phase 1 (SA enum): ~82 min
#   Sort:              ~1 min
#   Per s pass:        ~4-8s each (higher s than 32 built ~4s; s=8/16 may
#                      take slightly longer due to more kept candidates)
#   Total wall:        ~85 min
#
# Peak RSS estimate: ~22 GB (dominated by sorted candidate array, same as
# JR-016 build; s=8 and s=16 add a small amount over s=32 for the extra
# kept sd_vector + int_vector but that's O(100 MB), noise-level vs 22 GB).

set -euo pipefail

REPO="$HOME/pangenome-index-latest"
DATA="$HOME/mem-projection/pangenome-pipeline/runs/hprc-chr6-2026-06-02"
OUT="$DATA/tag_head_samples"

mkdir -p "$OUT"

RI="$DATA/hprcv1_chr6.ri"
LTAGS="$DATA/hprcv1_chr6.ltags"

# Sanity: refuse to run if s=32/64/128/256 don't already exist. If they're
# missing something changed since JR-016 and we want to fail fast rather
# than silently do a full rebuild.
for s in 32 64 128 256; do
    f="$OUT/hprcv1_chr6.tag_samples.s${s}"
    if [ ! -f "$f" ]; then
        echo "ERROR: expected existing file not found: $f"
        echo "Refusing to run; investigate before rebuilding."
        exit 1
    fi
done
echo "Existing s=32/64/128/256 files present; leaving untouched."

# Refuse to overwrite existing s=8 or s=16 outputs; force explicit rm first.
for s in 8 16; do
    f="$OUT/hprcv1_chr6.tag_samples.s${s}"
    if [ -f "$f" ]; then
        echo "ERROR: $f already exists. Delete it explicitly before rerunning."
        exit 1
    fi
done

LOG="$OUT/build_tag_head_samples_s8_s16.log"
echo "Running build_tag_head_samples for s=8,16..." | tee "$LOG"
echo "  ri:    $RI"    | tee -a "$LOG"
echo "  ltags: $LTAGS" | tee -a "$LOG"
echo "  out:   $OUT/hprcv1_chr6" | tee -a "$LOG"
echo "  s:     8, 16"  | tee -a "$LOG"
echo "  log:   $LOG"   | tee -a "$LOG"

"$HOME/.guix-profile/bin/time" -v "$REPO/bin/build_tag_head_samples" \
    "$RI" "$LTAGS" "$OUT/hprcv1_chr6" \
    --s 8,16 \
    2>&1 | tee -a "$LOG"

echo
echo "New output files (should be exactly two: s8, s16):"
ls -la "$OUT"/hprcv1_chr6.tag_samples.s8 "$OUT"/hprcv1_chr6.tag_samples.s16 | tee -a "$LOG"

echo
echo "Full s-curve now available:"
ls -la "$OUT"/hprcv1_chr6.tag_samples.s* | tee -a "$LOG"
