#!/usr/bin/env bash
# Run build_tag_head_samples against the HPRC chr6 index on Vesuvio.
#
# Prerequisites:
#   - This branch built: cd ~/pangenome-index-latest && make bin/build_tag_head_samples
#   - The .ri under $DATA must be post-JR-014 (has HAS_RUN_HEAD_C flag).
#     The current Vesuvio index is expected to be post-JR-014 already; if the
#     binary rejects the .ri with a "pre-JR-014 format" error, rebuild it:
#       $REPO/bin/build_rindex $DATA/hprcv1_chr6.rl_bwt > $DATA/hprcv1_chr6.ri
#
# Usage: bash run_hprc_chr6.sh
set -euo pipefail

REPO="$HOME/pangenome-index-latest"
DATA="$HOME/mem-projection/pangenome-pipeline/runs/hprc-chr6-2026-06-02"
OUT="$DATA/tag_head_samples"

mkdir -p "$OUT"

RI="$DATA/hprcv1_chr6.ri"
LTAGS="$DATA/hprcv1_chr6.ltags"

# Delete stale v1 (BWT-order) output files so they don't confuse anyone.
# v2 files are the same name; overwrite is fine but explicit cleanup makes
# it obvious in the log that we're starting fresh.
if ls "$OUT"/hprcv1_chr6.tag_samples.s* >/dev/null 2>&1; then
    echo "Removing stale output files (probably v1 BWT-order-rule)..."
    rm -v "$OUT"/hprcv1_chr6.tag_samples.s*
fi

# Log output goes both to stdout and to a file for later reference.
LOG="$OUT/build_tag_head_samples.log"
echo "Running build_tag_head_samples..." | tee "$LOG"
echo "  ri:    $RI"    | tee -a "$LOG"
echo "  ltags: $LTAGS" | tee -a "$LOG"
echo "  out:   $OUT/hprcv1_chr6" | tee -a "$LOG"

"$HOME/.guix-profile/bin/time" -v "$REPO/bin/build_tag_head_samples" \
    "$RI" "$LTAGS" "$OUT/hprcv1_chr6" \
    --s 32,64,128,256 \
    2>&1 | tee -a "$LOG"

echo
echo "Output files:"
ls -la "$OUT"/hprcv1_chr6.tag_samples.s* | tee -a "$LOG"
