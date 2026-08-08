#!/usr/bin/env bash
# Run build_tag_head_samples against the HPRC chr6 index on Vesuvio.
#
# Prerequisites on Vesuvio:
#   - This branch built: cd ~/pangenome-index-latest && make bin/build_tag_head_samples
#   - Post-JR-014 r-index: needs to be regenerated from .rl_bwt if the on-disk
#     .ri predates JR-014. Check the load error to know.
#
# Usage: bash run_hprc_chr6.sh
set -euo pipefail

REPO="$HOME/pangenome-index-latest"
DATA="$HOME/mem-projection/pangenome-pipeline/runs/hprc-chr6-2026-06-02"
OUT="$DATA/tag_head_samples"

mkdir -p "$OUT"

# Rebuild .ri from .rl_bwt if we have the .rl_bwt and the current .ri is stale.
# On Vesuvio the .rl_bwt should exist (it wasn't shipped to Mac).
if [ -f "$DATA/hprcv1_chr6.rl_bwt" ]; then
    if [ ! -f "$OUT/hprcv1_chr6.ri" ] || [ "$DATA/hprcv1_chr6.rl_bwt" -nt "$OUT/hprcv1_chr6.ri" ]; then
        echo "Rebuilding r-index from .rl_bwt (post-JR-014 format)..."
        "$REPO/bin/build_rindex" "$DATA/hprcv1_chr6.rl_bwt" > "$OUT/hprcv1_chr6.ri" 2> "$OUT/build_rindex.log"
        echo "Done. build_rindex.log tail:"
        tail -5 "$OUT/build_rindex.log"
    fi
    RI="$OUT/hprcv1_chr6.ri"
else
    echo "No .rl_bwt on this host; using existing .ri (may need regeneration)."
    RI="$DATA/hprcv1_chr6.ri"
fi

LTAGS="$DATA/hprcv1_chr6.ltags"

# Log output goes both to stdout and to a file for later reference.
LOG="$OUT/build_tag_head_samples.log"
echo "Running build_tag_head_samples..." | tee "$LOG"
echo "  ri:    $RI"    | tee -a "$LOG"
echo "  ltags: $LTAGS" | tee -a "$LOG"
echo "  out:   $OUT/hprcv1_chr6" | tee -a "$LOG"

/usr/bin/time -v "$REPO/bin/build_tag_head_samples" \
    "$RI" "$LTAGS" "$OUT/hprcv1_chr6" \
    --s 32,64,128,256 \
    2>&1 | tee -a "$LOG"

echo
echo "Output files:"
ls -la "$OUT"/hprcv1_chr6.tag_samples.s* | tee -a "$LOG"
