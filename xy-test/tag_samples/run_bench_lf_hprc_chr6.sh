#!/usr/bin/env bash
# Run bench_lf against the HPRC chr6 .ri on Vesuvio.
#
# Usage: bash run_bench_lf_hprc_chr6.sh
set -euo pipefail

REPO="$HOME/pangenome-index-latest"
DATA="$HOME/mem-projection/pangenome-pipeline/runs/hprc-chr6-2026-06-02"
OUT="$DATA/tag_head_samples"

mkdir -p "$OUT"

RI="$DATA/hprcv1_chr6.ri"
LOG="$OUT/bench_lf.log"

echo "Running bench_lf..." | tee "$LOG"
echo "  ri: $RI" | tee -a "$LOG"

# Prefer gtime (Guix) with -v; fall back to /usr/bin/time -v if gtime not
# available (matches JR-013's pattern for time-tooling on Vesuvio).
TIMEBIN=""
if command -v gtime >/dev/null 2>&1; then TIMEBIN="gtime -v"
elif [ -x /usr/bin/time ]; then TIMEBIN="/usr/bin/time -v"
fi

# Defaults tuned for HPRC chr6 -- larger N to average out cache noise on the
# 3.6 GB .ri. Chain uses 10k starts x 100 steps = 1M chain ops, matching
# random-mode N.
${TIMEBIN} "$REPO/bin/bench_lf" "$RI" \
    --n-random 1000000 \
    --m-starts 10000 \
    --k-steps  100 \
    --trials   5 \
    --verify   10000 \
    2>&1 | tee -a "$LOG"

echo
echo "Log written to: $LOG"
