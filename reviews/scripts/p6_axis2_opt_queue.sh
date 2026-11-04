#!/bin/sh
# Phase 6 axis 2: the BoTorch runs of Procedure 5, one after another (one process; the campaign uses the other four).
#   sh reviews/scripts/p6_axis2_opt_queue.sh
cd "$(dirname "$0")/../.." || exit 1
for args in "liu F_liu20 4" "liu F_liu20 16" "liu F_wg6 31" "p2 F_liu20 31" "p3 F_liu20 31"; do
  set -- $args
  [ -f "runs/p6_axis2/opt_$1_$2_k$3.json" ] && continue
  OMP_NUM_THREADS=1 uv run python reviews/scripts/p6_axis2_optimize.py --case "$1" --base "$2" --k "$3" \
    > "runs/p6_axis2/opt_$1_$2_k$3.log" 2>&1
done
