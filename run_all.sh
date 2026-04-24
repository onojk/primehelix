#!/usr/bin/env bash
set -e

echo "== Exact counts (<= 1e8) =="
python exact_lopsided_counts.py

echo
echo "== Convergence plot =="
python plot_lopsided_convergence.py

echo
echo "== Delta shape analysis =="
python analyze_delta_shape.py

echo
echo "Done."
