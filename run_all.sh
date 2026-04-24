#!/usr/bin/env bash
set -euo pipefail

echo "=== PRIMEHELIX: FULL ANALYSIS RUN ==="
echo

# Ensure we run from repo root
cd "$(dirname "$0")"

echo "== Exact counts (<= 1e8) =="
python exact_lopsided_counts.py
echo

echo "== Single-θ convergence plot =="
python plot_lopsided_convergence.py
echo

echo "== Multi-θ convergence =="
python plot_multi_theta.py
echo

echo "== Δ(N,θ) normalization =="
python plot_multi_theta_delta.py
echo

echo "== Curve collapse test =="
python plot_theta_curve_collapse.py
echo

echo "== Extended slope fitting (includes 1e3 + 1e9 sampled) =="
python fit_theta_slopes_extended.py
echo

echo "=== DONE ==="
