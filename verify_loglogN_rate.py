#!/usr/bin/env python3
"""
verify_loglogN_rate.py

For each theta_0 in {0.20, 0.25, 0.30, 0.35}, fit two models to the
deficit delta(N, theta_0) = 1 - P(N, theta_0) over N in {1e4 .. 1e8}:

  Model A: delta = a - b * log(N)          (local linear, original fit)
  Model B: delta = c / log(log(N))         (theoretical rate,
                                             Montgomery-Vaughan §7.4)

Reports R² for both, plots both fits, saves loglogN_fit.png.
"""
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pickle, pathlib
from bisect import bisect_left, bisect_right
from math import isqrt, ceil

# ── Data from exact enumeration (SPF sieve, no sampling) ─────────────────────
# P(N, theta_0) in percent — reproduced from exact_lopsided_counts.py output
# and plot_multi_theta.py.  Rows: theta = 0.20, 0.25, 0.30, 0.35
# Columns: N = 1e4, 1e5, 1e6, 1e7, 1e8

P_PCT = {
    0.20: [49.524, 51.933, 53.550, 55.533, 56.294],
    0.25: [59.276, 60.771, 62.468, 63.741, 64.947],
    0.30: [67.733, 69.403, 70.509, 71.825, 72.741],
    0.35: [76.457, 77.274, 78.392, 79.189, 79.950],
}

THETAS = [0.20, 0.25, 0.30, 0.35]
N_VALS = [1e4, 1e5, 1e6, 1e7, 1e8]

# ── Helper ────────────────────────────────────────────────────────────────────

def r2(y, yhat):
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0


def fit_model_A(log_n, delta):
    """delta = a - b*log(N); returns (a, b, R²)."""
    b_neg, a = np.polyfit(log_n, delta, 1)
    b = -b_neg
    yhat = a - b * log_n
    return a, b, r2(delta, yhat)


def fit_model_B(llN, delta):
    """delta = c / log(log(N)); 1-parameter fit via weighted least squares.
    Minimise ||delta - c/llN||² → c = sum(delta/llN) / sum(1/llN²)."""
    inv = 1.0 / llN
    c = np.dot(delta, inv) / np.dot(inv, inv)
    yhat = c / llN
    return c, r2(delta, yhat)

# ── Compute ───────────────────────────────────────────────────────────────────

log_n = np.array([math.log(n) for n in N_VALS])
llN   = np.array([math.log(math.log(n)) for n in N_VALS])

print("log(log(N)) values over the observed range:")
for n, ll in zip(N_VALS, llN):
    print(f"  N = {n:.0e}   log(N) = {math.log(n):.3f}   log(log(N)) = {ll:.3f}"
          f"   1/log(log(N)) = {1/ll:.4f}")

print()
print(f"Range of log(log(N)):  {llN.min():.3f} – {llN.max():.3f}  "
      f"(variation = {llN.max()-llN.min():.3f})")
print(f"Range of log(N):       {log_n.min():.3f} – {log_n.max():.3f}  "
      f"(variation = {log_n.max()-log_n.min():.3f})")
print()

print(f"{'theta':>6}  {'a':>8}  {'b':>8}  {'R²_A':>8}  "
      f"{'c':>8}  {'R²_B':>8}  {'winner'}")
print("─" * 72)

results = {}
for theta in THETAS:
    delta = 1.0 - np.array(P_PCT[theta]) / 100.0
    a, b, r2_A = fit_model_A(log_n, delta)
    c, r2_B    = fit_model_B(llN,   delta)
    winner = "Model A" if r2_A > r2_B else "Model B"
    diff   = abs(r2_A - r2_B)
    results[theta] = (a, b, r2_A, c, r2_B, winner, diff, delta)
    print(f"{theta:>6.2f}  {a:>8.5f}  {b:>8.6f}  {r2_A:>8.6f}  "
          f"{c:>8.5f}  {r2_B:>8.6f}  {winner} by {diff:.6f}")

print()
print("Key finding:")
r2_A_mean = np.mean([results[t][2] for t in THETAS])
r2_B_mean = np.mean([results[t][4] for t in THETAS])
print(f"  Mean R² Model A (linear in log N):    {r2_A_mean:.6f}")
print(f"  Mean R² Model B (c/log(log(N))):      {r2_B_mean:.6f}")
print()
if r2_A_mean > r2_B_mean:
    print("  Model A fits better within [1e4, 1e8] — consistent with the")
    print("  explanation that log(log(N)) varies by only "
          f"{llN.max()-llN.min():.3f} across this range,")
    print("  making it nearly linear in log(N) here.")
else:
    print("  Model B fits better — 1/log(log(N)) decay directly confirmed.")

# ── Plot ──────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(2, 2, figsize=(11, 8))
fig.suptitle(
    r"Deficit $\delta(N,\theta_0)$: Model A ($a - b\log N$) vs "
    r"Model B ($c/\log\log N$)",
    fontsize=12
)

N_smooth = np.logspace(4, 8, 400)
log_smooth = np.log(N_smooth)
llN_smooth = np.log(np.log(N_smooth))

for ax, theta in zip(axes.flat, THETAS):
    a, b, r2_A, c, r2_B, winner, diff, delta = results[theta]

    ax.scatter(N_VALS, delta, color="black", zorder=5, s=40, label="Data")
    ax.plot(N_smooth, a - b * log_smooth,
            color="#2166ac", lw=1.8,
            label=fr"A: $a-b\log N$  $R^2={r2_A:.4f}$")
    ax.plot(N_smooth, c / llN_smooth,
            color="#d6604d", lw=1.8, ls="--",
            label=fr"B: $c/\log\log N$  $R^2={r2_B:.4f}$")

    ax.set_xscale("log")
    ax.set_title(fr"$\theta_0 = {theta}$", fontsize=11)
    ax.set_xlabel(r"$N$", fontsize=10)
    ax.set_ylabel(r"$\delta(N,\theta_0)$", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(color="#dddddd", lw=0.5)
    ax.tick_params(labelsize=8)

plt.tight_layout()
plt.savefig("loglogN_fit.png", dpi=150, facecolor="white", bbox_inches="tight")
plt.close()
print("Saved: loglogN_fit.png")
