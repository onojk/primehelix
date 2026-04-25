#!/usr/bin/env python3
"""
verify_semiprime_tables.py

Six-task verification of paper results using an exact SPF (smallest prime
factor) sieve — no pair-counting, no bisect.  Every semiprime n = p·q
(p ≤ q, p = spf[n]) is enumerated directly.

Tasks
-----
1. Table 1: P(N, θ₀) proportions at N ∈ {1e4,1e5,1e6,1e7,1e8}
             θ₀ ∈ {0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45}
2. Table 2: Local-linear regression  P(N,θ₀) = a(θ₀) − b(θ₀)·log N
             over N ∈ {1e5,…,1e8}  (four points)
3. Table 3: b_eff = δ·log log N,  ratio to f(θ₀)
4. Activation-threshold plot: P vs N at θ₀ = 0.30, 400 log-spaced points
5. Extended Table 2: all θ₀ with R² values
6. Stabilisation test: max |P(N,θ₀) − P(N/2,θ₀)| across θ₀ at N=1e8

Validation checkpoints (exact SPF sieve counts)
------------------------------------------------
  N = 1e6  → 210,035 semiprimes
  N = 1e7  → 1,904,324 semiprimes
  N = 1e8  → 17,427,258 semiprimes

Output
------
  results/task1_table1.csv
  results/task2_table2.csv
  results/task3_table3.csv
  results/task5_extended_table2.csv
  results/task6_stabilisation.csv
  activation_spf_theta0.30.png
"""

import os
import csv
import math
import time
import pickle
import pathlib

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Constants ────────────────────────────────────────────────────────────────

THETA_VALS  = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45]
N_VALS      = [1e4, 1e5, 1e6, 1e7, 1e8]
N_REG       = [1e5, 1e6, 1e7, 1e8]   # regression points (Table 2)
CACHE_DIR   = pathlib.Path(".cache")
RESULTS_DIR = pathlib.Path("results")

GAMMA      = 0.5772156649015328
F_THETA    = {θ: math.log(1/θ) - θ - math.log(2) + 0.5 for θ in THETA_VALS}

CHECKPOINTS = {1_000_000: 210_035, 10_000_000: 1_904_324, 100_000_000: 17_427_258}

# ── SPF sieve ────────────────────────────────────────────────────────────────

def build_spf(N: int) -> np.ndarray:
    """Return smallest-prime-factor array for [0..N], dtype uint32."""
    CACHE_DIR.mkdir(exist_ok=True)
    path = CACHE_DIR / f"spf_{N}.pkl"
    if path.exists():
        print(f"  Loading SPF sieve (N={N:,}) from cache …")
        with open(path, "rb") as f:
            return pickle.load(f)

    print(f"  Building SPF sieve to N={N:,} …", end=" ", flush=True)
    t0 = time.time()
    spf = np.arange(N + 1, dtype=np.uint32)
    for p in range(2, int(N**0.5) + 1):
        if spf[p] == p:                          # p is prime
            composite = np.arange(p * p, N + 1, p, dtype=np.uint32)
            not_set   = spf[p * p :: p] == composite
            spf[p * p :: p][not_set] = p
    print(f"{time.time()-t0:.1f}s")

    with open(path, "wb") as f:
        pickle.dump(spf, f)
    return spf


def extract_semiprimes(spf: np.ndarray, N: int):
    """
    Return (ns_sp, theta) arrays for all semiprimes n ≤ N.

    n is semiprime iff spf[n//spf[n]] == n//spf[n]  (i.e. q = n/p is prime).
    θ(n) = log(p) / log(n)  where p = spf[n].
    """
    CACHE_DIR.mkdir(exist_ok=True)
    path = CACHE_DIR / f"semiprimes_{N}.pkl"
    if path.exists():
        print(f"  Loading semiprime arrays (N={N:,}) from cache …")
        with open(path, "rb") as f:
            return pickle.load(f)

    print(f"  Extracting semiprimes (N={N:,}) …", end=" ", flush=True)
    t0 = time.time()

    # candidates: composite n ≥ 4 with spf[n] < n  (i.e. not prime, not 1)
    n_arr   = np.arange(4, N + 1, dtype=np.uint32)
    p_arr   = spf[4 : N + 1]                        # spf[n]
    q_arr   = (n_arr // p_arr).astype(np.uint32)     # n / spf[n]

    is_sp   = (q_arr >= p_arr) & (spf[q_arr] == q_arr)
    ns_sp   = n_arr[is_sp].astype(np.uint64)
    ps_sp   = p_arr[is_sp].astype(np.float64)
    ns_f    = ns_sp.astype(np.float64)

    theta   = (np.log(ps_sp) / np.log(ns_f)).astype(np.float32)

    print(f"{len(ns_sp):,} semiprimes  ({time.time()-t0:.1f}s)")

    result = (ns_sp, theta)
    with open(path, "wb") as f:
        pickle.dump(result, f)
    return result

# ── P(N, θ₀) at a single N ───────────────────────────────────────────────────

def compute_P(ns_sp: np.ndarray, theta: np.ndarray, N: int, theta0: float) -> float:
    idx = int(np.searchsorted(ns_sp, N, side="right"))
    if idx == 0:
        return 0.0
    return float(np.sum(theta[:idx] <= np.float32(theta0))) / idx


def compute_P_grid(ns_sp: np.ndarray, theta: np.ndarray,
                   N_grid: np.ndarray, theta0: float) -> np.ndarray:
    """Vectorised P(N, θ₀) over a sorted N_grid using cumulative sums."""
    is_lop  = theta <= np.float32(theta0)
    cum_lop = np.cumsum(is_lop, dtype=np.int64)
    cum_tot = np.arange(1, len(theta) + 1, dtype=np.int64)

    P_vals = np.empty(len(N_grid))
    for j, N in enumerate(N_grid):
        idx = int(np.searchsorted(ns_sp, int(N), side="right"))
        P_vals[j] = float(cum_lop[idx - 1]) / idx if idx else 0.0
    return P_vals

# ── Validation ───────────────────────────────────────────────────────────────

def validate_checkpoints(ns_sp: np.ndarray, N_MAX: int):
    print("\nValidation checkpoints:")
    all_ok = True
    for N_chk, expected in sorted(CHECKPOINTS.items()):
        if N_chk > N_MAX:
            print(f"  N={N_chk:>12,}  expected {expected:>12,}  SKIP (beyond sieve limit)")
            continue
        count = int(np.searchsorted(ns_sp, N_chk, side="right"))
        status = "OK" if count == expected else f"FAIL (got {count:,})"
        print(f"  N={N_chk:>12,}  expected {expected:>12,}  got {count:>12,}  {status}")
        if count != expected:
            all_ok = False
    if not all_ok:
        raise RuntimeError("Semiprime count mismatch — check SPF sieve logic.")
    print("  All checkpoints passed.\n")

# ── Task 1: Table 1 ──────────────────────────────────────────────────────────

def task1_table1(ns_sp, theta):
    print("=" * 70)
    print("Task 1: Table 1 — P(N, θ₀)")
    print("=" * 70)

    header = ["N"] + [f"θ={t}" for t in THETA_VALS]
    rows   = []

    print(f"  {'N':>8}", end="")
    for t in THETA_VALS:
        print(f"  θ={t:.2f}", end="")
    print()
    print("  " + "-" * 65)

    for N in N_VALS:
        row = {"N": int(N)}
        print(f"  {N:>8.0e}", end="")
        for t in THETA_VALS:
            P = compute_P(ns_sp, theta, int(N), t)
            row[f"theta={t}"] = round(P, 5)
            print(f"  {P:.5f}", end="")
        print()
        rows.append(row)

    RESULTS_DIR.mkdir(exist_ok=True)
    fname = RESULTS_DIR / "task1_table1.csv"
    with open(fname, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"\n  Saved: {fname}\n")

# ── Task 2 & 5: Table 2 / Extended Table 2 ───────────────────────────────────

def task2_table2(ns_sp, theta):
    print("=" * 70)
    print("Task 2: Table 2 — Local-linear regression δ = a − b·log N")
    print("=" * 70)

    log_N = np.array([math.log(N) for N in N_REG])
    rows  = []

    print(f"  {'θ₀':>5}  {'a':>8}  {'b':>8}  {'R²':>7}")
    print("  " + "-" * 35)

    for t in THETA_VALS:
        P_vals    = np.array([compute_P(ns_sp, theta, int(N), t) for N in N_REG])
        delta_vals = 1.0 - P_vals                    # δ decreases with N
        b, a      = np.polyfit(log_N, delta_vals, 1) # δ = a + b·logN, b < 0
        b_pos     = -b                               # positive magnitude
        yhat      = a + b * log_N
        ss_res    = np.sum((delta_vals - yhat)**2)
        ss_tot    = np.sum((delta_vals - delta_vals.mean())**2)
        r2        = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        print(f"  {t:.2f}  {a:>8.5f}  {b_pos:>8.5f}  {r2:>7.4f}")
        rows.append({"theta0": t, "a": round(a, 6), "b_slope": round(b_pos, 6), "R2": round(r2, 5)})

    RESULTS_DIR.mkdir(exist_ok=True)
    fname = RESULTS_DIR / "task2_table2.csv"
    with open(fname, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"\n  Saved: {fname}\n")
    return rows

# ── Task 3: Table 3 — b_eff ──────────────────────────────────────────────────

def task3_table3(ns_sp, theta, reg_rows):
    print("=" * 70)
    print("Task 3: Table 3 — b_eff = mean(δ·log log N)  and  ratio to f(θ₀)")
    print("=" * 70)

    log_log_N = np.array([math.log(math.log(N)) for N in N_VALS])
    print(f"  log log N range: {log_log_N[0]:.4f} to {log_log_N[-1]:.4f}")
    print(f"  {'θ₀':>5}  {'b_eff':>8}  {'f(θ₀)':>8}  {'ratio':>7}  {'sym_ratio':>10}")
    print("  " + "-" * 48)

    rows = []
    max_sym_ratio = 0.0
    for t in THETA_VALS:
        bvals = []
        for i, N in enumerate(N_VALS):
            P     = compute_P(ns_sp, theta, int(N), t)
            delta = 1.0 - P
            bvals.append(delta * log_log_N[i])
        b_eff  = float(np.mean(bvals))
        f      = F_THETA[t]
        ratio  = b_eff / f if f > 0 else float("nan")
        sym    = ratio if ratio >= 1.0 else (1.0 / ratio)
        max_sym_ratio = max(max_sym_ratio, sym)
        print(f"  {t:.2f}  {b_eff:>8.3f}  {f:>8.3f}  {ratio:>7.3f}  {sym:>10.3f}")
        rows.append({
            "theta0":    t,
            "b_eff":     round(b_eff, 4),
            "f_theta":   round(f, 4),
            "ratio":     round(ratio, 4),
            "sym_ratio": round(sym, 4),
        })

    print(f"\n  Max symmetric ratio = {max_sym_ratio:.4f}  (paper bound: < 1.21×)")
    if max_sym_ratio < 1.21:
        print("  PASS")
    else:
        print("  FAIL — exceeds 1.21× bound")

    fname = RESULTS_DIR / "task3_table3.csv"
    with open(fname, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"\n  Saved: {fname}\n")

# ── Task 4: Activation-threshold plot ────────────────────────────────────────

def task4_activation_plot(ns_sp, theta):
    print("=" * 70)
    print("Task 4: Activation-threshold plot  θ₀ = 0.30")
    print("=" * 70)

    theta0   = 0.30
    N_grid   = np.logspace(4, 8, 400)
    P_vals   = compute_P_grid(ns_sp, theta, N_grid, theta0)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.set_facecolor("white")
    ax.grid(color="#dddddd", lw=0.5)
    ax.semilogx(N_grid, P_vals, color="#2266cc", lw=1.5,
                label=fr"$P(N,\,\theta_0={theta0})$  [SPF sieve]")
    ax.set_xlabel(r"$N$", fontsize=13)
    ax.set_ylabel(r"$P(N,\,\theta_0)$", fontsize=13)
    ax.set_title(fr"Activation thresholds  ($\theta_0 = {theta0}$, SPF sieve)", fontsize=12)
    ax.legend(fontsize=10)
    plt.tight_layout()
    outfile = "activation_spf_theta0.30.png"
    plt.savefig(outfile, dpi=155, facecolor="white", bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outfile}\n")

# ── Task 5: Extended Table 2 (all θ with R²) ─────────────────────────────────

def task5_extended_table2(ns_sp, theta):
    print("=" * 70)
    print("Task 5: Extended Table 2 — full δ regression + R²")
    print("=" * 70)

    log_N = np.array([math.log(N) for N in N_REG])
    rows  = []

    print(f"  {'θ₀':>5}  {'a(θ)':>9}  {'b(θ)':>9}  {'R²':>7}  "
          f"{'f(θ)':>7}  {'b_eff':>7}  {'beff/f':>7}")
    print("  " + "-" * 60)

    log_log_N = np.array([math.log(math.log(N)) for N in N_VALS])

    for t in THETA_VALS:
        P_vals     = np.array([compute_P(ns_sp, theta, int(N), t) for N in N_REG])
        delta_vals = 1.0 - P_vals
        b, a       = np.polyfit(log_N, delta_vals, 1)  # δ = a + b·logN, b < 0
        b_pos      = -b
        yhat       = a + b * log_N
        ss_res     = np.sum((delta_vals - yhat)**2)
        ss_tot     = np.sum((delta_vals - delta_vals.mean())**2)
        r2         = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        f          = F_THETA[t]
        bvals      = [(1 - compute_P(ns_sp, theta, int(N), t)) * log_log_N[i]
                      for i, N in enumerate(N_VALS)]
        b_eff      = float(np.mean(bvals))
        ratio_eff  = b_eff / f if f > 0 else float("nan")
        print(f"  {t:.2f}  {a:>9.5f}  {b_pos:>9.5f}  {r2:>7.4f}  {f:>7.5f}  "
              f"{b_eff:>7.3f}  {ratio_eff:>7.4f}")
        rows.append({
            "theta0":    t,
            "a":         round(a, 6),
            "b":         round(b_pos, 6),
            "R2":        round(r2, 5),
            "f_theta":   round(f, 6),
            "b_eff":     round(b_eff, 4),
            "beff_over_f": round(ratio_eff, 5),
        })

    fname = RESULTS_DIR / "task5_extended_table2.csv"
    with open(fname, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"\n  Saved: {fname}\n")

# ── Task 6: Stabilisation test ───────────────────────────────────────────────

def task6_stabilisation(ns_sp, theta):
    print("=" * 70)
    print("Task 6: Stabilisation test — max |P(N,θ) − P(N/2,θ)| at N=1e8")
    print("=" * 70)

    N_hi = int(1e8)
    N_lo = N_hi // 2

    rows = []
    diffs = []
    print(f"  {'θ₀':>5}  {'P(N/2)':>9}  {'P(N)':>9}  {'|diff|':>9}")
    print("  " + "-" * 40)

    for t in THETA_VALS:
        P_hi = compute_P(ns_sp, theta, N_hi, t)
        P_lo = compute_P(ns_sp, theta, N_lo, t)
        diff = abs(P_hi - P_lo)
        diffs.append(diff)
        print(f"  {t:.2f}  {P_lo:>9.6f}  {P_hi:>9.6f}  {diff:>9.6f}")
        rows.append({"theta0": t, "P_N_half": round(P_lo, 7), "P_N": round(P_hi, 7),
                     "abs_diff": round(diff, 7)})

    max_diff = max(diffs)
    print(f"\n  Max |diff| = {max_diff:.6f}  (expect < 0.01 for stabilisation)")
    if max_diff < 0.01:
        print("  PASS — P has stabilised")
    else:
        print("  NOTE — P still changing at rate > 0.01")

    fname = RESULTS_DIR / "task6_stabilisation.csv"
    with open(fname, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"\n  Saved: {fname}\n")

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    N_MAX = int(max(N_VALS))

    print("=" * 70)
    print("verify_semiprime_tables.py")
    print("SPF-sieve-based verification of paper results")
    print("=" * 70)
    print()

    t_start = time.time()

    spf     = build_spf(N_MAX)
    ns_sp, theta = extract_semiprimes(spf, N_MAX)

    validate_checkpoints(ns_sp, N_MAX)

    RESULTS_DIR.mkdir(exist_ok=True)

    task1_table1(ns_sp, theta)
    reg_rows = task2_table2(ns_sp, theta)
    task3_table3(ns_sp, theta, reg_rows)
    task4_activation_plot(ns_sp, theta)
    task5_extended_table2(ns_sp, theta)
    task6_stabilisation(ns_sp, theta)

    print("=" * 70)
    print(f"All tasks complete.  Total wall time: {time.time()-t_start:.1f}s")
    print("=" * 70)


if __name__ == "__main__":
    main()
