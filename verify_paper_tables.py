#!/usr/bin/env python3
"""
verify_paper_tables.py

Six tasks against the paper "On the Distribution of Smaller Factors in Semiprimes":
  1. Verify Table 1  — P(N, θ₀) values
  2. Reproduce Table 2 — local linear fit δ ≈ a - b·log N
  3. Reproduce Table 3 — beff vs naive log(1/θ₀) vs heuristic f(θ₀)
  4. Plot P(N, θ₀) at θ₀=0.10 with activation-threshold verticals
  5. Extend Table 2 to all 8 thresholds
  6. Stabilisation test — does δ·log log N flatten as N grows? (Conjecture 1)

All semiprime counts use exact pair enumeration over primes — equivalent to an SPF
sieve but memory-efficient: we iterate over the smaller factor p and count valid
larger factors q via bisect on the pre-built prime list.
"""

import math
import pickle
import pathlib
import time
from bisect import bisect_left, bisect_right
from math import isqrt, ceil, log

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Configuration ─────────────────────────────────────────────────────────────

THETAS_PAPER = [0.20, 0.25, 0.30, 0.35]          # Tables 1 and 2
THETAS_ALL   = [0.10, 0.15, 0.20, 0.25,           # Table 3 and extended Table 2
                0.30, 0.35, 0.40, 0.45]
N_VALS       = [1e4, 1e5, 1e6, 1e7, 1e8]
LOG2         = math.log(2)
CACHE        = pathlib.Path(".cache")

# Paper Table 1 reference values (percent)
TABLE1_REF = {
    0.20: [49.52, 51.93, 53.55, 55.53, 56.29],
    0.25: [59.28, 60.77, 62.47, 63.74, 64.95],
    0.30: [67.73, 69.40, 70.51, 71.83, 72.74],
    0.35: [76.46, 77.27, 78.39, 79.19, 79.95],
}

# Paper Table 2 reference (a, b, R²)
TABLE2_REF = {
    0.20: (0.5692, 0.007444, 0.9754),
    0.25: (0.4635, 0.006216, 0.9955),
    0.30: (0.3702, 0.005402, 0.9904),
    0.35: (0.2709, 0.003866, 0.9954),
}

# Paper Table 3 reference beff values
TABLE3_REF_BEFF = {
    0.10: 1.819, 0.15: 1.465, 0.20: 1.205, 0.25: 0.976,
    0.30: 0.763, 0.35: 0.562, 0.40: 0.373, 0.45: 0.187,
}

# ── Prime loading ──────────────────────────────────────────────────────────────

def load_primes(limit: int) -> list:
    CACHE.mkdir(exist_ok=True)
    path = CACHE / f"primes_{limit}.pkl"
    if path.exists():
        with open(path, "rb") as fh:
            return pickle.load(fh)
    print(f"Sieving primes to {limit:,} …")
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, isqrt(limit) + 1):
        if sieve[i]:
            sieve[i*i : limit+1 : i] = b"\x00" * len(sieve[i*i : limit+1 : i])
    primes = [i for i in range(2, limit + 1) if sieve[i]]
    with open(path, "wb") as fh:
        pickle.dump(primes, fh)
    return primes

# ── Core counting function ─────────────────────────────────────────────────────

def count_P(limit: int, theta: float, primes: list) -> tuple[float, int, int]:
    """
    Returns (P, lopsided, total) where P = proportion of semiprimes n ≤ limit
    with θ(n) ≤ theta.

    Exact enumeration: iterate over smallest prime factor p ≤ √limit, count
    valid larger factors q via bisect.  Equivalent to an SPF sieve on semiprimes.
    For θ ≤ θ₀: q ≥ p^((1-θ₀)/θ₀)  (derived from θ ≤ θ₀ definition).
    """
    total = lopsided = 0
    for p in primes:
        if p * p > limit:
            break
        q_max = limit // p
        # All valid q (primes in [p, q_max])
        cnt_all = bisect_right(primes, q_max) - bisect_left(primes, p)
        total += cnt_all
        # Lopsided q satisfy q ≥ p^((1-θ)/θ)
        q_min = max(p, ceil(p ** ((1.0 / theta) - 1.0)))
        if q_min <= q_max:
            lopsided += bisect_right(primes, q_max) - bisect_left(primes, q_min)
    P = lopsided / total if total else 0.0
    return P, lopsided, total

# ── Helpers ────────────────────────────────────────────────────────────────────

def heuristic_f(theta: float) -> float:
    return log(1.0 / theta) - theta - LOG2 + 0.5

def section(title: str):
    print()
    print("=" * 72)
    print(title)
    print("=" * 72)

# ── Task 1: Verify Table 1 ─────────────────────────────────────────────────────

def task1_verify_table1(primes):
    section("TASK 1 — Verify Table 1: P(N, θ₀) in percent")

    header = f"  {'N':>8}  " + "  ".join(f"θ={t}" for t in THETAS_PAPER)
    print(header)
    print("  " + "-" * 60)

    results = {t: [] for t in THETAS_PAPER}
    all_ok = True

    for n in N_VALS:
        limit = int(n)
        row = f"  {n:>8.0e}  "
        for theta in THETAS_PAPER:
            P, _, _ = count_P(limit, theta, primes)
            pct = P * 100
            results[theta].append(pct)
            row += f"{pct:7.2f}  "
        print(row)

    print()
    print("  Comparison with paper (|diff| ≤ 0.01 pp = PASS):")
    print(f"  {'theta':>6}  {'N':>8}  {'paper':>8}  {'computed':>9}  {'diff':>7}  status")
    print("  " + "-" * 55)
    for theta in THETAS_PAPER:
        for i, n in enumerate(N_VALS):
            ref = TABLE1_REF[theta][i]
            got = results[theta][i]
            diff = got - ref
            ok = abs(diff) <= 0.01
            if not ok:
                all_ok = False
            status = "PASS" if ok else "FAIL ←"
            print(f"  {theta:>6.2f}  {n:>8.0e}  {ref:>8.2f}  {got:>9.2f}  {diff:>+7.3f}  {status}")

    print()
    print(f"  Table 1 overall: {'ALL PASS' if all_ok else 'FAILURES DETECTED'}")
    return results

# ── Task 2: Reproduce Table 2 ─────────────────────────────────────────────────

def task2_table2(table1_results):
    section("TASK 2 — Reproduce Table 2: δ ≈ a(θ₀) − b(θ₀)·log N")

    logN = np.array([log(n) for n in N_VALS])
    print(f"  {'theta':>6}  {'a_paper':>9}  {'a_comp':>9}  "
          f"{'b_paper':>9}  {'b_comp':>9}  {'R²_paper':>9}  {'R²_comp':>9}")
    print("  " + "-" * 70)

    all_ok = True
    for theta in THETAS_PAPER:
        pct_vals  = np.array(table1_results[theta])
        delta_vals = 1.0 - pct_vals / 100.0
        b, a = np.polyfit(logN, delta_vals, 1)   # delta = a + b*logN  (b negative)
        fit   = a + b * logN
        ss_res = np.sum((delta_vals - fit) ** 2)
        ss_tot = np.sum((delta_vals - delta_vals.mean()) ** 2)
        r2    = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0

        a_ref, b_ref, r2_ref = TABLE2_REF[theta]
        # paper uses δ ≈ a - b·logN so their b is positive; polyfit gives b negative
        b_pos = -b

        ok_a  = abs(a - a_ref) < 0.001
        ok_b  = abs(b_pos - b_ref) < 0.0002
        ok_r2 = abs(r2 - r2_ref) < 0.002
        ok    = ok_a and ok_b and ok_r2
        if not ok:
            all_ok = False
        tag = "" if ok else "←"
        print(f"  {theta:>6.2f}  {a_ref:>9.4f}  {a:>9.4f}  "
              f"{b_ref:>9.6f}  {b_pos:>9.6f}  {r2_ref:>9.4f}  {r2:>9.4f}  {tag}")

    print()
    print(f"  Table 2 overall: {'ALL PASS (±tol)' if all_ok else 'FAILURES DETECTED'}")

# ── Task 3: Reproduce Table 3 ─────────────────────────────────────────────────

def task3_table3(primes):
    section("TASK 3 — Reproduce Table 3: beff vs naive vs heuristic f(θ₀)")

    llN = np.array([math.log(math.log(n)) for n in N_VALS])
    print(f"  log log N range: {llN[0]:.4f} to {llN[-1]:.4f} (variation {llN[-1]-llN[0]:.4f})")
    print()

    print(f"  {'θ₀':>5}  {'log(1/θ)':>9}  {'f(θ₀)':>8}  {'beff':>7}  "
          f"{'beff/log':>9}  {'beff/f':>7}  {'paper_beff':>10}  status")
    print("  " + "-" * 74)

    max_ratio_naive    = 0.0
    max_ratio_heuristic = 0.0
    all_ok = True

    for theta in THETAS_ALL:
        bvals = []
        for n in N_VALS:
            P, _, _ = count_P(int(n), theta, primes)
            delta = 1.0 - P
            bvals.append(delta * math.log(math.log(n)))
        b_eff = float(np.mean(bvals))

        log_inv = log(1.0 / theta)
        f_val   = heuristic_f(theta)
        r_naive = b_eff / log_inv if log_inv > 0 else float("nan")
        r_heur  = b_eff / f_val  if f_val  > 0 else float("nan")

        max_ratio_naive     = max(max_ratio_naive,     1.0 / r_naive if r_naive < 1 else r_naive)
        max_ratio_heuristic = max(max_ratio_heuristic, 1.0 / r_heur  if r_heur  < 1 else r_heur)

        ref_beff = TABLE3_REF_BEFF.get(theta)
        ok = ref_beff is None or abs(b_eff - ref_beff) < 0.002
        if not ok:
            all_ok = False
        tag = "" if ok else "←"
        ref_str = f"{ref_beff:10.3f}" if ref_beff is not None else "         —"
        print(f"  {theta:>5.2f}  {log_inv:>9.3f}  {f_val:>8.3f}  {b_eff:>7.3f}  "
              f"{r_naive:>9.3f}  {r_heur:>7.3f}  {ref_str}  {tag}")

    print()
    print(f"  Max relative error — naive:    {max_ratio_naive:.2f}×")
    print(f"  Max relative error — heuristic: {max_ratio_heuristic:.2f}×")
    claim_ok = max_ratio_heuristic < 1.21
    print(f"  Claim 'under 1.21×': {'CONFIRMED' if claim_ok else 'NOT CONFIRMED'}")
    print(f"  Table 3 beff overall: {'ALL PASS' if all_ok else 'FAILURES DETECTED'}")

# ── Task 4: Activation threshold plot ─────────────────────────────────────────

def task4_activation_plot(primes):
    section("TASK 4 — Plot P(N, θ₀=0.10) with activation threshold verticals")

    THETA = 0.10
    # Dense N grid on log scale
    N_grid = np.logspace(3, 8.3, 200)
    P_vals = []
    for n in N_grid:
        P, _, _ = count_P(int(n), THETA, primes)
        P_vals.append(P * 100)

    # Activation thresholds for θ₀=0.10: p^(1/0.10) = p^10
    threshold_primes = [2, 3, 5, 7, 11, 13]
    thresholds = [(p, p ** (1.0 / THETA)) for p in threshold_primes
                  if p ** (1.0 / THETA) <= 3e8]

    fig, ax = plt.subplots(figsize=(10, 5.5))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    ax.semilogx(N_grid, P_vals, color="#1f77b4", lw=2.0, label=r"$P(N,\theta_0=0.10)$")

    colors = plt.cm.tab10(np.linspace(0, 0.7, len(thresholds)))
    for (p, N_thresh), c in zip(thresholds, colors):
        ax.axvline(N_thresh, color=c, lw=1.2, ls="--", alpha=0.85,
                   label=fr"$p={p}$: $N={p}^{{10}}={N_thresh:.2g}$")
        ax.text(N_thresh * 1.06, ax.get_ylim()[0] + 0.5 if ax.get_ylim()[0] > 0 else 1,
                f"$p={p}$", fontsize=8, color=c, va="bottom", rotation=90)

    # Mark the N_VALS used in the paper
    for n in N_VALS:
        P_n, _, _ = count_P(int(n), THETA, primes)
        ax.scatter([n], [P_n * 100], color="black", s=40, zorder=6)

    ax.set_xlabel(r"$N$", fontsize=13)
    ax.set_ylabel(r"$P(N,\,\theta_0=0.10)$ (%)", fontsize=12)
    ax.set_title(
        r"$P(N,\theta_0{=}0.10)$ with prime activation thresholds $N = p^{1/\theta_0} = p^{10}$"
        "\n(dots = paper N values; dashed = threshold where prime p first contributes)",
        fontsize=10,
    )
    ax.legend(fontsize=8, loc="lower right", ncol=2)
    ax.grid(color="#dddddd", lw=0.5)

    out = "activation_thresholds_plot.png"
    plt.tight_layout()
    plt.savefig(out, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")

    print()
    print("  Activation thresholds at θ₀=0.10 (p^10):")
    print(f"  {'p':>4}  {'N = p^10':>14}  in paper range?")
    for p, Nt in thresholds:
        in_range = 1e4 <= Nt <= 1e8
        print(f"  {p:>4}  {Nt:>14,.0f}  {'YES ←' if in_range else 'no'}")

# ── Task 5: Extended Table 2 ───────────────────────────────────────────────────

def task5_extended_table2(primes):
    section("TASK 5 — Extended Table 2: all 8 thresholds")

    logN = np.array([log(n) for n in N_VALS])
    llN  = np.array([math.log(math.log(n)) for n in N_VALS])

    print(f"  {'θ₀':>5}  {'a':>8}  {'b':>9}  {'R²':>7}  "
          f"{'f(θ₀)':>8}  {'beff':>7}  {'note'}")
    print("  " + "-" * 65)

    for theta in THETAS_ALL:
        P_vals    = []
        beff_vals = []
        for i, n in enumerate(N_VALS):
            P, _, _ = count_P(int(n), theta, primes)
            delta = 1.0 - P
            P_vals.append(P)
            beff_vals.append(delta * llN[i])

        delta_arr = 1.0 - np.array(P_vals)
        b, a      = np.polyfit(logN, delta_arr, 1)
        fit       = a + b * logN
        ss_res    = np.sum((delta_arr - fit) ** 2)
        ss_tot    = np.sum((delta_arr - delta_arr.mean()) ** 2)
        r2        = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        b_eff     = float(np.mean(beff_vals))
        f_val     = heuristic_f(theta)
        b_pos     = -b

        note = ""
        if theta in TABLE2_REF:
            a_ref, b_ref, _ = TABLE2_REF[theta]
            note = "(paper)" if abs(a - a_ref) < 0.001 and abs(b_pos - b_ref) < 0.0002 else "(paper — mismatch)"
        print(f"  {theta:>5.2f}  {a:>8.4f}  {b_pos:>9.6f}  {r2:>7.4f}  "
              f"{f_val:>8.3f}  {b_eff:>7.3f}  {note}")

# ── Task 6: Stabilisation test ─────────────────────────────────────────────────

def task6_stabilisation(primes):
    section("TASK 6 — Stabilisation test: δ·log log N vs N (Conjecture 1)")

    # Use a denser N grid on log scale: 12 points from 1e4 to 1e8
    N_dense = np.logspace(4, 8, 12)
    llN_dense = np.array([math.log(math.log(n)) for n in N_dense])

    COLORS = plt.cm.tab10(np.linspace(0, 0.9, len(THETAS_ALL)))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    fig.suptitle(
        r"Stabilisation test: $\delta(N,\theta_0)\cdot\log\log N$ vs $N$"
        "\n(flat → Conjecture 1 holds at this range; rising → subleading corrections)",
        fontsize=11,
    )
    for ax in axes:
        ax.set_facecolor("white")
        ax.grid(color="#dddddd", lw=0.5)

    ax1, ax2 = axes

    print(f"  {'θ₀':>5}  {'f(θ₀)':>8}  " +
          "  ".join(f"N={n:.0e}" for n in [1e4, 1e5, 1e6, 1e7, 1e8]))
    print("  " + "-" * 72)

    for theta, color in zip(THETAS_ALL, COLORS):
        beff_dense = []
        for n in N_dense:
            P, _, _ = count_P(int(n), theta, primes)
            delta = 1.0 - P
            beff_dense.append(delta * math.log(math.log(n)))

        # Print values at the 5 paper N points (indices closest to paper N_VALS)
        paper_idx = [np.argmin(np.abs(N_dense - n)) for n in N_VALS]
        vals_str = "  ".join(f"{beff_dense[i]:7.4f}" for i in paper_idx)
        f_val = heuristic_f(theta)
        print(f"  {theta:>5.2f}  {f_val:>8.4f}  {vals_str}")

        ax1.semilogx(N_dense, beff_dense, color=color, lw=1.8, marker=".", ms=4,
                     label=fr"$\theta_0={theta}$")
        ax1.axhline(f_val, color=color, lw=0.8, ls="--", alpha=0.5)

        # Normalised: beff / f(theta)
        f_val_safe = f_val if abs(f_val) > 1e-9 else 1e-9
        norm = [b / f_val_safe for b in beff_dense]
        ax2.semilogx(N_dense, norm, color=color, lw=1.8, marker=".", ms=4,
                     label=fr"$\theta_0={theta}$")

    ax1.set_xlabel(r"$N$", fontsize=12)
    ax1.set_ylabel(r"$\delta\cdot\log\log N$", fontsize=12)
    ax1.set_title(r"Absolute $b_{\rm eff}(N)$ — dashed = $f(\theta_0)$ heuristic", fontsize=10)
    ax1.legend(fontsize=7, ncol=2, loc="upper left")

    ax2.axhline(1.0, color="black", lw=1.2, ls=":", label="ratio = 1")
    ax2.set_xlabel(r"$N$", fontsize=12)
    ax2.set_ylabel(r"$b_{\rm eff}(N)\,/\,f(\theta_0)$", fontsize=12)
    ax2.set_title(r"Normalised: ratio to heuristic $f(\theta_0)$", fontsize=10)
    ax2.legend(fontsize=7, ncol=2, loc="lower right")

    plt.tight_layout()
    out = "stabilisation_test.png"
    plt.savefig(out, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close()
    print(f"\n  Saved: {out}")

    # Verdict
    print()
    print("  Flatness verdict (CV of beff over N ∈ [1e4, 1e8]):")
    print(f"  {'θ₀':>5}  {'mean_beff':>10}  {'std_beff':>9}  {'CV':>7}  verdict")
    print("  " + "-" * 50)
    for theta in THETAS_ALL:
        bvals_paper = []
        for n in N_VALS:
            P, _, _ = count_P(int(n), theta, primes)
            delta = 1.0 - P
            bvals_paper.append(delta * math.log(math.log(n)))
        mean_b = np.mean(bvals_paper)
        std_b  = np.std(bvals_paper)
        cv     = std_b / mean_b if mean_b > 0 else float("nan")
        verdict = "flat (CV<0.05)" if cv < 0.05 else "drifting"
        print(f"  {theta:>5.2f}  {mean_b:>10.4f}  {std_b:>9.5f}  {cv:>7.4f}  {verdict}")

# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    print("=" * 72)
    print("verify_paper_tables.py — All six tasks")
    print("=" * 72)

    # Load primes — 50M cache covers all queries
    # max q needed: 1e8/2 = 5e7; max p needed: sqrt(1e8) = 1e4
    prime_limit = 50_000_000
    t0 = time.time()
    primes = load_primes(prime_limit)
    print(f"\nLoaded {len(primes):,} primes to {prime_limit:,} ({time.time()-t0:.2f}s)")
    print(f"log log N range over [1e4, 1e8]: "
          f"{math.log(math.log(1e4)):.4f} to {math.log(math.log(1e8)):.4f}")

    t1 = time.time()
    table1_results = task1_verify_table1(primes)
    print(f"  (elapsed: {time.time()-t1:.1f}s)")

    t1 = time.time()
    task2_table2(table1_results)
    print(f"  (elapsed: {time.time()-t1:.1f}s)")

    t1 = time.time()
    task3_table3(primes)
    print(f"  (elapsed: {time.time()-t1:.1f}s)")

    t1 = time.time()
    task4_activation_plot(primes)
    print(f"  (elapsed: {time.time()-t1:.1f}s)")

    t1 = time.time()
    task5_extended_table2(primes)
    print(f"  (elapsed: {time.time()-t1:.1f}s)")

    t1 = time.time()
    task6_stabilisation(primes)
    print(f"  (elapsed: {time.time()-t1:.1f}s)")

    print()
    print("=" * 72)
    print("All tasks complete.")
    print("=" * 72)


if __name__ == "__main__":
    main()
