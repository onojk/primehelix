#!/usr/bin/env python3
"""
literature_check.py

Compares the exact sum S(N, θ₀) = Σ_{p ≤ N^θ₀} π(N/p)  against the
Montgomery-Vaughan PNT approximation MV(N, θ₀) = Σ_{p ≤ N^θ₀} N/(p·log(N/p)).

Also checks:
  - MV2: second-order PNT term  N/(p·log²(N/p))  added as correction
  - e^γ Mertens factor (γ = Euler-Mascheroni ≈ 0.5772)
  - 1/log N secondary-term fit to the ratio
  - Whether ratio → 1 (finite-N effect) or stabilises (structural bias)

Saves table to literature_check_results.csv.
"""

import csv
import math
import pickle
import pathlib
import time
from bisect import bisect_left, bisect_right
from math import isqrt, ceil, log, exp

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Constants ─────────────────────────────────────────────────────────────────

GAMMA   = 0.5772156649015328   # Euler-Mascheroni
EXP_G   = exp(GAMMA)           # ≈ 1.7813

THETAS  = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45]
N_VALS  = [1e4, 1e5, 1e6, 1e7, 1e8]
CACHE   = pathlib.Path(".cache")

# ── Prime loading ─────────────────────────────────────────────────────────────

def load_primes(limit: int):
    CACHE.mkdir(exist_ok=True)
    path = CACHE / f"primes_{limit}.pkl"
    if path.exists():
        with open(path, "rb") as f:
            return pickle.load(f)
    print(f"Sieving to {limit:,} …")
    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, isqrt(limit) + 1):
        if is_prime[i]:
            is_prime[i*i:limit+1:i] = b"\x00" * len(is_prime[i*i:limit+1:i])
    primes = [i for i in range(2, limit + 1) if is_prime[i]]
    with open(path, "wb") as f:
        pickle.dump(primes, f)
    return primes

# ── Core sums ─────────────────────────────────────────────────────────────────

def compute_S_and_MV(limit: int, theta: float, primes: list):
    """
    S_exact  = Σ_{p ≤ limit^θ} π(limit/p)            (exact, q unrestricted)
    S_proper = Σ_{p ≤ limit^θ} (π(limit/p) - π(p-1)) (q > p, unordered pairs)
    MV       = Σ_{p ≤ limit^θ} limit/(p·log(limit/p)) (PNT leading term)
    MV2      = MV + Σ limit/(p·log²(limit/p))          (PNT with 2nd order term)
    """
    p_max  = int(limit ** theta)
    S_raw  = 0
    S_prop = 0
    MV     = 0.0
    MV2    = 0.0

    for p in primes:
        if p > p_max:
            break
        qmax = limit // p
        if qmax < 2:
            continue

        pi_qmax = bisect_right(primes, qmax)         # π(⌊limit/p⌋)
        pi_pm1  = bisect_left(primes, p)             # π(p-1) = #{primes < p}

        S_raw  += pi_qmax
        S_prop += pi_qmax - pi_pm1                   # unordered q > p

        x      = limit / p                           # continuous N/p
        logx   = log(x)
        MV    += x / logx
        MV2   += x / logx + x / (logx * logx)       # leading + 2nd order

    return S_raw, S_prop, MV, MV2


def compute_T_total(limit: int, primes: list):
    """Total unordered semiprimes ≤ limit."""
    sqrt_limit = isqrt(limit)
    T = 0
    for p in primes:
        if p > sqrt_limit:
            break
        qmax   = limit // p
        pi_pm1 = bisect_left(primes, p)
        T     += bisect_right(primes, qmax) - pi_pm1
    return T

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    limit_max  = int(max(N_VALS))
    prime_limit = limit_max // 2

    print("=" * 78)
    print("literature_check.py")
    print("S(N,θ)=Σπ(N/p)  vs  MV=ΣN/(p·logN/p)  — ratio and correction checks")
    print("=" * 78)
    print()

    t0 = time.time()
    primes = load_primes(prime_limit)
    print(f"Primes to {prime_limit:,}: {len(primes):,}  ({time.time()-t0:.2f}s)")
    print(f"γ = {GAMMA:.6f},  e^γ = {EXP_G:.6f}")
    print()

    rows = []

    for theta in THETAS:
        theory_pref = log(1.0 / theta)

        print(f"θ₀ = {theta}")
        print(f"  {'N':>8}  {'S_prop':>10}  {'MV':>12}  {'MV2':>12}  "
              f"{'ratio':>7}  {'ratio·e^γ':>10}  {'MV2/S':>8}  {'1/logN':>8}")
        print("  " + "-" * 78)

        ratios   = []
        log_N_inv = []

        for n in N_VALS:
            limit    = int(n)
            S_raw, S_prop, MV, MV2 = compute_S_and_MV(limit, theta, primes)

            ratio      = S_prop / MV    if MV   > 0 else float("nan")
            ratio_MV2  = S_prop / MV2   if MV2  > 0 else float("nan")
            ratio_egam = ratio * EXP_G
            inv_logN   = 1.0 / log(n)

            print(f"  {n:>8.0e}  {S_prop:>10,}  {MV:>12.1f}  {MV2:>12.1f}  "
                  f"{ratio:>7.5f}  {ratio_egam:>10.5f}  {ratio_MV2:>8.5f}  "
                  f"{inv_logN:>8.5f}")

            ratios.append(ratio)
            log_N_inv.append(inv_logN)
            rows.append({
                "theta0":       theta,
                "N":            limit,
                "S_prop":       S_prop,
                "MV":           round(MV, 3),
                "MV2":          round(MV2, 3),
                "ratio_S_MV":   round(ratio, 6),
                "ratio_MV2_S":  round(ratio_MV2, 6),
                "ratio_x_egamma": round(ratio_egam, 6),
                "inv_logN":     round(inv_logN, 6),
            })

        # Fit ratio = a + b/log N to check secondary-term structure
        ratios    = np.array(ratios)
        log_N_inv = np.array(log_N_inv)
        b, a      = np.polyfit(log_N_inv, ratios, 1)   # ratio = a + b*(1/logN)
        fit       = a + b * log_N_inv
        ss_res    = np.sum((ratios - fit) ** 2)
        ss_tot    = np.sum((ratios - ratios.mean()) ** 2)
        r2_fit    = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0

        ratio_inf = a    # extrapolated limit as N→∞
        drift     = ratios[-1] - ratios[0]

        print()
        print(f"  Fit  ratio = {a:.5f} + {b:.5f}/log N  (R²={r2_fit:.4f})")
        print(f"  Extrapolated ratio as N→∞:  {ratio_inf:.5f}")
        print(f"  Drift [N=1e4 → N=1e8]:  {drift:+.5f}")

        # Verdict
        if abs(ratio_inf - 1.0) < 0.015 and r2_fit > 0.90:
            verdict = "FINITE-N EFFECT — ratio → 1 with 1/log N correction"
        elif abs(ratio_inf - 1.0) < 0.015:
            verdict = "FINITE-N EFFECT — ratio consistent with → 1"
        elif r2_fit > 0.90:
            verdict = f"STRUCTURAL BIAS — ratio → {ratio_inf:.4f} ≠ 1 (well-fit by 1/log N)"
        else:
            verdict = f"STRUCTURAL BIAS — ratio → {ratio_inf:.4f} ≠ 1 (complex structure)"
        print(f"  Verdict: {verdict}")
        print()

    # ── Summary table ─────────────────────────────────────────────────────────
    print()
    print("SUMMARY: ratio S_proper / MV at N=1e8")
    print("=" * 78)
    print(f"  {'θ₀':>6}  {'ratio@1e8':>10}  {'ratio·e^γ@1e8':>14}  "
          f"{'ratio→∞':>10}  {'verdict'}")
    print("  " + "-" * 70)

    for theta in THETAS:
        theta_rows = [r for r in rows if r["theta0"] == theta]
        r_1e8  = theta_rows[-1]["ratio_S_MV"]
        rg_1e8 = theta_rows[-1]["ratio_x_egamma"]
        inv_lN = np.array([r["inv_logN"] for r in theta_rows])
        ratios = np.array([r["ratio_S_MV"] for r in theta_rows])
        b, a   = np.polyfit(inv_lN, ratios, 1)

        drift  = ratios[-1] - ratios[0]
        approaching_1 = abs(a - 1.0) < 0.02
        tag = "→ 1" if approaching_1 else f"→ {a:.4f}"

        print(f"  {theta:>6.2f}  {r_1e8:>10.5f}  {rg_1e8:>14.5f}  "
              f"{a:>10.5f}  {tag}")

    print()
    print(f"e^γ = {EXP_G:.5f}   (Mertens factor — does ratio·e^γ ≈ 1 anywhere? check above)")
    print()

    # ── Overall interpretation ────────────────────────────────────────────────
    ratios_inf = []
    for theta in THETAS:
        theta_rows = [r for r in rows if r["theta0"] == theta]
        inv_lN = np.array([r["inv_logN"] for r in theta_rows])
        ratios = np.array([r["ratio_S_MV"] for r in theta_rows])
        b, a   = np.polyfit(inv_lN, ratios, 1)
        ratios_inf.append(a)

    all_near_1 = all(abs(a - 1.0) < 0.03 for a in ratios_inf)
    print("OVERALL INTERPRETATION:")
    if all_near_1:
        print("  All extrapolated ratios → 1.  The gap between S_proper and MV is")
        print("  a pure finite-N effect, fully explained by the 1/log N secondary")
        print("  PNT term.  No structural bias in the approximation.")
        print()
        print("  CONSEQUENCE FOR PAPER: the prefactor deficit in b_eff comes from")
        print("  finite-N effects in the PNT, not from missing structure.  The")
        print("  paper statement 'consistent with O(1/log log N) with subleading")
        print("  corrections' is correct; the subleading term is ~1/log N.")
    else:
        mean_inf = np.mean(ratios_inf)
        print(f"  Extrapolated ratios do NOT converge to 1 (mean = {mean_inf:.4f}).")
        print("  There is a systematic bias in the MV approximation beyond 1/log N.")
        print()
        print("  CONSEQUENCE FOR PAPER: the prefactor gap is not purely a finite-N")
        print("  effect.  A constant factor persists in the asymptotic approximation.")
        print("  This may reflect Mertens-type corrections not explicit in MV §7.4.")

    # ── Plot ──────────────────────────────────────────────────────────────────
    COLORS = plt.cm.tab10(np.linspace(0, 0.9, len(THETAS)))
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    fig.suptitle(
        r"$S_{\mathrm{proper}}(N,\theta_0)$ / MV$(N,\theta_0)$:  "
        r"ratio and correction structure",
        fontsize=12
    )

    ax1, ax2 = axes
    for ax in axes:
        ax.set_facecolor("white")
        ax.grid(color="#dddddd", lw=0.5)

    # Panel 1: ratio vs 1/log N
    for theta, color in zip(THETAS, COLORS):
        theta_rows = [r for r in rows if r["theta0"] == theta]
        x    = [r["inv_logN"] for r in theta_rows]
        y    = [r["ratio_S_MV"] for r in theta_rows]
        b_c, a_c = np.polyfit(x, y, 1)
        x_s  = np.linspace(0, max(x) * 1.05, 100)
        ax1.scatter(x, y, color=color, s=30, zorder=5)
        ax1.plot(x_s, a_c + b_c * x_s, color=color, lw=1.2, ls="--",
                 label=fr"$\theta_0={theta}$  ($r_\infty={a_c:.3f}$)")

    ax1.axhline(1.0, color="black", lw=1.2, ls=":", label=r"ratio $= 1$")
    ax1.set_xlabel(r"$1/\log N$", fontsize=12)
    ax1.set_ylabel(r"$S_{\mathrm{proper}} \,/\, \mathrm{MV}$", fontsize=12)
    ax1.set_title(r"Ratio vs $1/\log N$  (extrapolate to $N\to\infty$ at $x=0$)", fontsize=10)
    ax1.legend(fontsize=7, ncol=2, loc="upper left")

    # Panel 2: ratio vs log(log(N))
    llN_vals = [log(log(n)) for n in N_VALS]
    for theta, color in zip(THETAS, COLORS):
        theta_rows = [r for r in rows if r["theta0"] == theta]
        y = [r["ratio_S_MV"] for r in theta_rows]
        ax2.plot(llN_vals, y, color=color, marker="o", ms=4, lw=1.5,
                 label=fr"$\theta_0={theta}$")

    ax2.axhline(1.0, color="black", lw=1.2, ls=":", label=r"ratio $= 1$")
    ax2.set_xlabel(r"$\log(\log(N))$", fontsize=12)
    ax2.set_ylabel(r"$S_{\mathrm{proper}} \,/\, \mathrm{MV}$", fontsize=12)
    ax2.set_title(r"Ratio vs $\log(\log(N))$", fontsize=10)

    # Top x-axis: N values
    ax2t = ax2.twiny()
    ax2t.set_xlim(ax2.get_xlim())
    ax2t.set_xticks(llN_vals)
    ax2t.set_xticklabels([r"$10^{%d}$" % int(math.log10(n)) for n in N_VALS], fontsize=8)
    ax2t.set_xlabel(r"$N$", fontsize=10)

    plt.tight_layout()
    plt.savefig("literature_check_plot.png", dpi=155, facecolor="white", bbox_inches="tight")
    plt.close()
    print()
    print("Saved: literature_check_plot.png")

    # ── CSV ───────────────────────────────────────────────────────────────────
    fname = "literature_check_results.csv"
    with open(fname, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"Saved: {fname}")


if __name__ == "__main__":
    main()
