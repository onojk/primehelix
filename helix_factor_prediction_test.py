#!/usr/bin/env python3
"""
helix_factor_prediction_test.py

HYPOTHESIS: Semiprimes with the same smaller factor p cluster at predictable
angular positions on the PrimeHelix.  If true, helix angle carries factor
information — position has predictive power over factorization.

Three experiments test this from different angles (pun intended).
"""

import array
import json
import math
import random
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Constants ─────────────────────────────────────────────────────────────────
SCAN_LIMIT  = 1_000_000
ANGLE_STEP  = 0.27          # radians per lopsided-index step (from lattice.py)
TWO_PI      = 2 * math.pi
N_SECTORS   = 120           # angular bins for Exp 3 sector table
N_PERMS     = 500           # permutation iterations for Exp 1
MIN_GROUP   = 30            # min semiprimes per p value for Exp 1 stats
TEST_N      = 1000          # test-set size for Exp 3

# ── Colour palette ────────────────────────────────────────────────────────────
DARK_BG    = "#0a0a12"
PANEL_BG   = "#12121e"
SPINE_COL  = "#2a2a40"
TEXT_COL   = "#ccccee"
PALETTE    = ["#9f8ff8", "#3ecf96", "#f07050", "#ffd080", "#aaaaaa",
              "#ff7799", "#44ccff", "#ffaa44", "#88ff88", "#ff88ff"]

# ═══════════════════════════════════════════════════════════════════════════════
# SIEVE  (mirrors primehelix_lattice.py exactly)
# ═══════════════════════════════════════════════════════════════════════════════

def _lopsided(p: int, q: int) -> bool:
    """balance = (q-p)/sqrt(pq) >= 10  ↔  (q-p)^2 >= 100·p·q"""
    d = q - p
    return d * d >= 100 * p * q


def build_sieve(limit: int = SCAN_LIMIT):
    """Return sorted list of (n, p, q) for every lopsided semiprime n≤limit."""
    print(f"Sieving to {limit:,} …", end=" ", flush=True)
    t0 = time.time()
    spf = array.array("i", range(limit + 1))
    i = 2
    while i * i <= limit:
        if spf[i] == i:
            for j in range(i * i, limit + 1, i):
                if spf[j] == j:
                    spf[j] = i
        i += 1

    sps = []
    for n in range(4, limit + 1):
        p = spf[n]
        if p == n:
            continue
        q = n // p
        if spf[q] != q:
            continue
        if _lopsided(p, q):
            sps.append((n, p, q))

    print(f"{len(sps):,} lopsided semiprimes  ({time.time()-t0:.2f}s)")
    return sps   # lop_idx = list index

# ═══════════════════════════════════════════════════════════════════════════════
# CIRCULAR STATISTICS
# ═══════════════════════════════════════════════════════════════════════════════

def circ_R(angles: np.ndarray) -> float:
    """Mean resultant length R ∈ [0,1].  R≈0: uniform; R≈1: concentrated."""
    if len(angles) == 0:
        return 0.0
    C = np.cos(angles).mean()
    S = np.sin(angles).mean()
    return float(np.sqrt(C * C + S * S))


def rayleigh_pval(n: int, R: float) -> float:
    """Approximate p-value for Rayleigh test (Mardia & Jupp 2000)."""
    Z = n * R * R
    if Z > 60:
        return 0.0
    pv = math.exp(-Z) * (1 + (2 * Z - Z * Z) / (4 * n))
    return max(0.0, min(1.0, pv))

# ═══════════════════════════════════════════════════════════════════════════════
# REGRESSION HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

def linreg(x: np.ndarray, y: np.ndarray):
    """OLS: y = a + b*x.  Returns (a, b, R²)."""
    xm, ym = x.mean(), y.mean()
    dx = x - xm
    cov = (dx * (y - ym)).sum()
    var = (dx * dx).sum()
    if var < 1e-30:
        return float(ym), 0.0, 0.0
    b = float(cov / var)
    a = float(ym - b * xm)
    pred = a + b * x
    ss_res = float(((y - pred) ** 2).sum())
    ss_tot = float(((y - ym) ** 2).sum())
    r2 = 1.0 - ss_res / max(ss_tot, 1e-30)
    return a, b, r2

# ═══════════════════════════════════════════════════════════════════════════════
# EXPERIMENT 1 — Angular clustering by small factor p
# ═══════════════════════════════════════════════════════════════════════════════

def experiment_1(lop_sps: list) -> dict:
    print("\n" + "═" * 70)
    print("EXPERIMENT 1  Angular clustering by small factor p")
    print("  Null hypothesis: helix angles for a given p are uniformly distributed.")
    print("  Test: compare mean R across p-groups to 500-iteration permutation baseline.")
    print("═" * 70)

    n_lop = len(lop_sps)
    # Assign each lopsided semiprime its helix angle by sequential index
    all_angles = np.array([(i * ANGLE_STEP) % TWO_PI for i in range(n_lop)])
    all_ps     = np.array([p for _, p, _ in lop_sps], dtype=np.int32)

    # ── Group by p ────────────────────────────────────────────────────────────
    p_groups: dict[int, np.ndarray] = {}
    for pv in np.unique(all_ps):
        mask = all_ps == pv
        if mask.sum() >= MIN_GROUP:
            p_groups[int(pv)] = all_angles[mask]

    print(f"\n  Lopsided semiprimes total : {n_lop:,}")
    print(f"  Distinct p values (n≥{MIN_GROUP}) : {len(p_groups)}")

    # ── Observed R per p ──────────────────────────────────────────────────────
    obs_R = {p: circ_R(a) for p, a in p_groups.items()}
    mean_obs = float(np.mean(list(obs_R.values())))

    # ── Permutation baseline ──────────────────────────────────────────────────
    print(f"  Running {N_PERMS}-iteration permutation test …", flush=True)
    counts = np.array([len(p_groups[p]) for p in sorted(p_groups)])
    buf    = all_angles.copy()
    rng    = np.random.default_rng(42)
    perm_means = np.empty(N_PERMS)
    for t in range(N_PERMS):
        rng.shuffle(buf)
        rs, pos = [], 0
        for c in counts:
            rs.append(circ_R(buf[pos: pos + c]))
            pos += c
        perm_means[t] = np.mean(rs)

    base_mean = float(perm_means.mean())
    base_std  = float(perm_means.std())
    z_score   = (mean_obs - base_mean) / max(base_std, 1e-12)

    print(f"\n  Permutation baseline mean R : {base_mean:.6f} ± {base_std:.6f}")
    print(f"  Observed  mean R            : {mean_obs:.6f}")
    print(f"  Z-score vs baseline         : {z_score:+.2f}  "
          f"{'← SIGNIFICANT (|z|>3)' if abs(z_score) > 3 else '(not significant)'}")

    # ── Per-p table ───────────────────────────────────────────────────────────
    sorted_p = sorted(obs_R.items(), key=lambda x: -x[1])
    print(f"\n  Top-20 most-clustered p values (highest R):")
    hdr = f"  {'p':>8}  {'count':>8}  {'R':>9}  {'circ_std':>10}  {'Rayleigh-p':>12}  sig"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    top_rows = []
    for pv, Rv in sorted_p[:20]:
        angs  = p_groups[pv]
        cs    = float(np.sqrt(max(-2 * math.log(max(Rv, 1e-15)), 0)))
        pval  = rayleigh_pval(len(angs), Rv)
        sig   = " *" if pval < 0.05 else ""
        print(f"  {pv:>8}  {len(angs):>8,}  {Rv:>9.6f}  {cs:>10.4f}  "
              f"{pval:>12.6f}{sig}")
        top_rows.append({"p": pv, "n": len(angs), "R": Rv,
                         "circ_std": cs, "rayleigh_p": pval})

    n_sig = sum(1 for pv, Rv in obs_R.items()
                if rayleigh_pval(len(p_groups[pv]), Rv) < 0.05)
    print(f"\n  Rayleigh-significant (p<0.05): {n_sig}/{len(p_groups)}")

    # ── Interpretation ────────────────────────────────────────────────────────
    if abs(z_score) <= 3 and n_sig < len(p_groups) * 0.10:
        interp = ("NO clustering detected.  p-values are uniformly interleaved on "
                  "the helix — angle carries NO information about small factor p.")
    elif z_score > 3:
        interp = (f"CLUSTERING DETECTED (z={z_score:.1f}).  "
                  f"{n_sig}/{len(p_groups)} p-groups are non-uniformly distributed.")
    else:
        interp = ("Marginal result — weak or inconsistent clustering signal.")
    print(f"\n  Interpretation: {interp}")

    return {
        "n_lopsided": n_lop,
        "n_qualifying_p": len(p_groups),
        "baseline_mean_R": base_mean,
        "baseline_std_R":  base_std,
        "observed_mean_R": mean_obs,
        "z_score": z_score,
        "significant": bool(abs(z_score) > 3),
        "n_rayleigh_sig": n_sig,
        "top_p": top_rows,
        "interpretation": interp,
    }, p_groups, all_angles, all_ps


# ═══════════════════════════════════════════════════════════════════════════════
# EXPERIMENT 2 — Factor ratio prediction from θ
# ═══════════════════════════════════════════════════════════════════════════════

def experiment_2(lop_sps: list) -> dict:
    print("\n" + "═" * 70)
    print("EXPERIMENT 2  Factor ratio prediction from θ = log(p)/log(n)")
    print("  Tests whether θ or helix angle predict the factor imbalance p/q.")
    print("═" * 70)

    ns = np.array([n for n, p, q in lop_sps], dtype=np.float64)
    ps = np.array([p for n, p, q in lop_sps], dtype=np.float64)
    qs = np.array([q for n, p, q in lop_sps], dtype=np.float64)
    N  = len(lop_sps)

    theta      = np.log(ps) / np.log(ns)     # log(p)/log(n) ∈ (0, 0.5)
    ratio      = ps / qs                      # p/q ∈ (0, 1); small = very lopsided
    log_ratio  = np.log(ratio)
    log_theta  = np.log(theta)
    helix_ang  = (np.arange(N) * ANGLE_STEP) % TWO_PI

    # A. log-log regression: log(ratio) ~ a + b·log(θ)
    a_ll, b_ll, r2_ll = linreg(log_theta, log_ratio)

    # B. Mathematical identity check: log(p/q) = (2θ-1)·log(n)  [EXACT]
    #    Proof: log(p)=θ·log n, log(q)=(1-θ)·log n, so log(p/q)=(2θ-1)·log n.
    taut = (2 * theta - 1) * np.log(ns)
    ss_res_t  = float(((log_ratio - taut) ** 2).sum())
    ss_tot_t  = float(((log_ratio - log_ratio.mean()) ** 2).sum())
    r2_taut   = 1.0 - ss_res_t / max(ss_tot_t, 1e-30)

    # C. Helix angle vs log(ratio) — the novel question
    a_ha, b_ha, r2_ha = linreg(helix_ang, log_ratio)

    # D. Helix angle vs θ — does angle proxy for the structural parameter?
    a_ht, b_ht, r2_ht = linreg(helix_ang, theta)

    # E. Multiple regression: log(ratio) ~ β₀ + β₁·angle + β₂·log(n)
    X    = np.column_stack([np.ones(N), helix_ang, np.log(ns)])
    beta = np.linalg.lstsq(X, log_ratio, rcond=None)[0]
    r2_multi = float(1 - ((log_ratio - X @ beta) ** 2).sum() /
                     max(((log_ratio - log_ratio.mean()) ** 2).sum(), 1e-30))

    resid_std = float((log_ratio - (a_ll + b_ll * log_theta)).std())

    print(f"\n  Samples: {N:,}")
    print(f"  θ range    : [{theta.min():.5f}, {theta.max():.5f}]")
    print(f"  ratio range: [{ratio.min():.6f}, {ratio.max():.5f}]")
    print()
    print(f"  A. log(ratio) = {a_ll:.4f} + {b_ll:.4f}·log(θ)   R² = {r2_ll:.5f}")
    print(f"     [θ alone, no log n — expect low R² since n varies 1–1M]")
    print()
    print(f"  B. Tautology: log(ratio) = (2θ−1)·log(n)          R² = {r2_taut:.6f}")
    print(f"     [Must be ≈1.000 — confirms p·q=n algebraic identity]")
    print()
    print(f"  C. Helix angle vs log(ratio)                       R² = {r2_ha:.5f}")
    ang_interp = ("← meaningful signal" if r2_ha > 0.05
                  else "← weak: angle is a noisy proxy for log n, not for θ")
    print(f"     {ang_interp}")
    print()
    print(f"  D. Helix angle vs θ                                R² = {r2_ht:.5f}")
    print()
    print(f"  E. Multiple reg: log(ratio) ~ 1 + angle + log(n)  R² = {r2_multi:.5f}")
    print(f"     β = [{beta[0]:.4f},  {beta[1]:.6f}·angle,  {beta[2]:.4f}·log n]")
    print()
    print(f"  Residual std (log-log fit): {resid_std:.4f} nats")
    print()
    print("  KEY INSIGHT: The identity in (B) is exact — θ and log(n) together")
    print("  determine ratio perfectly.  The question is whether the HELIX ANGLE")
    print("  (= lop_idx * 0.27 mod 2π) adds independent information.  R²(C)")
    print("  measures this: it captures only the angle's power, excluding log n.")

    return {
        "n_samples": N,
        "loglog_a": a_ll, "loglog_b": b_ll, "loglog_r2": r2_ll,
        "tautology_r2": float(r2_taut),
        "helix_angle_vs_ratio_r2": float(r2_ha),
        "helix_angle_vs_theta_r2": float(r2_ht),
        "multi_reg_r2": float(r2_multi),
        "multi_reg_beta": beta.tolist(),
        "residual_std": resid_std,
        # Sample points for plotting (all, downsampled for large arrays)
        "_plot_theta":   theta[::max(1, N//4000)].tolist(),
        "_plot_ratio":   ratio[::max(1, N//4000)].tolist(),
        "_plot_helix":   helix_ang[::max(1, N//4000)].tolist(),
        "_plot_lograt":  log_ratio[::max(1, N//4000)].tolist(),
    }, theta, ratio, helix_ang, log_ratio, log_theta


# ═══════════════════════════════════════════════════════════════════════════════
# EXPERIMENT 3 — Search space reduction
# ═══════════════════════════════════════════════════════════════════════════════

def experiment_3(lop_sps: list) -> dict:
    print("\n" + "═" * 70)
    print("EXPERIMENT 3  Search space reduction via helix angle")
    print("  Can helix sector narrow the list of candidate factors p?")
    print("  Baseline: top-K by global p-frequency.  Helix: top-K in local sector.")
    print("═" * 70)

    n_lop = len(lop_sps)
    ns_arr = np.array([n for n, p, q in lop_sps], dtype=np.float64)  # sorted

    # ── Helix angle for each lopsided semiprime ────────────────────────────────
    all_angles = (np.arange(n_lop) * ANGLE_STEP) % TWO_PI
    all_ps     = np.array([p for _, p, _ in lop_sps], dtype=np.int32)

    # ── Build sector → p distribution ─────────────────────────────────────────
    sector_w = TWO_PI / N_SECTORS
    sec_of   = (all_angles / sector_w).astype(int) % N_SECTORS
    sector_p: list[Counter] = [Counter() for _ in range(N_SECTORS)]
    for idx in range(n_lop):
        sector_p[sec_of[idx]][int(all_ps[idx])] += 1

    global_counts = Counter(all_ps.tolist())

    # ── Angle estimator from n (no factoring required) ─────────────────────────
    # Given n, count of lopsided semiprimes ≤ n = searchsorted in sorted ns_arr.
    # This is knowledge-free: we look up n in the pre-sieved table, not factor it.
    def angle_from_n(n_val: float) -> float:
        idx = int(np.searchsorted(ns_arr, n_val, side="left"))
        return (idx * ANGLE_STEP) % TWO_PI

    # Smooth approximation: no sieve table, only n's magnitude
    # lop_count(n) ≈ 0.73 · n · ln(ln(n)) / ln(n)   (asymptotic)
    def angle_from_magnitude(n_val: float) -> float:
        ln_n = math.log(n_val)
        lln  = math.log(max(ln_n, 1.1))
        est  = 0.73 * n_val * lln / ln_n
        return (est * ANGLE_STEP) % TWO_PI

    # ── Test set ───────────────────────────────────────────────────────────────
    rng       = random.Random(42)
    test_idxs = rng.sample(range(n_lop), min(TEST_N, n_lop))
    test      = [(lop_sps[i], i) for i in test_idxs]

    Ks = [1, 3, 5, 10, 20, 50]
    # Full trial-division search for n ≤ 1M: primes ≤ sqrt(1M) ≈ 1000 → 168 primes
    FULL = 168

    rows = []
    hit_A = {}; hit_B = {}; hit_glob = {}

    for K in Ks:
        g_top = {p for p, _ in global_counts.most_common(K)}
        cA = cB = cG = 0

        for (n_val, p_true, q_true), true_idx in test:
            # Strategy A: exact angle from sieve index (upper-bound, no factoring)
            ang_A = (true_idx * ANGLE_STEP) % TWO_PI
            sec_A = int(ang_A / sector_w) % N_SECTORS
            top_A = {p for p, _ in sector_p[sec_A].most_common(K)}

            # Strategy B: smooth magnitude approximation (truly knowledge-free)
            ang_B = angle_from_magnitude(float(n_val))
            sec_B = int(ang_B / sector_w) % N_SECTORS
            top_B = {p for p, _ in sector_p[sec_B].most_common(K)}

            if p_true in top_A: cA += 1
            if p_true in top_B: cB += 1
            if p_true in g_top: cG += 1

        n_t = len(test)
        hit_A[K]    = cA / n_t
        hit_B[K]    = cB / n_t
        hit_glob[K] = cG / n_t
        rows.append({
            "K": K,
            "helix_exact_hr": cA / n_t,
            "helix_approx_hr": cB / n_t,
            "global_hr": cG / n_t,
            "search_reduction_pct": round((1 - K / FULL) * 100, 1),
        })

    # ── Print results ──────────────────────────────────────────────────────────
    print(f"\n  Test set      : {len(test)} random lopsided semiprimes")
    print(f"  Sectors       : {N_SECTORS}  (width {math.degrees(sector_w):.1f}°)")
    print(f"  Full search   : {FULL} primes (trial division ≤ √1M)")
    print()
    print(f"  {'K':>4}  {'Helix-exact':>12}  {'Helix-approx':>13}  "
          f"{'Global':>8}  {'Reduction':>11}")
    print(f"  {'-'*57}")
    for row in rows:
        K = row["K"]
        print(f"  {K:>4}  {row['helix_exact_hr']:>11.1%}  "
              f"{row['helix_approx_hr']:>13.1%}  "
              f"{row['global_hr']:>8.1%}  "
              f"{row['search_reduction_pct']:>9.0f}%")

    # ── Lift analysis ──────────────────────────────────────────────────────────
    lifts_exact  = [hit_A[K] - hit_glob[K] for K in Ks]
    lifts_approx = [hit_B[K] - hit_glob[K] for K in Ks]
    mean_lift_exact  = float(np.mean(lifts_exact))
    mean_lift_approx = float(np.mean(lifts_approx))

    print(f"\n  Mean helix-exact  lift vs global : {mean_lift_exact:+.4f} "
          f"({mean_lift_exact*100:+.2f} pp)")
    print(f"  Mean helix-approx lift vs global : {mean_lift_approx:+.4f} "
          f"({mean_lift_approx*100:+.2f} pp)")

    # ── Hard cases: p > 50 (trial division finds nothing) ─────────────────────
    hard = [(sp, idx) for sp, idx in test if sp[1] > 50]
    hard_stats = {}
    if hard:
        print(f"\n  Hard semiprimes (p>50): {len(hard)}")
        for K in [10, 20]:
            g_top = {p for p, _ in global_counts.most_common(K)}
            cA = cG = 0
            for (n_val, p_true, _), true_idx in hard:
                ang_A = (true_idx * ANGLE_STEP) % TWO_PI
                sec_A = int(ang_A / sector_w) % N_SECTORS
                top_A = {p for p, _ in sector_p[sec_A].most_common(K)}
                if p_true in top_A: cA += 1
                if p_true in g_top: cG += 1
            n_h = len(hard)
            hr_A = cA / n_h; hr_G = cG / n_h
            print(f"    K={K:>2}: helix {hr_A:.1%}  global {hr_G:.1%}  "
                  f"lift {hr_A-hr_G:+.3f}")
            hard_stats[K] = {"helix": hr_A, "global": hr_G, "lift": hr_A - hr_G}
    else:
        print("  No hard semiprimes (p>50) in test set.")

    # ── Verdict ────────────────────────────────────────────────────────────────
    thr_meaningful = 0.05
    thr_modest     = 0.01
    if mean_lift_exact > thr_meaningful:
        verdict = (f"HELIX PROVIDES MEANINGFUL LIFT ({mean_lift_exact*100:+.2f} pp) — "
                   "angular position carries real factor information")
    elif mean_lift_exact > thr_modest:
        verdict = (f"HELIX PROVIDES MODEST LIFT ({mean_lift_exact*100:+.2f} pp) — "
                   "angular position weakly informative")
    elif abs(mean_lift_exact) <= thr_modest:
        verdict = ("HELIX MATCHES GLOBAL BASELINE — "
                   "p-values are uniformly distributed across sectors; "
                   "angle adds NO predictive power beyond global frequency")
    else:
        verdict = ("GLOBAL OUTPERFORMS HELIX — "
                   "sector table is noisier than global frequency")

    print(f"\n  VERDICT: {verdict}")

    # ── Key metric print ───────────────────────────────────────────────────────
    k10 = next(r for r in rows if r["K"] == 10)
    print()
    print(f"  KEY METRIC (K=10, {FULL}-prime full search):")
    print(f"    Helix-guided : {k10['helix_exact_hr']:.1%} hit rate  "
          f"using {10}/{FULL} = {10/FULL:.1%} of search space")
    print(f"    Global freq  : {k10['global_hr']:.1%} hit rate  "
          f"using {10}/{FULL} = {10/FULL:.1%} of search space")
    lift10 = k10["helix_exact_hr"] - k10["global_hr"]
    print(f"    Helix lift   : {lift10:+.1%}")
    if abs(lift10) < 0.02:
        print(f"    → Geometry is VISUALIZATION ONLY for this metric.")
    else:
        sign = "adds" if lift10 > 0 else "loses"
        print(f"    → Helix {sign} {abs(lift10)*100:.1f}pp vs global guess.")

    return {
        "n_test": len(test),
        "n_sectors": N_SECTORS,
        "full_search_size": FULL,
        "by_K": rows,
        "mean_lift_exact_pp": mean_lift_exact,
        "mean_lift_approx_pp": mean_lift_approx,
        "hard_cases": hard_stats,
        "verdict": verdict,
    }, hit_A, hit_glob


# ═══════════════════════════════════════════════════════════════════════════════
# PLOTS
# ═══════════════════════════════════════════════════════════════════════════════

def _style(ax):
    ax.set_facecolor(PANEL_BG)
    ax.tick_params(colors=TEXT_COL, labelsize=8)
    for s in ax.spines.values():
        s.set_color(SPINE_COL)
    ax.xaxis.label.set_color(TEXT_COL)
    ax.yaxis.label.set_color(TEXT_COL)
    ax.title.set_color(TEXT_COL)


def plot_theta_ratio(r2: dict, out_dir: Path):
    th  = np.array(r2["_plot_theta"])
    ra  = np.array(r2["_plot_ratio"])
    ha  = np.array(r2["_plot_helix"])
    lr  = np.array(r2["_plot_lograt"])
    r2_ll = r2["loglog_r2"]
    r2_ha = r2["helix_angle_vs_ratio_r2"]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.patch.set_facecolor(DARK_BG)

    # Panel A: log-log  theta vs ratio
    ax = axes[0]
    _style(ax)
    lth = np.log10(th)
    lra = np.log10(ra)
    ax.scatter(lth, lra, s=1.2, alpha=0.2, color="#3ecf96", rasterized=True)
    coef = np.polyfit(lth, lra, 1)
    xs = np.linspace(lth.min(), lth.max(), 300)
    ax.plot(xs, np.polyval(coef, xs), color="#ffd080", lw=1.8,
            label=f"slope={coef[0]:.3f}\nR²={r2_ll:.4f}")
    ax.set_xlabel("log₁₀(θ)   [= log₁₀(log p / log n)]")
    ax.set_ylabel("log₁₀(p / q)")
    ax.set_title("Factor ratio vs θ  (log-log)", pad=8)
    ax.legend(facecolor=PANEL_BG, edgecolor=SPINE_COL, labelcolor=TEXT_COL,
              fontsize=9)

    # Panel B: helix angle vs log(ratio)
    ax = axes[1]
    _style(ax)
    ax.scatter(ha, lr, s=1.2, alpha=0.2, color="#9f8ff8", rasterized=True)
    coef_h = np.polyfit(ha, lr, 1)
    xs_h = np.linspace(0, TWO_PI, 300)
    ax.plot(xs_h, np.polyval(coef_h, xs_h), color="#ffd080", lw=1.8,
            label=f"R²={r2_ha:.5f}")
    ax.set_xlabel("Helix angle (rad)")
    ax.set_ylabel("log(p / q)")
    ax.set_title("Helix angle vs log(factor ratio)", pad=8)
    ax.set_xticks([0, math.pi / 2, math.pi, 3 * math.pi / 2, TWO_PI])
    ax.set_xticklabels(["0", "π/2", "π", "3π/2", "2π"])
    ax.legend(facecolor=PANEL_BG, edgecolor=SPINE_COL, labelcolor=TEXT_COL,
              fontsize=9)

    plt.tight_layout(pad=0.6)
    p = out_dir / "theta_ratio_fit.png"
    plt.savefig(p, dpi=130, facecolor=DARK_BG)
    plt.close()
    print(f"  Saved: {p}")


def plot_angle_distributions(p_groups: dict, out_dir: Path):
    top10 = sorted(p_groups.items(), key=lambda x: -len(x[1]))[:10]
    N_BINS = 36
    edges  = np.linspace(0, TWO_PI, N_BINS + 1)
    ctrs   = (edges[:-1] + edges[1:]) / 2
    w      = TWO_PI / N_BINS

    fig, axes = plt.subplots(2, 5, figsize=(18, 7))
    fig.patch.set_facecolor(DARK_BG)
    fig.suptitle(
        "Helix angle distribution for the 10 most frequent small primes p\n"
        "White dashed = uniform expectation.  Peaks = angular clustering.",
        color=TEXT_COL, fontsize=11,
    )

    for k, ((pv, angs), col) in enumerate(zip(top10, PALETTE)):
        ax = axes[k // 5][k % 5]
        _style(ax)
        counts, _  = np.histogram(angs, bins=edges)
        expected   = len(angs) / N_BINS
        R          = circ_R(angs)
        pval       = rayleigh_pval(len(angs), R)

        ax.bar(ctrs, counts, width=w * 0.85, color=col, alpha=0.75, align="center")
        ax.axhline(expected, color="white", lw=0.9, ls="--", alpha=0.55)
        ax.set_title(f"p={pv}  n={len(angs):,}\nR={R:.4f}  ray-p={pval:.3g}",
                     color=TEXT_COL, fontsize=7.5)
        ax.set_xlim(0, TWO_PI)
        ax.set_xticks([0, math.pi, TWO_PI])
        ax.set_xticklabels(["0", "π", "2π"], fontsize=7)
        ax.set_xlabel("angle (rad)", fontsize=7)

    plt.tight_layout(pad=0.4)
    p = out_dir / "angle_distributions.png"
    plt.savefig(p, dpi=130, facecolor=DARK_BG)
    plt.close()
    print(f"  Saved: {p}")


def plot_search_reduction(hit_A: dict, hit_glob: dict, out_dir: Path):
    Ks  = sorted(hit_A.keys())
    hrs = np.array([hit_A[k] for k in Ks])
    hrg = np.array([hit_glob[k] for k in Ks])
    FULL = 168

    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor(DARK_BG)
    _style(ax)

    ax.plot(Ks, hrs * 100, "o-", color="#3ecf96", lw=2, ms=6,
            label="Helix sector (exact angle)")
    ax.plot(Ks, hrg * 100, "s--", color="#aaaaaa", lw=1.6, ms=5,
            label="Global frequency baseline")
    ax.plot(Ks, np.array(Ks) / FULL * 100, ":", color="#f07050", lw=1.2,
            label=f"Random  K/{FULL} primes")

    ax.fill_between(Ks, hrg * 100, hrs * 100,
                    where=(hrs >= hrg), alpha=0.15, color="#3ecf96",
                    label="Helix lift (positive)")
    ax.fill_between(Ks, hrg * 100, hrs * 100,
                    where=(hrs < hrg),  alpha=0.15, color="#f07050",
                    label="Helix deficit (negative)")

    ax.set_xlabel("K  (number of candidate p values considered)")
    ax.set_ylabel("Hit rate  (%)")
    ax.set_title("Experiment 3: Does helix sector narrow the factor search?",
                 color=TEXT_COL, pad=8)
    ax.set_ylim(0, 105)
    ax.grid(color=SPINE_COL, lw=0.5, alpha=0.7)
    ax.legend(facecolor=PANEL_BG, edgecolor=SPINE_COL, labelcolor=TEXT_COL,
              fontsize=8)

    plt.tight_layout()
    p = out_dir / "search_reduction.png"
    plt.savefig(p, dpi=130, facecolor=DARK_BG)
    plt.close()
    print(f"  Saved: {p}")


# ═══════════════════════════════════════════════════════════════════════════════
# JSON serialisation helper
# ═══════════════════════════════════════════════════════════════════════════════

class _Enc(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.integer,)):  return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray):     return obj.tolist()
        if isinstance(obj, bool):           return bool(obj)
        return super().default(obj)


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    out_dir = Path(".")
    random.seed(42)
    np.random.seed(42)

    print("=" * 70)
    print("helix_factor_prediction_test.py")
    print("Hypothesis: helix angle carries predictive information about p")
    print("=" * 70)

    lop_sps = build_sieve(SCAN_LIMIT)

    r1, p_groups, all_angles, all_ps = experiment_1(lop_sps)
    r2, theta, ratio, helix_a, log_r, log_th = experiment_2(lop_sps)
    r3, hit_A, hit_glob = experiment_3(lop_sps)

    print("\n" + "─" * 70)
    print("Saving plots …")
    plot_theta_ratio(r2, out_dir)
    plot_angle_distributions(p_groups, out_dir)
    plot_search_reduction(hit_A, hit_glob, out_dir)

    # Scrub private plot keys before serialising
    r2_json = {k: v for k, v in r2.items() if not k.startswith("_")}
    blob = {"experiment_1": r1, "experiment_2": r2_json, "experiment_3": r3}
    jpath = out_dir / "clustering_results.json"
    with open(jpath, "w") as f:
        json.dump(blob, f, indent=2, cls=_Enc)
    print(f"  Saved: {jpath}")

    # ── Final summary ──────────────────────────────────────────────────────────
    print("\n" + "═" * 70)
    print("FINAL SUMMARY")
    print("═" * 70)

    z  = r1["z_score"]
    sig = "SIGNIFICANT" if abs(z) > 3 else "not significant"
    print(f"\n  Exp 1 — Clustering:  mean R={r1['observed_mean_R']:.6f}  "
          f"baseline={r1['baseline_mean_R']:.6f}  z={z:+.2f}  ({sig})")
    print(f"    {r1['n_rayleigh_sig']}/{r1['n_qualifying_p']} p-groups "
          f"pass Rayleigh test (p<0.05)")

    print(f"\n  Exp 2 — Ratio fit:   R²(θ only)={r2['loglog_r2']:.4f}  "
          f"R²(angle)={r2['helix_angle_vs_ratio_r2']:.5f}  "
          f"R²(tautology)={r2['tautology_r2']:.5f}")
    print(f"    Helix angle explains {r2['helix_angle_vs_ratio_r2']*100:.2f}% "
          f"of variance in log(p/q).")

    k10 = next(r for r in r3["by_K"] if r["K"] == 10)
    lift = k10["helix_exact_hr"] - k10["global_hr"]
    print(f"\n  Exp 3 — Search (K=10): helix={k10['helix_exact_hr']:.1%}  "
          f"global={k10['global_hr']:.1%}  lift={lift:+.1%}  "
          f"({k10['search_reduction_pct']:.0f}% of search space used)")
    print(f"    {r3['verdict']}")

    print()
    print("  ┌─────────────────────────────────────────────────────────────────┐")
    if abs(lift) < 0.02 and abs(z) <= 3:
        print("  │ CONCLUSION: Helix geometry is VISUALIZATION ONLY.              │")
        print("  │ p-values interleave uniformly on the helix — no angular        │")
        print("  │ clustering by factor.  The geometry encodes lopsidedness        │")
        print("  │ beautifully but carries no additional factoring signal.         │")
    elif abs(lift) >= 0.05 or abs(z) > 3:
        print("  │ CONCLUSION: Helix geometry IS PREDICTIVE.                       │")
        print("  │ Angular position carries statistically significant information  │")
        print("  │ about the smaller factor p.  Further investigation warranted.  │")
    else:
        print("  │ CONCLUSION: MARGINAL SIGNAL.  Helix provides weak but          │")
        print("  │ inconsistent predictive lift — geometry is primarily            │")
        print("  │ visualization with hints of structure worth exploring further.  │")
    print("  └─────────────────────────────────────────────────────────────────┘")
    print()


if __name__ == "__main__":
    main()
