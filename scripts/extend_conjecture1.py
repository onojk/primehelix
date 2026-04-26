#!/usr/bin/env python3
"""
extend_conjecture1.py

Task A: Extend Conjecture 1 evidence to N=1e9 using memory-efficient
        pair-counting (equivalent to SPF sieve but avoids the ~4 GB
        uint32 array; see verify_paper_tables.py docstring).
        Prime list up to N_max/2 = 5e8 via odd-only numpy sieve (~238 MB
        peak during build, ~105 MB retained).

Task B: Fit six single-parameter subleading correction models to the
        residuals  residual(N, θ₀) = δ(N,θ₀)·log log N − f(θ₀).

Outputs
-------
  results/conjecture1_extended.csv
  results/subleading_models.csv
  stabilisation_extended.png
"""

import csv
import math
import pickle
import pathlib
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Constants ──────────────────────────────────────────────────────────────────

THETA_VALS  = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45]
THETA_PLOT  = [0.20, 0.25, 0.30, 0.35]
CHECKPOINTS = [int(1e8), int(2e8), int(5e8), int(1e9)]
CACHE_DIR   = pathlib.Path(".cache")
RESULTS_DIR = pathlib.Path("results")

F_THETA = {t: math.log(1.0/t) - t - math.log(2) + 0.5 for t in THETA_VALS}
EXPECTED_1E8 = 17_427_258          # validated semiprime count at N=1e8

# ── Prime list builder ─────────────────────────────────────────────────────────

def build_prime_list(limit: int) -> np.ndarray:
    """Sorted numpy uint32 array of all primes ≤ limit via odd-only sieve."""
    CACHE_DIR.mkdir(exist_ok=True)
    path = CACHE_DIR / f"primes_np_{limit}.pkl"
    if path.exists():
        print(f"  Loading prime list (limit={limit:,}) from cache …")
        with open(path, "rb") as f:
            return pickle.load(f)

    print(f"  Sieving primes to {limit:,} (odd-only numpy) …", end=" ", flush=True)
    t0 = time.time()

    # Odd index i represents number 2i + 1; index 0 → 1, index 1 → 3, ...
    size = limit // 2 + 1
    is_prime_odd = np.ones(size, dtype=bool)   # ~238 MB for limit=5e8
    is_prime_odd[0] = False                    # 1 is not prime
    for p in range(3, int(limit**0.5) + 1, 2):
        if is_prime_odd[p // 2]:
            is_prime_odd[p * p // 2 :: p] = False

    # Extract: is_prime_odd[i] True → prime = 2i+1; skip index 0 (= 1)
    odd_primes = (np.flatnonzero(is_prime_odd[1:]) + 1) * 2 + 1
    primes = np.concatenate([np.array([2], dtype=np.uint32),
                              odd_primes.astype(np.uint32)])
    del is_prime_odd   # free ~238 MB
    elapsed = time.time() - t0
    print(f"{len(primes):,} primes  ({elapsed:.1f}s)")

    with open(path, "wb") as f:
        pickle.dump(primes, f)
    return primes

# ── Pair-counting P(N, θ₀) ────────────────────────────────────────────────────

def count_P(N: int, theta: float, primes: np.ndarray):
    """
    Exact pair-counting:  iterate p over primes ≤ √N, count valid q via
    numpy searchsorted.  Returns (P, lopsided_count, total_sp_count).

    total_sp_count is independent of theta (all unordered semiprime pairs).
    """
    sqrt_N  = int(N**0.5)
    p_end   = int(np.searchsorted(primes, sqrt_N, side='right'))

    total    = 0
    lopsided = 0

    for p in primes[:p_end]:
        q_max = N // p
        if q_max < p:
            break                                  # all subsequent p also > √N
        q_max_idx = int(np.searchsorted(primes, q_max, side='right'))
        p_idx     = int(np.searchsorted(primes, p,     side='left'))
        cnt_all   = q_max_idx - p_idx
        total    += cnt_all

        # θ(pq) ≤ θ₀  ⟺  p ≤ (pq)^θ₀  ⟺  q ≥ p^(1/θ₀ − 1)
        q_min_raw = p ** (1.0 / theta - 1.0)
        q_min     = max(int(p), math.ceil(q_min_raw))
        if q_min <= q_max:
            q_min_idx  = int(np.searchsorted(primes, q_min, side='left'))
            lopsided  += q_max_idx - q_min_idx

    P = lopsided / total if total else 0.0
    return P, lopsided, total

# ── Checkpoint validation ──────────────────────────────────────────────────────

def validate_n1e8(primes: np.ndarray):
    """Verify total semiprime count at N=1e8 matches SPF-sieve reference."""
    N = int(1e8)
    _, _, total = count_P(N, 0.10, primes)   # theta doesn't affect total
    status = "OK" if total == EXPECTED_1E8 else f"FAIL (got {total:,})"
    print(f"  N=1e8 semiprime count: expected {EXPECTED_1E8:,}  got {total:,}  {status}")
    if total != EXPECTED_1E8:
        raise RuntimeError("Semiprime count mismatch at N=1e8 — debug pair-counting.")

# ── Task A ─────────────────────────────────────────────────────────────────────

def task_A(primes: np.ndarray):
    print("=" * 70)
    print("Task A: Extended Conjecture 1 data (N up to 1e9)")
    print("=" * 70)

    rows = []
    print(f"  {'N':>11}  {'θ₀':>5}  {'total_sp':>12}  {'P':>9}  {'δ':>9}  "
          f"{'llogN':>7}  {'b_eff':>8}")
    print("  " + "-" * 70)

    for N in CHECKPOINTS:
        llog_N = math.log(math.log(N))
        # total_sp is theta-independent; compute once for this N
        _, _, total_sp = count_P(N, 0.10, primes)

        for t in THETA_VALS:
            P, lopsided, _ = count_P(N, t, primes)
            delta = 1.0 - P
            b_eff = delta * llog_N
            f_val = F_THETA[t]
            print(f"  {N:>11,}  {t:>5.2f}  {total_sp:>12,}  {P:>9.6f}  {delta:>9.6f}  "
                  f"{llog_N:>7.4f}  {b_eff:>8.4f}")
            rows.append({
                "N":        N,
                "theta0":   t,
                "total_sp": total_sp,
                "lopsided": lopsided,
                "P":        round(P, 7),
                "delta":    round(delta, 7),
                "llogN":    round(llog_N, 6),
                "product":  round(b_eff, 6),
            })
        print()

    RESULTS_DIR.mkdir(exist_ok=True)
    fname = RESULTS_DIR / "conjecture1_extended.csv"
    with open(fname, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"  Saved: {fname}\n")
    return rows

# ── Extended stabilisation plot ────────────────────────────────────────────────

def plot_extended_stabilisation(rows):
    COLORS = plt.cm.tab10(np.linspace(0, 0.5, len(THETA_PLOT)))
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.set_facecolor("white")
    ax.grid(color="#dddddd", lw=0.5)

    for t, color in zip(THETA_PLOT, COLORS):
        t_rows = sorted([r for r in rows if r["theta0"] == t], key=lambda r: r["N"])
        N_plot = [r["N"] for r in t_rows]
        b_plot = [r["product"] for r in t_rows]
        f_val  = F_THETA[t]
        ax.semilogx(N_plot, b_plot, "o-", color=color, lw=1.8, ms=8, zorder=5,
                    label=fr"$\theta_0={t}$")
        ax.axhline(f_val, color=color, lw=0.9, ls="--", alpha=0.5)

    ax.set_xlabel(r"$N$", fontsize=13)
    ax.set_ylabel(r"$b_{\mathrm{eff}} = \delta\cdot\log\log N$", fontsize=13)
    ax.set_title(
        r"$b_{\mathrm{eff}}(N,\theta_0)$ extended to $N=10^9$  "
        r"(dashed = heuristic $f(\theta_0)$)",
        fontsize=11)
    ax.legend(fontsize=10, ncol=2)
    plt.tight_layout()
    outfile = "stabilisation_extended.png"
    plt.savefig(outfile, dpi=155, facecolor="white", bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outfile}")

# ── Task B ─────────────────────────────────────────────────────────────────────

def task_B(task_A_rows):
    print("=" * 70)
    print("Task B: Subleading correction model fitting")
    print("=" * 70)

    # Build residual dict: (N, theta) → residual = b_eff - f(theta)
    residual_dict = {}
    for r in task_A_rows:
        residual_dict[(r["N"], r["theta0"])] = r["product"] - F_THETA[r["theta0"]]

    # Flatten to arrays for cross-theta fitting
    N_arr     = []
    theta_arr = []
    resid_arr = []
    for N in CHECKPOINTS:
        for t in THETA_VALS:
            v = residual_dict.get((N, t))
            if v is not None and not math.isnan(v):
                N_arr.append(N)
                theta_arr.append(t)
                resid_arr.append(v)

    N_arr     = np.array(N_arr,     dtype=np.float64)
    theta_arr = np.array(theta_arr, dtype=np.float64)
    resid_arr = np.array(resid_arr, dtype=np.float64)
    llogN_arr = np.log(np.log(N_arr))
    f_arr     = np.array([F_THETA[t] for t in theta_arr])

    def fit_model(predictor, target):
        """Single-parameter OLS: target ≈ A * predictor (no intercept)."""
        A     = np.dot(predictor, target) / np.dot(predictor, predictor)
        fitted = A * predictor
        ss_res = np.sum((target - fitted)**2)
        ss_tot = np.sum((target - target.mean())**2)
        r2     = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        max_ae = float(np.max(np.abs(target - fitted)))
        return float(A), r2, max_ae

    models = {}

    # S1: residual = A · θ
    A, r2, mae = fit_model(theta_arr, resid_arr)
    models["S1"] = {"form": "A·θ₀", "A": A, "R2": r2, "max_abs": mae}

    # S2: residual = A · θ²
    A, r2, mae = fit_model(theta_arr**2, resid_arr)
    models["S2"] = {"form": "A·θ₀²", "A": A, "R2": r2, "max_abs": mae}

    # S3: residual = A · θ · log(1/θ)
    A, r2, mae = fit_model(theta_arr * np.log(1.0/theta_arr), resid_arr)
    models["S3"] = {"form": "A·θ₀·log(1/θ₀)", "A": A, "R2": r2, "max_abs": mae}

    # S4: residual = A · (log 2 − θ)
    A, r2, mae = fit_model(math.log(2) - theta_arr, resid_arr)
    models["S4"] = {"form": "A·(log2−θ₀)", "A": A, "R2": r2, "max_abs": mae}

    # S5: residual = A / log log N  (fit globally)
    A, r2, mae = fit_model(1.0 / llogN_arr, resid_arr)
    models["S5"] = {"form": "A/loglogN", "A": A, "R2": r2, "max_abs": mae}

    # S6: residual = A · f(θ)
    A, r2, mae = fit_model(f_arr, resid_arr)
    models["S6"] = {"form": "A·f(θ₀)", "A": A, "R2": r2, "max_abs": mae}

    print(f"  {'Model':>6}  {'Form':>18}  {'A':>8}  {'R²':>7}  {'maxAbsRes':>10}")
    print("  " + "-" * 58)

    best_model = None
    best_r2    = -np.inf
    rows_out   = []

    for name, m in models.items():
        cand = " ← candidate" if m["R2"] > 0.97 else ""
        print(f"  {name:>6}  {m['form']:>18}  {m['A']:>8.4f}  {m['R2']:>7.4f}  "
              f"{m['max_abs']:>10.4f}{cand}")
        if m["R2"] > best_r2:
            best_r2    = m["R2"]
            best_model = name
        rows_out.append({
            "model":     name,
            "form":      m["form"],
            "A":         round(m["A"], 6),
            "R2":        round(m["R2"], 5),
            "max_abs_residual": round(m["max_abs"], 6),
            "candidate_R2gt097": m["R2"] > 0.97,
        })

    print()
    print(f"  Best: {best_model}  (R²={best_r2:.4f})")
    if best_r2 > 0.97:
        m = models[best_model]
        print(f"  Conjecture 1b candidate:")
        print(f"    δ·loglogN = f(θ) + ({m['A']:.4f})·[{m['form']}] + o(1/loglogN)")
    else:
        print("  No single-parameter model achieves R² > 0.97.")
        print("  Subleading correction structure remains open.")

    fname = RESULTS_DIR / "subleading_models.csv"
    with open(fname, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()))
        w.writeheader(); w.writerows(rows_out)
    print(f"\n  Saved: {fname}")

    return models, best_model, best_r2

# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    N_prime_max = max(CHECKPOINTS) // 2   # primes up to 5e8 for pair-counting at N=1e9

    print("=" * 70)
    print("extend_conjecture1.py — Tasks A + B")
    print("=" * 70)
    print()

    t0 = time.time()

    primes = build_prime_list(N_prime_max)
    print(f"  Prime list: {len(primes):,} primes, "
          f"{primes.nbytes // (1024*1024):.0f} MB\n")

    print("Checkpoint validation:")
    validate_n1e8(primes)
    print()

    rows_A = task_A(primes)
    plot_extended_stabilisation(rows_A)
    print()

    models, best_model, best_r2 = task_B(rows_A)

    print()
    print("=" * 70)
    print(f"Done.  Wall time: {time.time()-t0:.1f}s")
    print("=" * 70)
    return models, best_model, best_r2


if __name__ == "__main__":
    main()
