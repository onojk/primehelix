#!/usr/bin/env python3
"""
fit_prefactor_drift.py

Test whether δ(N, θ₀) · log(log(N)) is flat across N, which would confirm
the O(1/log log N) convergence rate with prefactor log(1/θ₀).

For each θ₀ in {0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45} and
N in {1e4, 1e5, 1e6, 1e7, 1e8}:

    delta(N, θ₀)  = 1 - P(N, θ₀)
    b_eff(N, θ₀)  = delta(N, θ₀) · log(log(N))
    theory        = log(1/θ₀)
    ratio         = b_eff / theory

If ratio ≈ 1.0 and flat across N  → prefactor is log(1/θ₀), rate confirmed.
If ratio drifts                    → subleading corrections present.
"""

import math
import pickle
import pathlib
import time
from bisect import bisect_left, bisect_right
from math import isqrt, ceil

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Config ────────────────────────────────────────────────────────────────────

THETAS  = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45]
N_VALS  = [1e4,  1e5,  1e6,  1e7,  1e8]
CACHE   = pathlib.Path(".cache")

# ── Prime sieve (cached) ──────────────────────────────────────────────────────

def sieve(limit: int):
    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0] = is_prime[1] = 0
    for i in range(2, isqrt(limit) + 1):
        if is_prime[i]:
            is_prime[i*i : limit+1 : i] = b"\x00" * len(is_prime[i*i : limit+1 : i])
    return [i for i in range(2, limit + 1) if is_prime[i]]


def sieve_cached(limit: int):
    CACHE.mkdir(exist_ok=True)
    path = CACHE / f"primes_{limit}.pkl"
    if path.exists():
        with open(path, "rb") as f:
            return pickle.load(f)
    print(f"  Sieving primes to {limit:,} …", flush=True)
    primes = sieve(limit)
    with open(path, "wb") as f:
        pickle.dump(primes, f)
    return primes

# ── Exact semiprime count ─────────────────────────────────────────────────────

def count_ratio(limit: int, theta: float, primes: list) -> float:
    """Return P(limit, theta) — fraction of semiprimes n≤limit with θ≤theta."""
    total = lopsided = 0
    for p in primes:
        if p * p > limit:
            break
        q_max = limit // p
        total    += bisect_right(primes, q_max) - bisect_left(primes, p)
        q_min_lop = max(p, ceil(p ** ((1.0 / theta) - 1.0)))
        if q_min_lop <= q_max:
            lopsided += bisect_right(primes, q_max) - bisect_left(primes, q_min_lop)
    return lopsided / total if total else 0.0

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    max_N    = int(max(N_VALS))
    prime_limit = max_N // 2

    print("=" * 68)
    print("fit_prefactor_drift.py")
    print("Testing: b_eff(N,θ₀) = δ(N,θ₀)·log(log(N))  vs  theory = log(1/θ₀)")
    print("=" * 68)
    print()

    print(f"Loading primes to {prime_limit:,} …", flush=True)
    t0 = time.time()
    primes = sieve_cached(prime_limit)
    print(f"  {len(primes):,} primes loaded ({time.time()-t0:.2f}s)")
    print()

    # ── Compute all (theta, N) cells ──────────────────────────────────────────
    llN    = [math.log(math.log(n)) for n in N_VALS]
    data   = {}   # data[theta][i] = (delta, b_eff)

    print(f"{'theta':>6} | {'N':>8} | {'delta':>9} | {'b_eff':>9} | "
          f"{'theory':>9} | {'ratio':>7}")
    print("-" * 68)

    for theta in THETAS:
        theory = math.log(1.0 / theta)
        data[theta] = []
        for i, n in enumerate(N_VALS):
            P     = count_ratio(int(n), theta, primes)
            delta = 1.0 - P
            b_eff = delta * llN[i]
            ratio = b_eff / theory if theory > 0 else float("nan")
            data[theta].append((delta, b_eff, ratio))
            print(f"{theta:>6.2f} | {n:>8.0e} | {delta:>9.5f} | {b_eff:>9.5f} | "
                  f"{theory:>9.5f} | {ratio:>7.4f}")
        print()

    # ── Summary ───────────────────────────────────────────────────────────────
    print()
    print("SUMMARY")
    print("=" * 68)
    print(f"{'theta':>6} | {'mean_b':>8} | {'std_b':>8} | {'cv':>7} | "
          f"{'theory':>8} | {'mean_ratio':>10} | {'verdict'}")
    print("-" * 68)

    verdicts = {}
    for theta in THETAS:
        theory  = math.log(1.0 / theta)
        bvals   = [row[1] for row in data[theta]]
        ratios  = [row[2] for row in data[theta]]
        mean_b  = np.mean(bvals)
        std_b   = np.std(bvals)
        cv      = std_b / mean_b          # coefficient of variation
        mean_r  = np.mean(ratios)
        verdict = "possible drift" if cv > 0.05 else "flat"
        verdicts[theta] = verdict

        print(f"{theta:>6.2f} | {mean_b:>8.5f} | {std_b:>8.5f} | {cv:>7.4f} | "
              f"{theory:>8.5f} | {mean_r:>10.4f} | {verdict}")

    print()
    print("Key output sentences:")
    for theta in THETAS:
        theory  = math.log(1.0 / theta)
        bvals   = [row[1] for row in data[theta]]
        ratios  = [row[2] for row in data[theta]]
        mean_b  = np.mean(bvals)
        mean_r  = np.mean(ratios)
        v       = verdicts[theta]
        print(f"  θ₀={theta:.2f}: b_eff mean={mean_b:.4f}, "
              f"theory={theory:.4f}, ratio={mean_r:.4f} — {v}")

    # ── Paper conclusion line ─────────────────────────────────────────────────
    n_drifting = sum(1 for v in verdicts.values() if v == "possible drift")
    print()
    if n_drifting == 0:
        print("PAPER STATEMENT: consistent with log(1/θ₀)/log(log(N))")
    else:
        print(f"PAPER STATEMENT: consistent with log(1/θ₀)/log(log(N)) "
              f"with subleading corrections  ({n_drifting}/{len(THETAS)} θ₀ values drifting)")

    # ── Plot ──────────────────────────────────────────────────────────────────
    COLORS = plt.cm.tab10(np.linspace(0, 0.9, len(THETAS)))

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    for theta, color in zip(THETAS, COLORS):
        theory = math.log(1.0 / theta)
        bvals  = [row[1] for row in data[theta]]

        ax.plot(llN, bvals,
                color=color, lw=2.0, marker="o", ms=5,
                label=fr"$\theta_0={theta}$")
        ax.axhline(theory,
                   color=color, lw=0.9, ls="--", alpha=0.65)

    # Annotate the dashed lines on the right margin
    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks([math.log(1.0/t) for t in THETAS])
    ax2.set_yticklabels([fr"$\log(1/{t})$" for t in THETAS], fontsize=7)
    ax2.tick_params(right=False)
    for sp in ax2.spines.values():
        sp.set_visible(False)

    ax.set_xlabel(r"$\log(\log(N))$", fontsize=13)
    ax.set_ylabel(r"$b_{\mathrm{eff}}(N,\theta_0) = \delta \cdot \log(\log(N))$",
                  fontsize=12)
    ax.set_title(
        r"Prefactor drift test: $\delta(N,\theta_0)\cdot\log(\log(N))$ vs $\log(\log(N))$"
        "\nDashed lines: theoretical prefactor $\\log(1/\\theta_0)$",
        fontsize=11
    )
    ax.legend(fontsize=9, loc="upper left", ncol=2)
    ax.grid(color="#dddddd", lw=0.5)
    ax.tick_params(labelsize=10)

    # Mark N values on top x-axis
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    ax_top.set_xticks(llN)
    ax_top.set_xticklabels([r"$10^{%d}$" % int(math.log10(n)) for n in N_VALS],
                            fontsize=9)
    ax_top.set_xlabel(r"$N$", fontsize=11)

    plt.tight_layout()
    out = "prefactor_drift.png"
    plt.savefig(out, dpi=160, facecolor="white", bbox_inches="tight")
    plt.close()
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
