#!/usr/bin/env python3
"""
spectral_fingerprint.py

Computes C(t) = |Σ n^(it)| / C(0) for semiprimes, primes, and random integers.
Plots all three and marks the first 10 nontrivial Riemann zeta zeros.
"""
import array
import math
import random
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Config ────────────────────────────────────────────────────────────────────
N_MAX  = 500_000
N_T    = 2000
T_MAX  = 50.0
CHUNK  = 50       # t-values per batch (memory: CHUNK × N × 8 bytes)

ZETA_ZEROS = [
    14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
    37.586178, 40.918719, 43.327073, 48.005151, 49.773832,
]

# ── SPF sieve ─────────────────────────────────────────────────────────────────

def spf_sieve(limit: int):
    spf = array.array("i", range(limit + 1))
    i = 2
    while i * i <= limit:
        if spf[i] == i:
            for j in range(i * i, limit + 1, i):
                if spf[j] == j:
                    spf[j] = i
        i += 1
    return spf

# ── Build number sets ─────────────────────────────────────────────────────────

print("Building SPF sieve …", flush=True)
t0  = time.time()
spf = spf_sieve(N_MAX)
print(f"  done in {time.time()-t0:.2f}s")

print("Collecting semiprimes …", flush=True)
semiprimes = []
for n in range(4, N_MAX + 1):
    p = spf[n]
    if p == n:
        continue
    q = n // p
    if spf[q] != q:
        continue
    semiprimes.append(n)

print("Collecting primes …", flush=True)
primes = [n for n in range(2, N_MAX + 1) if spf[n] == n]

N_semi  = len(semiprimes)
N_prime = len(primes)

print("Building random set …", flush=True)
rng        = random.Random(42)
random_set = sorted(rng.sample(range(4, N_MAX + 1), N_semi))
N_rand     = len(random_set)

print(f"\n  Semiprimes : {N_semi:,}")
print(f"  Primes     : {N_prime:,}")
print(f"  Random     : {N_rand:,}  (uniform in [4, {N_MAX:,}])")

# ── Spectral function ─────────────────────────────────────────────────────────

def compute_C(log_n: np.ndarray, t_arr: np.ndarray) -> np.ndarray:
    """C(t) = |Σ exp(i·t·log n)|, computed in memory-efficient chunks."""
    C = np.empty(len(t_arr))
    for i in range(0, len(t_arr), CHUNK):
        t_chunk = t_arr[i : i + CHUNK]           # (chunk,)
        phases  = np.outer(t_chunk, log_n)        # (chunk, N)
        cs      = np.cos(phases).sum(axis=1)
        ss      = np.sin(phases).sum(axis=1)
        C[i : i + CHUNK] = np.sqrt(cs * cs + ss * ss)
    return C

log_semi   = np.log(np.array(semiprimes,  dtype=np.float64))
log_prime  = np.log(np.array(primes,      dtype=np.float64))
log_random = np.log(np.array(random_set,  dtype=np.float64))

t_arr = np.linspace(0.0, T_MAX, N_T)

print("\nComputing C(t) …")
t0 = time.time()

print("  semiprimes …", end=" ", flush=True)
C_semi   = compute_C(log_semi,   t_arr)
print(f"{time.time()-t0:.1f}s")

t1 = time.time()
print("  primes …",     end=" ", flush=True)
C_prime  = compute_C(log_prime,  t_arr)
print(f"{time.time()-t1:.1f}s")

t1 = time.time()
print("  random …",     end=" ", flush=True)
C_random = compute_C(log_random, t_arr)
print(f"{time.time()-t1:.1f}s")

# Normalise: C(0) = N exactly (all e^0 = 1)
C_semi_n   = C_semi   / N_semi
C_prime_n  = C_prime  / N_prime
C_random_n = C_random / N_rand

print(f"\nTotal compute time: {time.time()-t0:.1f}s")

# ── Print values at zeta zeros ────────────────────────────────────────────────

print()
print(f"{'Zeta zero':>12}  {'C_semi':>10}  {'C_prime':>10}  {'C_random':>10}")
print("─" * 50)
for z in ZETA_ZEROS:
    idx = min(int(np.searchsorted(t_arr, z)), N_T - 1)
    print(f"{z:>12.6f}  {C_semi_n[idx]:>10.6f}  "
          f"{C_prime_n[idx]:>10.6f}  {C_random_n[idx]:>10.6f}")

# ── Summary stats ─────────────────────────────────────────────────────────────

# Which set has higher C at zeta zeros vs generic t?
zero_idxs    = [min(int(np.searchsorted(t_arr, z)), N_T-1) for z in ZETA_ZEROS]
nozero_mask  = np.ones(N_T, bool)
nozero_mask[zero_idxs] = False

print()
print("Mean C(t) at zeta zeros vs rest-of-range:")
for label, C_n in [("semiprimes", C_semi_n), ("primes", C_prime_n), ("random", C_random_n)]:
    at_zeros = C_n[zero_idxs].mean()
    elsewhere = C_n[nozero_mask & (t_arr > 10)].mean()   # t>10 to skip the large-t=0 peak
    print(f"  {label:<12}  at zeros={at_zeros:.5f}  elsewhere={elsewhere:.5f}  "
          f"ratio={at_zeros/elsewhere:.3f}")

# ── Plot ──────────────────────────────────────────────────────────────────────

DARK_BG   = "#0a0a12"
PANEL_BG  = "#12121e"
SPINE_COL = "#2a2a40"
TEXT_COL  = "#ccccee"

fig, ax = plt.subplots(figsize=(15, 6))
fig.patch.set_facecolor(DARK_BG)
ax.set_facecolor(PANEL_BG)
ax.tick_params(colors=TEXT_COL, labelsize=9)
for s in ax.spines.values():
    s.set_color(SPINE_COL)

ax.plot(t_arr, C_semi_n,   color="#3ecf96", lw=1.1, label=f"Semiprimes (N={N_semi:,})", alpha=0.95)
ax.plot(t_arr, C_prime_n,  color="#9f8ff8", lw=1.1, label=f"Primes     (N={N_prime:,})", alpha=0.95)
ax.plot(t_arr, C_random_n, color="#888888", lw=0.8, label=f"Random     (N={N_rand:,})",  alpha=0.65)

# Zeta zero verticals
xform = ax.get_xaxis_transform()   # x=data, y=axes fraction
for i, z in enumerate(ZETA_ZEROS):
    ax.axvline(z, color="#f07050", lw=0.9, ls="--", alpha=0.70)
    ax.text(z + 0.18, 0.98, f"ζ{i+1}", color="#f07050",
            fontsize=7, va="top", transform=xform)

ax.set_xlabel("t", color=TEXT_COL, fontsize=11)
ax.set_ylabel("C(t) / C(0)", color=TEXT_COL, fontsize=11)
ax.set_title(
    r"Spectral fingerprint:  $C(t) = |\sum n^{it}| \,/\, N$   "
    "for semiprimes, primes, and random integers up to 500 000\n"
    "Dashed orange lines: first 10 nontrivial Riemann zeta zeros",
    color=TEXT_COL, fontsize=11, pad=10,
)
ax.legend(facecolor=PANEL_BG, edgecolor=SPINE_COL, labelcolor=TEXT_COL,
          fontsize=9, loc="upper right")
ax.set_xlim(0, T_MAX)
ax.set_ylim(bottom=0)
ax.grid(color=SPINE_COL, lw=0.45, alpha=0.55)

plt.tight_layout()
fname = "spectral_fingerprint.png"
plt.savefig(fname, dpi=140, facecolor=DARK_BG)
plt.close()
print(f"\nSaved: {fname}")
