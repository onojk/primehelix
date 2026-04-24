#!/usr/bin/env python3
"""
Multiplicative-phase helix experiment.
Tests whether phi_t(n) = t*log(n) mod 2pi clusters semiprimes
by small factor p — especially at zeta-zero t values.
"""
import math
import csv
from collections import defaultdict

N_MAX = 500_000

# t values to test: generic + first 5 nontrivial Riemann zeta zeros
T_VALUES = [
    0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0,
    14.134725,  # zeta zero 1
    21.022040,  # zeta zero 2
    25.010858,  # zeta zero 3
    30.424876,  # zeta zero 4
    32.935062,  # zeta zero 5
]

def spf_sieve(limit):
    spf = list(range(limit + 1))
    for i in range(2, int(limit**0.5) + 1):
        if spf[i] == i:
            for j in range(i*i, limit + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf

def circ_R(angles):
    n = len(angles)
    if n == 0: return 0.0
    s = sum(math.sin(a) for a in angles)
    c = sum(math.cos(a) for a in angles)
    return math.sqrt(s*s + c*c) / n

def shuffled_baseline(groups, trials=500):
    """Mean R under random permutation of p labels."""
    all_angles = [a for angles in groups.values() for a in angles]
    import random
    sizes = [len(v) for v in groups.values()]
    rs = []
    for _ in range(trials):
        random.shuffle(all_angles)
        idx = 0
        r_vals = []
        for sz in sizes:
            r_vals.append(circ_R(all_angles[idx:idx+sz]))
            idx += sz
        rs.append(sum(r_vals)/len(r_vals))
    return sum(rs)/len(rs), (sum((r-sum(rs)/len(rs))**2 for r in rs)/len(rs))**0.5

print("Building SPF sieve...")
spf = spf_sieve(N_MAX)

print("Collecting semiprimes...")
semis = []
for n in range(4, N_MAX + 1):
    p = spf[n]
    if p == n: continue
    q = n // p
    if spf[q] != q: continue  # q not prime
    if p == q: continue
    theta = math.log(p) / math.log(n)
    semis.append((n, p, q, theta))

print(f"Found {len(semis):,} semiprimes up to {N_MAX:,}")
print()

rows = []
print(f"{'t':>12}  {'groups':>6}  {'obs_R':>8}  {'base_R':>8}  {'z-score':>8}  {'verdict'}")
print("-" * 65)

for t in T_VALUES:
    groups = defaultdict(list)
    for n, p, q, theta in semis:
        phi = (t * math.log(n)) % (2 * math.pi)
        groups[p].append(phi)

    valid = {p: a for p, a in groups.items() if len(a) >= 20}
    if not valid: continue

    obs_R = sum(circ_R(a) for a in valid.values()) / len(valid)
    base_mean, base_std = shuffled_baseline(valid, trials=200)
    z = (obs_R - base_mean) / base_std if base_std > 0 else 0

    # per-p breakdown for top clustered p values
    per_p = sorted([(p, circ_R(a), len(a)) for p, a in valid.items()],
                   key=lambda x: x[1], reverse=True)[:5]

    verdict = "SIGNAL" if z > 3 else ("weak" if z > 1.5 else "noise")
    print(f"{t:>12.6f}  {len(valid):>6}  {obs_R:>8.5f}  {base_mean:>8.5f}  {z:>+8.2f}  {verdict}")
    rows.append((t, len(valid), obs_R, base_mean, z, verdict))

print()
print("TOP t VALUES BY Z-SCORE:")
for row in sorted(rows, key=lambda x: x[4], reverse=True)[:5]:
    print(f"  t={row[0]:.6f}  z={row[4]:+.2f}  {row[5]}")

print()
print("ZETA ZEROS vs GENERIC:")
zeta_ts = {14.134725, 21.022040, 25.010858, 30.424876, 32.935062}
zeta_rows = [r for r in rows if r[0] in zeta_ts]
generic_rows = [r for r in rows if r[0] not in zeta_ts]
if zeta_rows and generic_rows:
    zeta_z = sum(r[4] for r in zeta_rows)/len(zeta_rows)
    generic_z = sum(r[4] for r in generic_rows)/len(generic_rows)
    print(f"  Mean z-score at zeta zeros:  {zeta_z:+.3f}")
    print(f"  Mean z-score at generic t:   {generic_z:+.3f}")
    print(f"  Difference: {zeta_z - generic_z:+.3f}")
    if zeta_z > generic_z + 1.0:
        print("  >> ZETA ZEROS OUTPERFORM — potential resonance signal")
    else:
        print("  >> No significant difference between zeta and generic t")

with open("multiplicative_phase_scan.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["t", "num_p_groups", "obs_R", "base_R", "z_score", "verdict"])
    for row in rows:
        w.writerow(row)

print()
print("Results saved to multiplicative_phase_scan.csv")
