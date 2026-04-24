import math

EMPIRICAL_DATA = {
    1_000_000:   {"lopsided_pct": 73.18, "source": "primehelix scan [1, 1M)"},
    10_000_000:  {"lopsided_pct": 78.54, "source": "primehelix scan [1, 10M)"},
    100_000_000: {"lopsided_pct": 82.22, "source": "primehelix scan [1, 100M)"},
}

def pi_approx(x):
    if x < 2: return 0.0
    lx = math.log(x)
    return x / lx * (1 + 1/lx + 2/lx**2)

def semiprime_density(p, N):
    if p <= 1: return 0.0
    q_max = N / p
    q_min = max(p, 2.0)
    if q_max <= q_min: return 0.0
    return max(0.0, pi_approx(q_max) - pi_approx(q_min))

def integrate(lo, hi, N, n=2000):
    if lo >= hi or lo < 2: return 0.0
    log_lo, log_hi = math.log(max(lo,2.0)), math.log(hi)
    if log_hi <= log_lo: return 0.0
    total = 0.0
    prev_p = math.exp(log_lo)
    prev_f = semiprime_density(prev_p, N) / math.log(prev_p)
    for i in range(1, n+1):
        p = math.exp(log_lo + (log_hi - log_lo) * i / n)
        f = semiprime_density(p, N) / math.log(p)
        total += 0.5 * (prev_f + f) * (p - prev_p)
        prev_p, prev_f = p, f
    return total

def theoretical_lopsided(N, bit_gap=8):
    sqrt_N = math.sqrt(N)
    T = sqrt_N / (2 ** ((bit_gap + 1) / 2))
    lop = integrate(2.0, T, N)
    total = integrate(2.0, sqrt_N, N)
    if total <= 0: return 0.0
    return lop / total * 100

print()
print("=" * 72)
print("  LOPSIDED SEMIPRIME: THEORY vs EMPIRICAL")
print("=" * 72)
print(f"  {'N':>12}  {'Theory %':>10}  {'Empirical %':>12}  {'Delta':>8}  {'Verdict'}")
print("-" * 72)

results = []
for N, emp in sorted(EMPIRICAL_DATA.items()):
    theory = theoretical_lopsided(float(N))
    delta = emp["lopsided_pct"] - theory
    verdict = "theory holds" if abs(delta) < 1.0 else f"{'↑' if delta>0 else '↓'} +{abs(delta):.2f}pp offset"
    print(f"  {N:>12,}  {theory:>9.2f}%  {emp['lopsided_pct']:>11.2f}%  {delta:>+7.2f}pp  {verdict}")
    results.append((N, theory, emp["lopsided_pct"], delta))

print("=" * 72)

if len(results) >= 2:
    (N1,_,p1,d1),(N2,_,p2,d2) = results[0], results[-1]
    print()
    print("  CONVERGENCE MODEL: P(lopsided) ~ 1 - C / log(N)^alpha")
    r1, r2 = 1-p1/100, 1-p2/100
    if r2 > 0:
        alpha = math.log(r1/r2) / math.log(math.log(N2)/math.log(N1))
        C = r1 * math.log(N1)**alpha
        print(f"  alpha = {alpha:.3f},  C = {C:.3f}")
        print()
        print("  Predictions:")
        for Np in [1e9, 1e10, 1e12]:
            pred = (1 - C/math.log(Np)**alpha)*100
            print(f"    N = {Np:.0e}  =>  {pred:.1f}%")

print()
print("  DELTA TREND (should shrink toward 0):")
for N,th,emp,d in results:
    print(f"    N={N:>12,}  delta={d:+.2f}pp")
print()
