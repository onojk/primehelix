from theta_sweep_generate import generate_semiprimes
from lopsided import summarize_lopsided


LIMITS = [10_000, 100_000, 1_000_000]
THETA = 0.25


print("\n=== LOPSIDED CONVERGENCE ===")
print(f"Definition: n = p*q is lopsided if min(p,q) <= n^{THETA}")
print()
print(f"{'N':>12} {'total':>12} {'lopsided':>12} {'balanced':>12} {'ratio':>10} {'percent':>10}")
print("-" * 74)

for limit in LIMITS:
    semiprimes = generate_semiprimes(limit)
    s = summarize_lopsided(semiprimes, theta=THETA)

    print(
        f"{limit:>12,} "
        f"{s['total']:>12,} "
        f"{s['lopsided']:>12,} "
        f"{s['balanced']:>12,} "
        f"{s['lopsided_ratio']:>10.4f} "
        f"{100 * s['lopsided_ratio']:>9.2f}%"
    )
