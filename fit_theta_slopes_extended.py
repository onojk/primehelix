import random
from bisect import bisect_left, bisect_right
from math import ceil, isqrt, log

import matplotlib.pyplot as plt
import numpy as np
from sympy import isprime, primepi


THETAS = [0.20, 0.25, 0.30, 0.35]
EXACT_LIMITS = [1_000, 10_000, 100_000, 1_000_000, 10_000_000, 100_000_000]
SAMPLE_LIMIT = 1_000_000_000
SAMPLES = 20_000


def sieve(limit):
    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"

    for i in range(2, isqrt(limit) + 1):
        if is_prime[i]:
            is_prime[i * i : limit + 1 : i] = b"\x00" * (
                ((limit - i * i) // i) + 1
            )

    return [i for i in range(2, limit + 1) if is_prime[i]]


def q_threshold_for_theta(p, theta):
    # p <= (p*q)^theta
    # log(p) <= theta(log(p)+log(q))
    # log(q) >= (1/theta - 1) log(p)
    exponent = (1.0 / theta) - 1.0
    return ceil(p ** exponent)


def exact_ratio(limit, theta, primes):
    total = 0
    lopsided = 0

    for p in primes:
        if p * p > limit:
            break

        q_min_total = p
        q_max = limit // p

        total += bisect_right(primes, q_max) - bisect_left(primes, q_min_total)

        q_min_lop = max(p, q_threshold_for_theta(p, theta))

        if q_min_lop <= q_max:
            lopsided += bisect_right(primes, q_max) - bisect_left(primes, q_min_lop)

    return lopsided / total if total else 0.0


def build_weighted_p_table(limit):
    p_primes = sieve(isqrt(limit))

    cumulative = []
    running_total = 0

    for p in p_primes:
        q_min = p
        q_max = limit // p

        if q_max < q_min:
            continue

        weight = int(primepi(q_max) - primepi(q_min - 1))

        if weight <= 0:
            continue

        running_total += weight
        cumulative.append((running_total, p, q_min, q_max))

    return cumulative, running_total


def random_prime_between(lo, hi):
    while True:
        q = random.randint(lo, hi)
        if q >= 2 and isprime(q):
            return q


def sample_pairs(limit, samples):
    cumulative, total_weight = build_weighted_p_table(limit)
    cumulative_keys = [row[0] for row in cumulative]

    pairs = []

    for i in range(samples):
        r = random.randint(1, total_weight)
        idx = bisect_left(cumulative_keys, r)

        _, p, q_min, q_max = cumulative[idx]
        q = random_prime_between(q_min, q_max)

        pairs.append((p, q))

        if (i + 1) % 5000 == 0:
            print(f"  sampled {i + 1:,}/{samples:,}")

    return pairs, total_weight


def pair_is_lopsided(p, q, theta):
    return log(p) <= theta * (log(p) + log(q))


def main():
    print("Sieving primes up to 50,000,000 for exact counts...")
    primes = sieve(max(EXACT_LIMITS) // 2)
    print(f"Primes found: {len(primes):,}")

    results = {theta: [] for theta in THETAS}

    print("\n=== EXACT RATIOS ===")
    for theta in THETAS:
        print(f"\nθ = {theta:.2f}")
        for limit in EXACT_LIMITS:
            ratio = exact_ratio(limit, theta, primes)
            results[theta].append((limit, ratio, "exact"))
            print(f"  N <= {limit:>12,}: {100 * ratio:6.2f}%")

    print(f"\nSampling N <= {SAMPLE_LIMIT:,}...")
    pairs, total_pairs = sample_pairs(SAMPLE_LIMIT, SAMPLES)
    print(f"Estimated total semiprime pairs: {total_pairs:,}")

    print("\n=== SAMPLED 1e9 RATIOS ===")
    for theta in THETAS:
        lopsided = sum(1 for p, q in pairs if pair_is_lopsided(p, q, theta))
        ratio = lopsided / SAMPLES
        results[theta].append((SAMPLE_LIMIT, ratio, "sample"))
        print(f"  θ = {theta:.2f}: {100 * ratio:6.2f}%")

    print("\n=== EXTENDED FIT Δ(N,θ) = a(θ) - b(θ)log(N) ===")
    print(f"{'theta':>8} {'a(theta)':>12} {'b(theta)':>12}")
    print("-" * 36)

    theta_values = []
    b_values = []

    plt.figure(figsize=(9, 5))

    for theta in THETAS:
        ns = np.array([row[0] for row in results[theta]], dtype=float)
        ps = np.array([row[1] for row in results[theta]], dtype=float)
        delta = 1.0 - ps
        log_n = np.log(ns)

        m, c = np.polyfit(log_n, delta, 1)
        a = c
        b = -m

        theta_values.append(theta)
        b_values.append(b)

        print(f"{theta:>8.2f} {a:>12.6f} {b:>12.6f}")

        plt.plot(ns, 100 * delta, marker="o", linewidth=2, label=f"θ = {theta:.2f}")

    plt.xscale("log")
    plt.xlabel("N cutoff")
    plt.ylabel("Δ(N,θ) = non-lopsided (%)")
    plt.title("Extended Δ(N,θ): exact 1e3–1e8 plus sampled 1e9")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("extended_theta_delta.png", dpi=200)
    print("\nSaved: extended_theta_delta.png")

    plt.figure(figsize=(8, 5))
    plt.plot(theta_values, b_values, marker="o", linewidth=2)
    plt.xlabel("θ")
    plt.ylabel("b(θ)")
    plt.title("Extended decay-rate parameter b(θ)")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("extended_theta_slope_fit.png", dpi=200)
    print("Saved: extended_theta_slope_fit.png")


if __name__ == "__main__":
    main()
