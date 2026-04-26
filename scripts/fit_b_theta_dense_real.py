import os
import random
from bisect import bisect_left
from math import ceil, isqrt, log

import matplotlib.pyplot as plt
import numpy as np
from sympy import isprime, primepi

from exact_lopsided_counts import sieve_cached


# Deterministic
os.environ["MPLBACKEND"] = "Agg"
random.seed(0)
np.random.seed(0)


THETAS = [0.15, 0.20, 0.25, 0.30, 0.35, 0.40]

EXACT_LIMITS = [
    1_000,
    10_000,
    100_000,
    1_000_000,
    10_000_000,
    100_000_000,
]

SAMPLE_LIMIT = 1_000_000_000
SAMPLES = 20_000


def q_threshold_for_theta(p, theta):
    exponent = (1.0 / theta) - 1.0
    return ceil(p ** exponent)


def exact_ratio(limit, theta, primes):
    from bisect import bisect_left, bisect_right

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


def sieve_small(limit):
    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"

    for i in range(2, isqrt(limit) + 1):
        if is_prime[i]:
            is_prime[i * i : limit + 1 : i] = b"\x00" * (
                ((limit - i * i) // i) + 1
            )

    return [i for i in range(2, limit + 1) if is_prime[i]]


def build_weighted_p_table(limit):
    p_primes = sieve_small(isqrt(limit))

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


def r_squared(y_true, y_pred):
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    return 1.0 - ss_res / ss_tot


def main():
    print("\n=== DENSE REAL θ ANALYSIS ===")
    print(f"θ values: {THETAS}")
    print(f"Exact limits: {EXACT_LIMITS}")
    print(f"Sample limit: {SAMPLE_LIMIT:,}")
    print(f"Samples: {SAMPLES:,}\n")

    print("Loading / building prime cache for exact counts...")
    primes = sieve_cached(max(EXACT_LIMITS) // 2)
    print(f"Primes loaded: {len(primes):,}")

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

    theta_values = []
    b_values = []
    a_values = []

    print("\n=== FIT Δ(N,θ) = a(θ) - b(θ)log(N) ===")
    print(f"{'theta':>8} {'a(theta)':>12} {'b(theta)':>12}")
    print("-" * 36)

    for theta in THETAS:
        ns = np.array([row[0] for row in results[theta]], dtype=float)
        ps = np.array([row[1] for row in results[theta]], dtype=float)

        delta = 1.0 - ps
        log_n = np.log(ns)

        m, c = np.polyfit(log_n, delta, 1)

        a = c
        b = -m

        theta_values.append(theta)
        a_values.append(a)
        b_values.append(b)

        print(f"{theta:>8.2f} {a:>12.6f} {b:>12.6f}")

    theta_arr = np.array(theta_values)
    b_arr = np.array(b_values)

    # Model A: b = k(1-theta)^alpha
    x_a = np.log(1.0 - theta_arr)
    y = np.log(b_arr)

    alpha, log_k = np.polyfit(x_a, y, 1)
    k = np.exp(log_k)
    b_fit_a = k * (1.0 - theta_arr) ** alpha
    r2_a = r_squared(b_arr, b_fit_a)

    # Model B: b = c / theta^gamma
    x_b = np.log(theta_arr)
    m_b, log_c = np.polyfit(x_b, y, 1)
    gamma = -m_b
    c = np.exp(log_c)
    b_fit_b = c / (theta_arr ** gamma)
    r2_b = r_squared(b_arr, b_fit_b)

    print("\n=== FIT b(θ) MODELS ===\n")

    print("Model A: b(θ) = k * (1 - θ)^α")
    print(f"k     = {k:.10f}")
    print(f"alpha = {alpha:.6f}")
    print(f"R²    = {r2_a:.6f}")
    print()

    print("Model B: b(θ) = c / θ^γ")
    print(f"c     = {c:.10f}")
    print(f"gamma = {gamma:.6f}")
    print(f"R²    = {r2_b:.6f}")
    print()

    print("θ * b(θ) diagnostic:")
    for theta, b in zip(theta_arr, b_arr):
        print(f"  θ={theta:.2f} | b={b:.8f} | θ*b={theta*b:.8f}")

    # Plot 1: dense b(theta)
    theta_smooth = np.linspace(0.14, 0.42, 400)
    b_smooth_a = k * (1.0 - theta_smooth) ** alpha
    b_smooth_b = c / (theta_smooth ** gamma)

    plt.figure(figsize=(8, 5))
    plt.scatter(theta_arr, b_arr, label="Observed b(θ)", zorder=5)
    plt.plot(theta_smooth, b_smooth_a, label="k(1-θ)^α")
    plt.plot(theta_smooth, b_smooth_b, linestyle="--", label="c/θ^γ")
    plt.xlabel("θ")
    plt.ylabel("b(θ)")
    plt.title("Dense real fit for decay-rate parameter b(θ)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("b_theta_dense_real_fit.png", dpi=200)

    # Plot 2: theta * b(theta)
    plt.figure(figsize=(8, 5))
    plt.plot(theta_arr, theta_arr * b_arr, marker="o", linewidth=2)
    plt.xlabel("θ")
    plt.ylabel("θ · b(θ)")
    plt.title("Diagnostic: θ · b(θ)")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("theta_times_b_real.png", dpi=200)

    # Plot 3: log-log b vs theta
    plt.figure(figsize=(8, 5))
    plt.plot(np.log(theta_arr), np.log(b_arr), marker="o", linewidth=2)
    plt.xlabel("log(θ)")
    plt.ylabel("log(b(θ))")
    plt.title("Log-log diagnostic: b(θ) vs θ")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("loglog_b_theta_real.png", dpi=200)

    print("\nSaved:")
    print(" - b_theta_dense_real_fit.png")
    print(" - theta_times_b_real.png")
    print(" - loglog_b_theta_real.png")


if __name__ == "__main__":
    main()
