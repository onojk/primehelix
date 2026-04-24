import os
import random
import numpy as np
import matplotlib.pyplot as plt

from bisect import bisect_left, bisect_right
from math import ceil, isqrt


# Ensure deterministic behavior (important for CI)
os.environ["MPLBACKEND"] = "Agg"
random.seed(0)
np.random.seed(0)


# N values (exact computation)
N = [10_000, 100_000, 1_000_000, 10_000_000, 100_000_000]

THETAS = [0.20, 0.25, 0.30, 0.35]


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
    exponent = (1.0 / theta) - 1.0
    return ceil(p ** exponent)


def compute_ratio(limit, theta, primes):
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

    return 100.0 * lopsided / total


def main():
    max_limit = max(N)
    prime_limit = max_limit // 2

    print(f"Sieve up to {prime_limit:,}...")
    primes = sieve(prime_limit)
    print(f"Primes found: {len(primes):,}")

    plt.figure(figsize=(9, 5))

    for theta in THETAS:
        print(f"\nComputing θ = {theta:.2f}")
        y_values = []

        for limit in N:
            ratio = compute_ratio(limit, theta, primes)
            y_values.append(ratio)
            print(f"  N <= {limit:,}: {ratio:.2f}%")

        plt.plot(
            N,
            y_values,
            marker="o",
            linewidth=2,
            label=f"θ = {theta:.2f}",
        )

    plt.xscale("log")

    plt.xlabel("N cutoff (log scale)")
    plt.ylabel("Lopsided semiprimes (%)")
    plt.title("Lopsided Semiprime Convergence Across θ")

    plt.grid(True, alpha=0.3)
    plt.legend()

    plt.tight_layout()

    output_file = "multi_theta_convergence.png"
    plt.savefig(output_file, dpi=200)

    print(f"\nSaved: {output_file}")


if __name__ == "__main__":
    main()
