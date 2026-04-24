from bisect import bisect_left, bisect_right
from math import ceil, isqrt
import os


THETA = 0.25

# FAST_MODE for CI (smaller limits)
FAST_MODE = os.getenv("FAST_MODE") == "1"

LIMITS = (
    [10_000, 100_000, 1_000_000]
    if FAST_MODE
    else [10_000, 100_000, 1_000_000, 10_000_000, 100_000_000]
)


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


def count_semiprimes(limit, primes, theta):
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

    balanced = total - lopsided
    ratio = lopsided / total if total else 0.0

    return total, lopsided, balanced, ratio


def main():
    max_limit = max(LIMITS)
    prime_limit = max_limit // 2

    mode = "FAST_MODE" if FAST_MODE else "FULL_MODE"
    print(f"Sieve up to {prime_limit:,}... ({mode})")

    primes = sieve(prime_limit)
    print(f"Primes found: {len(primes):,}")

    print("\n=== EXACT LOPSIDED COUNTS ===")
    print(f"Definition: n = p*q is lopsided if min(p,q) <= n^{THETA}")
    print()
    print(f"{'N':>14} {'total':>14} {'lopsided':>14} {'balanced':>14} {'percent':>10}")
    print("-" * 72)

    for limit in LIMITS:
        total, lop, bal, ratio = count_semiprimes(limit, primes, THETA)
        print(
            f"{limit:>14,} "
            f"{total:>14,} "
            f"{lop:>14,} "
            f"{bal:>14,} "
            f"{100 * ratio:>9.2f}%"
        )


if __name__ == "__main__":
    main()
