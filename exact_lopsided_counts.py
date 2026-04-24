from bisect import bisect_left, bisect_right
from math import isqrt


THETA = 0.25
LIMITS = [10_000, 100_000, 1_000_000, 10_000_000, 100_000_000]


def sieve(limit):
    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"

    for i in range(2, isqrt(limit) + 1):
        if is_prime[i]:
            is_prime[i * i : limit + 1 : i] = b"\x00" * (
                ((limit - i * i) // i) + 1
            )

    return [i for i in range(2, limit + 1) if is_prime[i]]


def count_semiprimes(limit, primes):
    total = 0
    lopsided = 0

    # p <= q and p*q <= limit
    for p in primes:
        if p * p > limit:
            break

        q_min_total = p
        q_max = limit // p

        total += bisect_right(primes, q_max) - bisect_left(primes, q_min_total)

        # For theta = 1/4:
        # p <= (p*q)^1/4
        # p^4 <= p*q
        # p^3 <= q
        q_min_lopsided = max(p, p ** 3)

        if q_min_lopsided <= q_max:
            lopsided += bisect_right(primes, q_max) - bisect_left(primes, q_min_lopsided)

    balanced = total - lopsided
    ratio = lopsided / total if total else 0

    return total, lopsided, balanced, ratio


def main():
    max_limit = max(LIMITS)
    prime_limit = max_limit // 2

    print(f"Sieve up to {prime_limit:,}...")
    primes = sieve(prime_limit)
    print(f"Primes found: {len(primes):,}")

    print("\n=== EXACT LOPSIDED COUNTS ===")
    print("Definition: n = p*q is lopsided if min(p,q) <= n^0.25")
    print()
    print(f"{'N':>14} {'total':>14} {'lopsided':>14} {'balanced':>14} {'percent':>10}")
    print("-" * 72)

    for limit in LIMITS:
        total, lop, bal, ratio = count_semiprimes(limit, primes)
        print(
            f"{limit:>14,} "
            f"{total:>14,} "
            f"{lop:>14,} "
            f"{bal:>14,} "
            f"{100 * ratio:>9.2f}%"
        )


if __name__ == "__main__":
    main()
