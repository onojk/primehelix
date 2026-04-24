import random
from bisect import bisect_left
from math import isqrt

from sympy import isprime, primepi


THETA = 0.25
SAMPLES = 20_000
LIMITS = [100_000_000, 1_000_000_000]


def sieve(limit):
    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"

    for i in range(2, isqrt(limit) + 1):
        if is_prime[i]:
            is_prime[i * i : limit + 1 : i] = b"\x00" * (
                ((limit - i * i) // i) + 1
            )

    return [i for i in range(2, limit + 1) if is_prime[i]]


def random_prime_between(lo, hi):
    """
    Return a random prime q with lo <= q <= hi using rejection sampling.
    """
    while True:
        q = random.randint(lo, hi)
        if q < 2:
            continue
        if isprime(q):
            return q


def build_weighted_p_table(limit):
    """
    Build weighted table for choosing p.

    Semiprime pair condition:
        p <= q
        p*q <= limit

    For each prime p <= sqrt(limit), valid q are primes:
        p <= q <= limit // p

    Weight of p = number of valid q values.
    """
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


def sample_one_semiprime_pair(cumulative, total_weight):
    """
    Uniformly sample one semiprime factor pair (p, q)
    from all prime pairs p <= q and p*q <= N.
    """
    r = random.randint(1, total_weight)
    idx = bisect_left([row[0] for row in cumulative], r)

    _, p, q_min, q_max = cumulative[idx]
    q = random_prime_between(q_min, q_max)

    return p, q


def is_lopsided_pair_theta_quarter(p, q):
    """
    For n = p*q and p <= q:

        p <= (p*q)^0.25

    is equivalent to:

        p^3 <= q

    This avoids floating-point error.
    """
    return p ** 3 <= q


def estimate(limit, samples=SAMPLES):
    cumulative, total_weight = build_weighted_p_table(limit)

    lopsided = 0

    for i in range(samples):
        p, q = sample_one_semiprime_pair(cumulative, total_weight)

        if is_lopsided_pair_theta_quarter(p, q):
            lopsided += 1

        if (i + 1) % 5000 == 0:
            print(f"  progress: {i + 1:,}/{samples:,}")

    ratio = lopsided / samples

    return {
        "limit": limit,
        "samples": samples,
        "total_pairs": total_weight,
        "lopsided": lopsided,
        "ratio": ratio,
    }


def main():
    print("\n=== HIGH-N UNIFORM SEMIPRIME-PAIR SAMPLER ===")
    print("Definition: n = p*q, p <= q, lopsided iff p^3 <= q")
    print(f"theta = {THETA}")
    print(f"samples per N = {SAMPLES:,}\n")

    for limit in LIMITS:
        print(f"Sampling N <= {limit:,}...")
        result = estimate(limit)

        print(
            f"N <= {result['limit']:,} | "
            f"estimated total pairs={result['total_pairs']:,} | "
            f"samples={result['samples']:,} | "
            f"lopsided={result['lopsided']:,} | "
            f"{100 * result['ratio']:.2f}%"
        )
        print()


if __name__ == "__main__":
    main()
