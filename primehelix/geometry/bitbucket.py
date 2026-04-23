"""
Bit-bucket analysis. Groups integers by bit-length, computes offsets,
prime density, and asymptotic comparisons.
Unified from geom_factor bit-bucket scripts.
"""
from __future__ import annotations
import math
from dataclasses import dataclass


@dataclass
class BitBucket:
    k: int          # bit-length (bucket index)
    lo: int         # 2^(k-1)
    hi: int         # 2^k - 1
    width: int      # 2^(k-1)
    n: int          # the number being analyzed
    offset: int     # n - 2^(k-1)
    offset_norm: float  # offset / width  (position in [0,1])


@dataclass
class BucketDensity:
    k: int
    prime_count: int
    bucket_width: int
    density: float          # primes / width
    asymptotic: float       # 1 / (k * ln2)
    ratio: float            # density / asymptotic


def bit_bucket(n: int) -> BitBucket:
    """Return the bit-bucket record for n."""
    if n < 1:
        raise ValueError("n must be >= 1")
    k = n.bit_length()
    lo = 1 << (k - 1)
    hi = (1 << k) - 1
    width = lo
    offset = n - lo
    return BitBucket(
        k=k, lo=lo, hi=hi, width=width,
        n=n, offset=offset,
        offset_norm=offset / width if width else 0.0,
    )


def bucket_density(k: int, prime_count: int) -> BucketDensity:
    """Compute prime density in bucket k, compare to PNT asymptotic."""
    width = 1 << (k - 1)
    density = prime_count / width
    asymptotic = 1.0 / (k * math.log(2))
    ratio = density / asymptotic if asymptotic else 0.0
    return BucketDensity(
        k=k,
        prime_count=prime_count,
        bucket_width=width,
        density=density,
        asymptotic=asymptotic,
        ratio=ratio,
    )


def prime_density_table(max_k: int = 20) -> list[BucketDensity]:
    """
    Build a prime density table for buckets k=2..max_k using a sieve.
    """
    limit = (1 << max_k) - 1
    sieve = bytearray([1]) * (limit + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(sieve[i*i::i]))

    rows = []
    for k in range(2, max_k + 1):
        lo = 1 << (k - 1)
        hi = (1 << k) - 1
        count = sum(sieve[lo:hi + 1])
        rows.append(bucket_density(k, count))
    return rows
