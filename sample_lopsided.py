import random
from math import isqrt

from lopsided import is_lopsided


PRIME_CAP = 200_000
THETA = 0.25
ACCEPTED_TARGET = 3_000


def sieve(limit):
    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"

    for i in range(2, isqrt(limit) + 1):
        if is_prime[i]:
            is_prime[i * i : limit + 1 : i] = b"\x00" * (
                ((limit - i * i) // i) + 1
            )

    return [i for i in range(2, limit + 1) if is_prime[i]]


primes = sieve(PRIME_CAP)


def sample_semiprimes(limit, accepted_target=ACCEPTED_TARGET, theta=THETA):
    accepted = 0
    attempts = 0
    lopsided = 0

    while accepted < accepted_target:
        p = random.choice(primes)
        q = random.choice(primes)

        if p > q:
            p, q = q, p

        n = p * q
        attempts += 1

        if n > limit:
            continue

        accepted += 1

        if is_lopsided(n, p, q, theta=theta):
            lopsided += 1

    return {
        "limit": limit,
        "accepted": accepted,
        "attempts": attempts,
        "lopsided": lopsided,
        "ratio": lopsided / accepted,
    }


if __name__ == "__main__":
    print(f"\n=== SAMPLED LOPSIDED SEMIPRIMES ===")
    print(f"theta = {THETA}")
    print(f"prime cap = {PRIME_CAP:,}")
    print(f"accepted target = {ACCEPTED_TARGET:,}\n")

    for N in [1_000_000, 10_000_000, 100_000_000]:
        s = sample_semiprimes(N)

        print(
            f"N <= {s['limit']:,} | "
            f"accepted={s['accepted']:,} | "
            f"attempts={s['attempts']:,} | "
            f"lopsided={s['lopsided']:,} | "
            f"{100 * s['ratio']:.2f}%"
        )
