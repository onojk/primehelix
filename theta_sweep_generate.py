from math import isqrt
from lopsided import summarize_lopsided


THETAS = [0.20, 0.25, 0.30, 0.35, 0.40]


def sieve(limit):
    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"

    for i in range(2, isqrt(limit) + 1):
        if is_prime[i]:
            start = i * i
            step = i
            is_prime[start:limit + 1:step] = b"\x00" * (((limit - start) // step) + 1)

    return [i for i in range(2, limit + 1) if is_prime[i]]


def generate_semiprimes(limit):
    primes = sieve(limit // 2)
    rows = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            n = p * q
            if n > limit:
                break
            rows.append((n, p, q))

    return rows


def run(limit):
    semiprimes = generate_semiprimes(limit)

    print(f"\n=== THETA SWEEP GENERATED: N <= {limit:,} ===")
    print(f"{'theta':>8} {'total':>10} {'lopsided':>12} {'balanced':>12} {'ratio':>10}")
    print("-" * 60)

    for theta in THETAS:
        s = summarize_lopsided(semiprimes, theta=theta)
        print(
            f"{s['theta']:>8.2f} "
            f"{s['total']:>10} "
            f"{s['lopsided']:>12} "
            f"{s['balanced']:>12} "
            f"{s['lopsided_ratio']:>10.4f}"
        )


if __name__ == "__main__":
    for limit in [10_000, 100_000, 1_000_000]:
        run(limit)
