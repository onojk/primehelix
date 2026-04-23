"""
Pollard Rho (Brent variant) with batch-GCD acceleration.
Uses gmpy2.gcd when available for speed on large integers.
"""
from __future__ import annotations
import math
import random
import time

try:
    import gmpy2
    _gcd = lambda a, b: int(gmpy2.gcd(a, b))
except ImportError:
    _gcd = math.gcd


def pollard_rho(n: int, timeout_ms: int = 3000, seed: int | None = None) -> int | None:
    """
    Return a nontrivial factor of n, or None on timeout/failure.
    Uses Brent's improvement with batch product-of-differences GCD.
    """
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3

    rng = random.Random(seed if seed is not None else random.randrange(2**32))
    deadline = time.monotonic() + timeout_ms / 1000.0

    while time.monotonic() < deadline:
        y = rng.randrange(1, n - 1)
        c = rng.randrange(1, n - 1)
        batch = 128
        g = r = q = 1

        f = lambda x: (x * x + c) % n

        while time.monotonic() < deadline and g == 1:
            x = y
            for _ in range(r):
                y = f(y)
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(batch, r - k)):
                    y = f(y)
                    diff = x - y if x >= y else y - x
                    q = (q * (diff % n)) % n
                g = _gcd(q, n)
                k += batch
            r <<= 1

        if g == n:
            while True:
                ys = f(ys)
                g = _gcd(abs(x - ys), n)
                if g > 1:
                    break

        if 1 < g < n:
            return g

    return None
