"""
Lenstra ECM (Elliptic Curve Method), Stage 1.
Pure Python + gmpy2. Finds factors where p has small group order over the curve.
"""
from __future__ import annotations
import math
import random
import time

try:
    import gmpy2
    _gcd = lambda a, b: int(gmpy2.gcd(a, b))
    _powmod = lambda b, e, m: int(gmpy2.powmod(b, e, m))
except ImportError:
    _gcd = math.gcd
    _powmod = pow

from .primes import small_primes_up_to


def _build_stage1_exponent(B1: int) -> int:
    M = 1
    for p in small_primes_up_to(B1):
        pk = p
        while pk * p <= B1:
            pk *= p
        M *= pk
    return M


def _inv_mod(a: int, n: int) -> int | None:
    g = _gcd(a % n, n)
    if g > 1:
        return None
    return pow(a % n, -1, n)


def _ec_add(P, Q, a: int, n: int):
    """Add two points on y^2 = x^3 + ax + b (mod n). Returns (point, factor_or_None)."""
    if P is None:
        return Q, None
    if Q is None:
        return P, None
    x1, y1 = P
    x2, y2 = Q

    if (x1 - x2) % n == 0:
        if (y1 + y2) % n == 0:
            return None, None  # point at infinity
        # doubling
        num = (3 * x1 * x1 + a) % n
        den = (2 * y1) % n
    else:
        num = (y2 - y1) % n
        den = (x2 - x1) % n

    inv = _inv_mod(den, n)
    if inv is None:
        g = _gcd(den % n, n)
        return None, g if 1 < g < n else None

    lam = (num * inv) % n
    x3 = (lam * lam - x1 - x2) % n
    y3 = (lam * (x1 - x3) - y1) % n
    return (x3, y3), None


def _ec_mul(k: int, P, a: int, n: int):
    """Scalar multiplication k*P. Returns (point, factor_or_None)."""
    R = None
    Q = P
    while k > 0:
        if k & 1:
            R, g = _ec_add(R, Q, a, n)
            if g:
                return None, g
        k >>= 1
        if k:
            Q, g = _ec_add(Q, Q, a, n)
            if g:
                return None, g
    return R, None


def ecm(n: int, B1: int = 50_000, curves: int = 200,
        seed: int | None = None, timeout_ms: int | None = None) -> int | None:
    """
    Lenstra ECM stage 1. Returns a nontrivial factor of n, or None.
    """
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3

    rng = random.Random(seed if seed is not None else 0xEC11)
    M = _build_stage1_exponent(B1)
    deadline = time.monotonic() + timeout_ms / 1000.0 if timeout_ms else None

    for _ in range(curves):
        if deadline and time.monotonic() > deadline:
            return None

        # Suyama parametrization: random curve via (x, y, a); b derived
        x = rng.randrange(2, n - 1)
        y = rng.randrange(2, n - 1)
        a = rng.randrange(2, n - 1)
        b = (y * y - x * x * x - a * x) % n

        disc = (-16 * (4 * pow(a, 3, n) + 27 * pow(b, 2, n))) % n
        g = _gcd(disc, n)
        if 1 < g < n:
            return g

        _, g = _ec_mul(M, (x, y), a, n)
        if g and 1 < g < n:
            return g

    return None
