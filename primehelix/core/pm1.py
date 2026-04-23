"""
Pollard p-1 and Williams p+1 methods.
Both use gmpy2.powmod for fast modular exponentiation.
"""
from __future__ import annotations
import math

try:
    import gmpy2
    _powmod = lambda b, e, m: int(gmpy2.powmod(b, e, m))
    _gcd = lambda a, b: int(gmpy2.gcd(a, b))
except ImportError:
    _powmod = pow
    _gcd = math.gcd

from .primes import small_primes_up_to


def pollard_pm1(n: int, B1: int = 100_000) -> int | None:
    """
    Pollard p-1 stage 1. Finds a factor p where p-1 is B1-smooth.
    Returns a nontrivial factor or None.
    """
    if n % 2 == 0:
        return 2
    a = 2
    for p in small_primes_up_to(B1):
        # raise a to the highest power of p <= B1
        pk = p
        while pk * p <= B1:
            pk *= p
        a = _powmod(a, pk, n)

    g = _gcd(a - 1, n)
    if 1 < g < n:
        return g
    return None


def williams_pp1(n: int, B1: int = 100_000) -> int | None:
    """
    Williams p+1 stage 1 (Lucas sequence approach, two starting seeds).
    Finds a factor p where p+1 is B1-smooth.
    """
    if n % 2 == 0:
        return 2

    def _lucas_mul(V, k, P, n):
        """Multiply Lucas V sequence: V_{km} from V_m using doubling."""
        if k == 0:
            return 2
        if k == 1:
            return V % n
        bits = bin(k)[3:]
        Vk = V
        Vk1 = (V * V - 2) % n
        for bit in bits:
            if bit == '0':
                Vk1 = (Vk * Vk1 - V) % n
                Vk = (Vk * Vk - 2) % n
            else:
                Vk = (Vk * Vk1 - V) % n
                Vk1 = (Vk1 * Vk1 - 2) % n
        return Vk

    for seed in (2, 5):
        V = seed
        for p in small_primes_up_to(B1):
            pk = p
            while pk * p <= B1:
                pk *= p
            V = _lucas_mul(V, pk, seed, n)

        g = _gcd(V - 2, n)
        if 1 < g < n:
            return g

    return None
