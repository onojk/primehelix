"""
Primality testing: Baillie-PSW (Miller-Rabin base-2 + strong Lucas PRP).
Deterministic for all 64-bit integers; no known counterexamples beyond.
Uses gmpy2 for fast modular arithmetic on large integers.
"""
from __future__ import annotations
import math
try:
    import gmpy2
    _USE_GMPY2 = True
except ImportError:
    _USE_GMPY2 = False

_SMALL_PRIMES = [
    2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
    73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,
]

def _powmod(base: int, exp: int, mod: int) -> int:
    if _USE_GMPY2:
        return int(gmpy2.powmod(base, exp, mod))
    return pow(base, exp, mod)

def _isqrt(n: int) -> int:
    if _USE_GMPY2:
        return int(gmpy2.isqrt(n))
    return math.isqrt(n)

def _is_square(n: int) -> bool:
    if n < 0:
        return False
    r = _isqrt(n)
    return r * r == n

def _jacobi(a: int, n: int) -> int:
    """Jacobi symbol (a/n), n odd positive."""
    a %= n
    result = 1
    while a:
        t = a & -a
        v2 = (t.bit_length() - 1)
        if v2:
            if n % 8 in (3, 5):
                result = -result
            a >>= v2
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a, n = n % a, a
    return result if n == 1 else 0

def _miller_rabin(n: int, a: int) -> bool:
    d = n - 1
    s = (d & -d).bit_length() - 1
    d >>= s
    x = _powmod(a % n, d, n)
    if x == 1 or x == n - 1:
        return True
    for _ in range(s - 1):
        x = (x * x) % n
        if x == n - 1:
            return True
    return False

def _lucas_selfridge(n: int):
    D = 5
    while True:
        j = _jacobi(D, n)
        if j == -1:
            return D, 1, (1 - D) // 4
        D = -D - 2 if D > 0 else -D + 2

def _strong_lucas_prp(n: int) -> bool:
    if n < 2: return False
    if n % 2 == 0: return n == 2
    if _is_square(n): return False

    D, P, Q = _lucas_selfridge(n)
    d = n + 1
    s = (d & -d).bit_length() - 1
    d >>= s

    # inv(2) mod n — always exists since n is odd
    inv2 = (n + 1) >> 1

    # Start at (U_1, V_1, Q^1), then process remaining bits of d
    U, V, Qk = 1, P % n, Q % n
    for bit in bin(d)[3:]:   # skip leading '1b1'
        # Double: (U_{2k}, V_{2k}, Q^{2k})
        U, V, Qk = (U * V) % n, (V * V - 2 * Qk) % n, (Qk * Qk) % n
        # Add 1 if this bit is set: use correct step-add formulas with inv(2)
        if bit == '1':
            new_U = (P * U + V) * inv2 % n
            new_V = (D * U + P * V) * inv2 % n
            U, V = new_U, new_V
            Qk = Qk * Q % n

    # U == U_d, V == V_d
    if U == 0 or V == 0:
        return True
    # Check V_{d * 2^r} for r = 1 .. s-1
    for _ in range(s - 1):
        V = (V * V - 2 * Qk) % n
        Qk = (Qk * Qk) % n
        if V == 0:
            return True
    return False

def is_prime(n: int) -> bool:
    """Baillie-PSW probable prime test. Deterministic for n < 2^64."""
    if n < 2: return False
    for p in _SMALL_PRIMES:
        if n == p: return True
        if n % p == 0: return False
    if not _miller_rabin(n, 2):
        return False
    return _strong_lucas_prp(n)

def small_primes_up_to(limit: int) -> list[int]:
    """Sieve of Eratosthenes."""
    if limit < 2:
        return []
    sieve = bytearray([1]) * (limit + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(limit**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    return [i for i, v in enumerate(sieve) if v]

# Mod-210 wheel (coprime to 2,3,5,7)
_WHEEL_MOD = 210
_WHEEL_RESIDUES = tuple(r for r in range(_WHEEL_MOD) if all(r % p for p in (2, 3, 5, 7)))
