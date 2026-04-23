"""
Tangent-split diagnostics. Two geometric interpretations of n:
  1. Equal split: L = 2*sqrt(n) — the "balanced" reference chord length
  2. Prime test: L = n+1 — discriminant check that reveals factor structure
     when n = p*q: discriminant = (p+q)^2 - 4pq = (p-q)^2, a perfect square.
Ported from onojk123/tangent_prime_test.py.
"""
from __future__ import annotations
import math
from dataclasses import dataclass

try:
    from decimal import Decimal, getcontext
    _HAS_DECIMAL = True
except ImportError:
    _HAS_DECIMAL = False


@dataclass
class EqualSplitInfo:
    n: int
    L: str          # 2*sqrt(n) as high-precision decimal string
    half: str       # sqrt(n)
    product: str    # half^2 ≈ n
    remainder: str  # product - n (should be near 0 for perfect square)


@dataclass
class TangentSplitInfo:
    n: int
    L: int              # n + 1
    discriminant: int   # L^2 - 4n = (n-1)^2
    sqrt_disc_exact: bool
    root1: int | None   # (L - sqrt(disc)) / 2 = 1
    root2: int | None   # (L + sqrt(disc)) / 2 = n


@dataclass
class IdealSplitInfo:
    n: int
    p: int
    q: int
    L: int              # p + q
    discriminant: int   # (p-q)^2  — always a perfect square
    roots: tuple[int, int]  # (p, q)


def equal_split(n: int, precision: int = 80) -> EqualSplitInfo:
    """Geometric reference: L = 2*sqrt(n)."""
    if _HAS_DECIMAL:
        getcontext().prec = precision
        nD = Decimal(n)
        half = nD.sqrt()
        L = 2 * half
        product = half * half
        remainder = product - nD
        return EqualSplitInfo(
            n=n, L=str(L), half=str(half),
            product=str(product), remainder=str(remainder),
        )
    sq = math.sqrt(n)
    return EqualSplitInfo(
        n=n, L=str(2 * sq), half=str(sq),
        product=str(sq * sq), remainder=str(sq * sq - n),
    )


def tangent_split(n: int) -> TangentSplitInfo:
    """
    Prime-test split: L = n+1, discriminant = L^2 - 4n = (n-1)^2.
    For a semiprime n=p*q, the "ideal" chord length is p+q,
    and its discriminant (p-q)^2 is a perfect square whose roots reveal p, q.
    """
    L = n + 1
    disc = L * L - 4 * n  # = (n-1)^2
    sqrt_disc = math.isqrt(disc)
    exact = (sqrt_disc * sqrt_disc == disc)
    r1 = (L - sqrt_disc) // 2 if exact else None
    r2 = (L + sqrt_disc) // 2 if exact else None
    return TangentSplitInfo(
        n=n, L=L, discriminant=disc,
        sqrt_disc_exact=exact, root1=r1, root2=r2,
    )


def ideal_split(n: int, p: int, q: int) -> IdealSplitInfo:
    """
    Given the actual factors p, q of n=p*q, compute the ideal chord.
    L = p+q; disc = (p-q)^2; roots = (p, q).
    This is the "target" that a geometric search would aim to find.
    """
    if p > q:
        p, q = q, p
    L = p + q
    disc = (q - p) ** 2
    return IdealSplitInfo(n=n, p=p, q=q, L=L, discriminant=disc, roots=(p, q))
