"""
Quadratic Sieve with large-prime relations and GF(2) nullspace.
Ported from geom_factor with fixes: multiple nullspace vectors tried,
corrected congruence reconstruction, gmpy2-accelerated where possible.
"""
from __future__ import annotations
import math
from collections import defaultdict

try:
    import gmpy2
    _isqrt = lambda n: int(gmpy2.isqrt(n))
    _gcd = lambda a, b: int(gmpy2.gcd(a, b))
except ImportError:
    _isqrt = math.isqrt
    _gcd = math.gcd

from .primes import small_primes_up_to


def _tonelli_shanks(n: int, p: int) -> int | None:
    n %= p
    if n == 0:
        return 0
    if p == 2:
        return n % 2
    if pow(n, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = next(z for z in range(2, p) if pow(z, (p - 1) // 2, p) == p - 1)
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 1, pow(t, 2, p)
        while t2 != 1:
            t2 = pow(t2, 2, p)
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        r, t, c, m = (r * b) % p, (t * b * b) % p, (b * b) % p, i
    return r


def _build_factor_base(N: int, B: int):
    FB = []
    for p in small_primes_up_to(B):
        if p == 2:
            FB.append((2, 1, 1))
            continue
        if pow(N % p, (p - 1) // 2, p) != 1:
            continue
        r = _tonelli_shanks(N % p, p)
        if r is None:
            continue
        FB.append((p, r, (p - r) % p))
    return FB


def _trial_smooth(val: int, FB):
    fac = defaultdict(int)
    tmp = abs(val)
    for p, _, _ in FB:
        while tmp % p == 0:
            fac[p] += 1
            tmp //= p
    return fac, tmp


def _gf2_nullspace(mat):
    """
    GF(2) row reduction; returns all nullspace vectors found.
    Tries every free column, not just the first.
    """
    rows = [row[:] for row in mat]
    m, n = len(rows), len(rows[0]) if rows else 0
    pivcol = {}
    r = 0
    for c in range(n):
        pivot = next((rr for rr in range(r, m) if rows[rr][c]), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        for rr in range(m):
            if rr != r and rows[rr][c]:
                rows[rr] = [rows[rr][j] ^ rows[r][j] for j in range(n)]
        pivcol[r] = c
        r += 1

    free_cols = [c for c in range(n) if c not in pivcol.values()]
    vectors = []
    for fc in free_cols:
        v = [0] * n
        v[fc] = 1
        for row_idx in range(r - 1, -1, -1):
            pc = pivcol[row_idx]
            s = sum(rows[row_idx][c] * v[c] for c in range(pc + 1, n)) % 2
            v[pc] = s
        vectors.append(v)
    return vectors


def quadratic_sieve(N: int, B_scale: float = 5.0,
                    want_extra: int = 50, max_M: int = 500_000,
                    LPB: int = 10**12) -> int | None:
    """
    Quadratic Sieve. Returns a nontrivial factor of N or None.

    Parameters
    ----------
    N        : number to factor (composite, odd)
    B_scale  : multiply heuristic smoothness bound by this factor
    want_extra: relations to collect beyond len(FB)
    max_M    : sieve half-width
    LPB      : large prime bound for partial relations
    """
    if N % 2 == 0:
        return 2

    logN = math.log(N)
    loglogN = math.log(logN)
    B0 = math.exp(0.5 * math.sqrt(logN * loglogN))
    B = max(10_000, int(B0 * B_scale))
    M = min(int(2.5 * B), max_M)

    A = _isqrt(N)
    FB = _build_factor_base(N, B)
    if not FB:
        return None

    FB_index = {p: i for i, (p, _, _) in enumerate(FB)}
    target = len(FB) + want_extra

    # Sieve
    rels = []
    singles: dict[int, tuple] = {}
    ln_p = {p: math.log(p) for p, _, _ in FB}

    for x in range(-M, M + 1):
        Ax = A + x
        Qx = Ax * Ax - N
        if Qx == 0:
            continue

        score = math.log(abs(Qx))
        for p, r1, r2 in FB:
            for r in (r1, r2):
                if Ax % p == r:
                    score -= ln_p[p]
        if score > 0.75 * math.log(abs(Qx)):
            continue

        fac, leftover = _trial_smooth(Qx, FB)

        if leftover == 1:
            rels.append((Ax, fac, Qx < 0))
        elif leftover <= LPB:
            if leftover in singles:
                Ax2, fac2, neg2 = singles.pop(leftover)
                merged = fac.copy()
                for p, e in fac2.items():
                    merged[p] += e
                rels.append((Ax * Ax2 % N, merged, (Qx < 0) ^ neg2))
            else:
                singles[leftover] = (Ax, fac, Qx < 0)

        if len(rels) >= target:
            break

    if len(rels) < max(8, len(FB) // 4):
        return None

    # Build GF(2) matrix (include sign as extra column)
    all_primes = sorted(FB_index.keys())
    n_cols = len(all_primes) + 1  # +1 for sign
    mat = []
    for (Ax, fac, neg) in rels:
        row = [0] * n_cols
        row[-1] = int(neg)
        for p, e in fac.items():
            idx = FB_index.get(p)
            if idx is not None:
                row[idx] = e & 1
        mat.append(row)

    nullvecs = _gf2_nullspace(mat)
    if not nullvecs:
        return None

    for vec in nullvecs:
        X = 1
        exp_sum = defaultdict(int)
        neg_count = 0
        for take, (Ax, fac, neg) in zip(vec, rels):
            if not take:
                continue
            X = (X * Ax) % N
            for p, e in fac.items():
                exp_sum[p] += e
            if neg:
                neg_count += 1

        if neg_count % 2 != 0:
            continue

        Y = 1
        for p, e in exp_sum.items():
            Y = (Y * pow(p, e // 2, N)) % N

        for delta in (X - Y, X + Y):
            g = _gcd(delta % N, N)
            if 1 < g < N:
                return g

    return None
