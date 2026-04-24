"""
Theoretical semiprime structure predictions via prime-sum model.

For a semiprime n = p·q (p ≤ q), the 'lopsided' tier is triggered when
balance ≥ 10 OR bit_gap > 8. The balance condition (balance = (q-p)/√n ≥ 10)
is the dominant threshold: solving (r-1)/√r = 10 gives r = q/p ≥ 101.99.

So: lopsided ↔ p < √(N / 101.99) for semiprimes in [1, N].

The theoretical fraction is estimated by summing π(N/p) - π(p) over qualifying
primes using the offset logarithmic integral li(x) as the prime-counting
approximation. No external dependencies required.
"""
from __future__ import annotations

import math

_EULER_MASCHERONI = 0.5772156649015328606

# q/p ratio threshold at which balance = 10: (r-1)/sqrt(r) = 10 → r ≈ 101.99
_LOPSIDED_RATIO = ((10 + math.sqrt(104)) / 2) ** 2


def _small_sieve(limit: int) -> list[int]:
    if limit < 2:
        return []
    sieve = bytearray([1]) * (limit + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, math.isqrt(limit) + 1):
        if sieve[i]:
            sieve[i * i :: i] = bytearray(len(sieve[i * i :: i]))
    return [i for i in range(2, limit + 1) if sieve[i]]


def _li(x: float) -> float:
    """Logarithmic integral li(x) = Ei(ln x) via convergent series.

    Series: Ei(t) = γ + ln(t) + Σ_{k=1}^∞ t^k / (k · k!)
    Recurrence: a_{k+1} = a_k · k·t / (k+1)²
    Converges for all t > 0; ~56 terms for t = ln(5×10^5).
    """
    if x <= 1.0:
        return 0.0
    t = math.log(x)
    a = t  # a_1 = t^1 / (1 · 1!)
    total = _EULER_MASCHERONI + math.log(t) + a
    for k in range(1, 350):
        a = a * k * t / (k + 1) ** 2
        total += a
        if a < 1e-14:
            break
    return total


_LI2 = _li(2.0)


def _pi_approx(x: float) -> float:
    """π(x) ≈ li(x) − li(2). Returns 0 for x < 2."""
    if x < 2.0:
        return 0.0
    return max(0.0, _li(x) - _LI2)


def theoretical_lopsided_fraction(N: int) -> float:
    """
    Theoretical fraction of semiprimes n = p·q ≤ N that are lopsided,
    derived from the prime distribution (π(x) ≈ li(x)).

    A semiprime p·q (p ≤ q) is lopsided when q/p > 101.99 — the ratio at
    which balance = (q−p)/√(pq) crosses 10.0.

    For each prime p ≤ √N:
      total pairs:    q prime in (p,    N/p]  →  π(N/p) − π(p)
      lopsided pairs: q prime in (102p, N/p]  →  π(N/p) − π(102p)

    Parameters
    ----------
    N : upper bound of the integer range (exclusive)

    Returns
    -------
    float in [0, 1]
    """
    sqrt_N = math.isqrt(N)
    primes = _small_sieve(sqrt_N + 1)

    total = 0.0
    lopsided = 0.0

    for p in primes:
        if p * p > N:
            break
        q_max = N / p
        total += max(0.0, _pi_approx(q_max) - _pi_approx(p))
        q_lopsided_min = _LOPSIDED_RATIO * p
        if q_lopsided_min < q_max:
            lopsided += max(0.0, _pi_approx(q_max) - _pi_approx(q_lopsided_min))

    return lopsided / total if total > 0 else 0.0


def compare_dickman_table(
    N_values: list[int],
    measured: dict[int, float],
) -> list[dict]:
    """
    Build the Dickman comparison table: theoretical vs measured lopsided %.

    Parameters
    ----------
    N_values : list of range upper bounds
    measured : dict mapping N → measured lopsided fraction (0.0–1.0)

    Returns
    -------
    List of dicts with keys:
        N, label, theoretical_pct, measured_pct, delta_pp
    """
    rows = []
    for N in N_values:
        theory = theoretical_lopsided_fraction(N)
        meas = measured.get(N)
        exp = f"1e{int(round(math.log10(N)))}"
        rows.append({
            "N": N,
            "label": exp,
            "theoretical_pct": round(theory * 100, 2),
            "measured_pct": round(meas * 100, 2) if meas is not None else None,
            "delta_pp": round((meas - theory) * 100, 2) if meas is not None else None,
        })
    return rows
