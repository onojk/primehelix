"""
Unified factoring pipeline.
Trial → p-1 → p+1 → Rho → ECM → QS
Each stage is tried with a time budget; escalates only if needed.
"""
from __future__ import annotations
import time
from dataclasses import dataclass, field

from .primes import is_prime, _SMALL_PRIMES
from .rho import pollard_rho
from .pm1 import pollard_pm1, williams_pp1
from .ecm import ecm
from .qs import quadratic_sieve


@dataclass
class FactorResult:
    n: int
    factors: dict[int, int]   # prime -> exponent
    method: str
    steps: list[str] = field(default_factory=list)
    elapsed_ms: float = 0.0
    complete: bool = True     # False if timed out with remaining composite cofactor


def _small_trial(n: int) -> int | None:
    for p in _SMALL_PRIMES:
        if n % p == 0:
            return p
    return None


def _factor_one(n: int, budget_ms: int, steps: list[str]) -> int | None:
    """Find one nontrivial factor of n using the full pipeline."""
    start = time.monotonic()
    remaining = lambda: max(50, budget_ms - int((time.monotonic() - start) * 1000))

    # Trial division (small primes)
    f = _small_trial(n)
    if f:
        steps.append(f"trial: {f}")
        return f

    digits = len(str(n))

    # p-1
    B1_pm1 = min(500_000, 200 * digits)
    f = pollard_pm1(n, B1=B1_pm1)
    if f:
        steps.append(f"p-1 (B1={B1_pm1}): {f}")
        return f

    # p+1
    f = williams_pp1(n, B1=B1_pm1)
    if f:
        steps.append(f"p+1 (B1={B1_pm1}): {f}")
        return f

    # Rho
    rho_budget = min(2000, remaining() // 2)
    f = pollard_rho(n, timeout_ms=rho_budget)
    if f:
        steps.append(f"rho: {f}")
        return f

    # ECM
    ecm_curves = 30 if digits < 40 else (60 if digits < 60 else 120)
    ecm_budget = min(5000, remaining() // 2)
    f = ecm(n, B1=50_000, curves=ecm_curves, timeout_ms=ecm_budget)
    if f:
        steps.append(f"ecm: {f}")
        return f

    # QS (last resort — for hard semiprimes)
    if digits <= 80 and remaining() > 500:
        f = quadratic_sieve(n)
        if f:
            steps.append(f"qs: {f}")
            return f

    return None


def factor(n: int, budget_ms: int = 10_000) -> FactorResult:
    """
    Fully factor n. Returns a FactorResult with prime factorization.
    """
    start = time.monotonic()
    steps: list[str] = []
    factors: dict[int, int] = {}

    def _record(p: int):
        factors[p] = factors.get(p, 0) + 1

    def _recurse(n: int, depth: int = 0) -> bool:
        if n == 1:
            return True
        if is_prime(n):
            _record(n)
            return True
        if depth > 40:
            _record(n)  # give up, treat as prime
            return False

        elapsed = int((time.monotonic() - start) * 1000)
        remaining = budget_ms - elapsed
        if remaining < 50:
            if is_prime(n):
                _record(n)
                return True
            _record(n)
            return False

        f = _factor_one(n, remaining, steps)
        if f is None:
            if is_prime(n):
                _record(n)
                return True
            steps.append(f"gave up on {n} ({len(str(n))} digits)")
            _record(n)
            return False

        ok1 = _recurse(f, depth + 1)
        ok2 = _recurse(n // f, depth + 1)
        return ok1 and ok2

    complete = _recurse(n)
    elapsed_ms = (time.monotonic() - start) * 1000

    method = "trial"
    for step in reversed(steps):
        for m in ("qs", "ecm", "rho", "p+1", "p-1", "trial"):
            if m in step:
                method = m
                break
        else:
            continue
        break

    return FactorResult(
        n=n,
        factors=factors,
        method=method,
        steps=steps,
        elapsed_ms=elapsed_ms,
        complete=complete,
    )


def classify(n: int, budget_ms: int = 10_000) -> tuple[str, FactorResult]:
    """
    Classify n as 'prime', 'semiprime', or 'composite'.
    Returns (classification, FactorResult).
    """
    if n < 2:
        dummy = FactorResult(n=n, factors={}, method="none", complete=True)
        return "invalid", dummy

    if is_prime(n):
        r = FactorResult(n=n, factors={n: 1}, method="bpsw", complete=True)
        return "prime", r

    r = factor(n, budget_ms=budget_ms)
    total_prime_factors = sum(r.factors.values())
    if total_prime_factors == 2 and r.complete:
        return "semiprime", r
    return "composite", r
