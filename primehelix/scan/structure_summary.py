from __future__ import annotations

from collections import Counter

from primehelix.display.json_output import structure_summary
from primehelix.geometry.residue import residue_profile
from primehelix.geometry.coil import coil_footprint


def _coil_for_result(n: int, classification: str, factors: dict[int, int]):
    if classification.lower() != "semiprime":
        return None

    primes = sorted(factors.keys())
    if len(primes) == 2:
        return coil_footprint(n, primes[0], primes[1])
    if len(primes) == 1:
        p = primes[0]
        return coil_footprint(n, p, p)
    return None


def summarize_range(start: int, stop: int, classify_fn, budget_ms: int = 2000):
    """
    Aggregate structure labels over [start, stop).

    Returns a dict:
      {
        "total": int,
        "counts": Counter,
      }
    """
    counts = Counter()
    total = 0

    for n in range(start, stop):
        classification, result = classify_fn(n, budget_ms=budget_ms)
        residue = residue_profile(n, result.factors, classification=classification)
        coil = _coil_for_result(n, classification, result.factors)
        label = structure_summary(classification, coil=coil, residue=residue)

        if label is None:
            label = classification.lower()

        counts[label] += 1
        total += 1

    return {
        "total": total,
        "counts": counts,
    }
