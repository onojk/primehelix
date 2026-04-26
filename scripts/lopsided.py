"""
Lopsided semiprime classification

Definition:
For a semiprime n = p * q with p <= q,
n is θ-lopsided if:

    min(p, q) <= n**θ

Baseline:
    θ = 0.25
"""

from math import pow


def is_lopsided(n: int, p: int, q: int, theta: float = 0.25) -> bool:
    """
    Return True if semiprime n = p*q is θ-lopsided.
    """
    return min(p, q) <= pow(n, theta)


def classify_semiprime(n: int, p: int, q: int, theta: float = 0.25) -> str:
    """
    Classify semiprime as 'lopsided' or 'balanced'
    """
    if is_lopsided(n, p, q, theta):
        return "lopsided"
    return "balanced"


def lopsided_ratio(semiprimes: list, theta: float = 0.25) -> float:
    """
    Compute fraction of lopsided semiprimes.

    semiprimes: list of tuples [(n, p, q), ...]
    """
    if not semiprimes:
        return 0.0

    count = 0

    for n, p, q in semiprimes:
        if is_lopsided(n, p, q, theta):
            count += 1

    return count / len(semiprimes)


def summarize_lopsided(semiprimes: list, theta: float = 0.25) -> dict:
    """
    Return summary stats.
    """
    total = len(semiprimes)
    lopsided_count = 0

    for n, p, q in semiprimes:
        if is_lopsided(n, p, q, theta):
            lopsided_count += 1

    balanced_count = total - lopsided_count

    return {
        "theta": theta,
        "total": total,
        "lopsided": lopsided_count,
        "balanced": balanced_count,
        "lopsided_ratio": lopsided_count / total if total else 0.0,
    }


if __name__ == "__main__":
    # quick sanity test
    test_data = [
        (91, 7, 13),     # not lopsided
        (202, 2, 101),   # lopsided
        (15, 3, 5),      # borderline small
        (77, 7, 11),     # not lopsided
        (26, 2, 13),     # lopsided
    ]

    summary = summarize_lopsided(test_data, theta=0.25)

    print("\n=== LOPSIDED TEST ===")
    for k, v in summary.items():
        print(f"{k}: {v}")
