import json
import math


def _factor_dict(result):
    return {str(p): e for p, e in sorted(result.factors.items())}


def _prime_factors_flat(result):
    flat = []
    for p, e in sorted(result.factors.items()):
        flat.extend([p] * e)
    return flat


def _factorization_string(result, ascii_only=True):
    parts = []
    for p, e in sorted(result.factors.items()):
        if e == 1:
            parts.append(str(p))
        else:
            parts.append(f"{p}^{e}")
    sep = " * " if ascii_only else " × "
    return sep.join(parts)


def _coil_insight(coil):
    if coil is None:
        return None

    bit_gap = getattr(coil, "bit_gap", None)
    balance = getattr(coil, "balance", None)

    if bit_gap is None or balance is None:
        return None

    if bit_gap <= 1 and balance < 0.15:
        return "nearly equal factors; structure is strongly balanced"
    if bit_gap <= 3 and balance < 1.0:
        return "fairly balanced semiprime; neither factor strongly dominates"
    if bit_gap <= 8 and balance < 10.0:
        return "moderately lopsided semiprime; larger factor has visible structural weight"
    if bit_gap <= 16 and balance < 100.0:
        return "lopsided semiprime; structure is dominated by the larger factor"
    return "extremely lopsided semiprime; small prime is nearly invisible structurally"


def _serialize_coil(coil):
    if coil is None:
        return None

    data = {}
    for key, value in vars(coil).items():
        if isinstance(value, (int, float, str, bool)) or value is None:
            data[key] = value
        else:
            data[key] = str(value)

    insight = _coil_insight(coil)
    if insight is not None:
        data["insight"] = insight

    return data


def structure_summary(classification, coil=None, residue=None):
    """
    Public compact structural identity string for scripting and aggregation.

    Label format (parts joined by " | "):
      Part 1 — classification: "prime", "semiprime", "composite", "invalid"
      Part 2 — balance (semiprimes with 2 distinct primes only):
                "balanced"  bit_gap ≤ 1 and balance < 0.15
                "moderate"  bit_gap ≤ 8 and balance < 10.0
                "lopsided"  otherwise
                omitted for prime squares (p²) — only one distinct prime
      Part 3 — residue family:
                primes:      "pythagorean" (p≡1 mod4), "gaussian" (p≡3 mod4)
                             omitted for p=2
                semiprimes:  "mod4_1x1", "mod4_1x3", "mod4_3x3",
                             "mod4_2x1", "mod4_2x3", "mod4_2x2" (factor of 2)

    Examples:
      prime | pythagorean
      prime | gaussian
      prime                           (p=2)
      semiprime | balanced | mod4_1x1
      semiprime | lopsided | mod4_1x3
      semiprime | mod4_3x3            (p² form — no balance tier)
      composite
      invalid
    """
    parts = []

    if classification:
        parts.append(classification.lower())

    if coil is not None:
        bit_gap = getattr(coil, "bit_gap", None)
        balance = getattr(coil, "balance", None)

        if bit_gap is not None and balance is not None:
            if bit_gap <= 1 and balance < 0.15:
                parts.append("balanced")
            elif bit_gap <= 8 and balance < 10.0:
                parts.append("moderate")
            else:
                parts.append("lopsided")

    if residue:
        pair = residue.get("semiprime_mod4_pair")
        if pair:
            parts.append(f"mod4_{pair}")

        family = residue.get("prime_family")
        if family:
            parts.append(family)

    return " | ".join(parts) if parts else None


def build_json_result(
    result,
    command: str,
    classification: str | None = None,
    coil=None,
    residue=None,
):
    payload = {
        "command": command,
        "n": result.n,
        "classification": classification.lower() if classification else None,
        "factors": _factor_dict(result),
        "prime_factors": _prime_factors_flat(result),
        "factorization": _factorization_string(result, ascii_only=True),
        "method": result.method,
        "elapsed_ms": round(result.elapsed_ms, 3),
        "complete": result.complete,
        "steps": list(getattr(result, "steps", []) or []),
    }

    coil_data = _serialize_coil(coil)
    if coil_data is not None:
        payload["coil"] = coil_data

    if residue is not None:
        payload["residue"] = residue

    structure = structure_summary(classification, coil=coil, residue=residue)
    if structure is not None:
        payload["structure"] = structure

    return payload


def label_entropy(counts: dict, total: int) -> float | None:
    """Shannon entropy (bits) of a structure-label distribution.

    Returns None if total is 0. A single-label distribution has entropy 0;
    a uniform distribution over k labels has entropy log2(k).
    """
    if not total:
        return None
    return round(
        -sum((c / total) * math.log2(c / total) for c in counts.values() if c > 0),
        4,
    )


def print_json(payload):
    print(json.dumps(payload, indent=2, sort_keys=False))
