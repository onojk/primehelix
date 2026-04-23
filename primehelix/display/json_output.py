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
    """
    Plain-English interpretation of coil geometry for JSON output.
    Uses bit_gap and balance heuristics only, so it stays stable and simple.
    """
    if coil is None:
        return None

    bit_gap = getattr(coil, "bit_gap", None)
    balance = getattr(coil, "balance", None)

    if bit_gap is None or balance is None:
        return None

    # Keep this simple and stable.
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


def build_json_result(result, command: str, classification: str | None = None, coil=None):
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

    return payload


def print_json(payload):
    print(json.dumps(payload, indent=2, sort_keys=False))
