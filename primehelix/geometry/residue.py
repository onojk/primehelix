from __future__ import annotations


def _prime_family_mod4(p: int) -> str | None:
    """
    For odd primes:
      1 mod 4 -> Pythagorean prime
      3 mod 4 -> Gaussian prime

    Returns None for non-odd-prime-style cases like 2.
    """
    if p == 2:
        return None
    r = p % 4
    if r == 1:
        return "pythagorean"
    if r == 3:
        return "gaussian"
    return None


def _flatten_factors(factors: dict[int, int]) -> list[int]:
    flat = []
    for p, e in sorted(factors.items()):
        flat.extend([p] * e)
    return flat


def residue_profile(n: int, factors: dict[int, int] | None = None, classification: str | None = None) -> dict:
    """
    Residue/arithmetic-family profile for a number and its factors.

    Returns a JSON-friendly dict.
    """
    data: dict = {
        "n_mod_4": n % 4,
        "n_mod_6": n % 6,
        "n_mod_30": n % 30,
        "classification": classification.lower() if classification else None,
    }

    if not factors:
        return data

    flat = _flatten_factors(factors)

    data["factor_mod_4"] = [p % 4 for p in flat]
    data["factor_mod_6"] = [p % 6 for p in flat]
    data["factor_mod_30"] = [p % 30 for p in flat]
    data["factor_families_mod4"] = [_prime_family_mod4(p) for p in flat]

    if classification and classification.lower() == "prime" and len(flat) == 1:
        data["prime_family"] = _prime_family_mod4(flat[0])

    if classification and classification.lower() == "semiprime" and len(flat) == 2:
        a, b = flat

        # Handle factor 2 explicitly
        if a == 2 or b == 2:
            other = b if a == 2 else a
            data["semiprime_mod4_pair"] = f"2x{other % 4}"
            data["semiprime_mod4_note"] = "includes factor 2; mod4 structure is degenerate"
        else:
            ra = a % 4
            rb = b % 4
            pair = f"{min(ra, rb)}x{max(ra, rb)}"
            data["semiprime_mod4_pair"] = pair

            if pair == "1x1":
                data["semiprime_mod4_note"] = "both factors are 1 mod 4"
            elif pair == "1x3":
                data["semiprime_mod4_note"] = "mixed 1 mod 4 and 3 mod 4 factor families"
            elif pair == "3x3":
                data["semiprime_mod4_note"] = "both factors are 3 mod 4"

    return data
