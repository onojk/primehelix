"""
Conical helix model. Maps integers to 3D points on a conical helix and
computes geometric footprints for semiprimes.
Ported and unified from rsacrack/coil_classifier.py.
"""
from __future__ import annotations
import math
import hashlib
import json
from dataclasses import dataclass


@dataclass
class CoilPoint:
    n: int
    x: float
    y: float
    z: float

    def dist_to(self, other: "CoilPoint") -> float:
        return math.sqrt(
            (self.x - other.x) ** 2 +
            (self.y - other.y) ** 2 +
            (self.z - other.z) ** 2
        )


@dataclass
class CoilFootprint:
    n: int
    p: int
    q: int
    d_n_to_q: float
    d_q_to_p: float
    d_p_to_1: float
    f1: float  # normalized d1
    f2: float  # normalized d2
    f3: float  # normalized d3
    balance: float      # |p-q| / sqrt(n)
    bit_gap: int        # |bitlen(p) - bitlen(q)|
    s1: float           # d1 per integer step
    s2: float           # d2 per integer step
    s3: float           # d3 per integer step
    sig_geom: str       # SHA-256 including coil params
    sig_invariant: str  # SHA-256 of factor pair only


def coil_point(n: int, r0: float = 1.0, alpha: float = 0.0125,
               beta: float = 0.005, L: float = 360.0) -> CoilPoint:
    """Map integer n to a point on the conical helix."""
    r = r0 + alpha * n
    theta = (2.0 * math.pi / L) * n
    return CoilPoint(
        n=n,
        x=r * math.cos(theta),
        y=r * math.sin(theta),
        z=beta * n,
    )


def coil_footprint(n: int, p: int, q: int,
                   r0: float = 1.0, alpha: float = 0.0125,
                   beta: float = 0.005, L: float = 360.0) -> CoilFootprint:
    """
    Compute the geometric footprint of semiprime n = p*q (p <= q)
    on the conical helix with given parameters.
    """
    if p > q:
        p, q = q, p

    Pn = coil_point(n, r0, alpha, beta, L)
    Pq = coil_point(q, r0, alpha, beta, L)
    Pp = coil_point(p, r0, alpha, beta, L)
    P1 = coil_point(1, r0, alpha, beta, L)

    d1 = Pn.dist_to(Pq)
    d2 = Pq.dist_to(Pp)
    d3 = Pp.dist_to(P1)

    total = d1 + d2 + d3 or 1.0
    f1, f2, f3 = d1 / total, d2 / total, d3 / total

    s1 = d1 / max(1, n - q)
    s2 = d2 / max(1, q - p)
    s3 = d3 / max(1, p - 1)

    balance = abs(q - p) / math.sqrt(n)
    bit_gap = abs(p.bit_length() - q.bit_length())

    def _sha(payload: dict) -> str:
        data = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
        return hashlib.sha256(data).hexdigest()

    sig_geom = _sha({
        "n": n, "primes": [p, q],
        "geom": {"r0": r0, "alpha": alpha, "beta": beta, "L": L},
        "norm": [round(f1, 6), round(f2, 6), round(f3, 6)],
    })
    sig_invariant = _sha({
        "n": n, "primes": [p, q],
        "bit_gap": bit_gap,
        "balance": round(balance, 8),
        "log_ratio": round(math.log(q) - math.log(p), 8),
    })

    return CoilFootprint(
        n=n, p=p, q=q,
        d_n_to_q=d1, d_q_to_p=d2, d_p_to_1=d3,
        f1=f1, f2=f2, f3=f3,
        balance=balance, bit_gap=bit_gap,
        s1=s1, s2=s2, s3=s3,
        sig_geom=sig_geom, sig_invariant=sig_invariant,
    )
