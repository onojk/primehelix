"""
Wheel-accelerated range scanner. Uses a mod-210 wheel (coprime to 2,3,5,7)
to skip ~77% of candidates. Resumable, supports gzip CSV output.
Unified from onojk123/coil_adaptive_wheel.py.
"""
from __future__ import annotations
import csv
import gzip
import io
import os
import sys
from bisect import bisect_right
from collections import Counter
from typing import Iterator

from ..core.primes import is_prime
from ..core.factor import factor as do_factor

# Mod-210 wheel
_MOD = 210
_WHEEL = tuple(sorted(r for r in range(_MOD) if all(r % p for p in (2, 3, 5, 7))))
_OFFSETS: list[int] = []
_prev = 0
for _r in _WHEEL:
    if _r == 0:
        continue
    _OFFSETS.append(_r - _prev)
    _prev = _r
_OFFSETS.append(_MOD - _prev)


def _next_wheel_candidate(n: int) -> int:
    r = n % _MOD
    i = bisect_right(_WHEEL, r)
    if i < len(_WHEEL) and _WHEEL[i] == r:
        return n
    delta = (_WHEEL[i] - r) if i < len(_WHEEL) else (_MOD - r + _WHEEL[0])
    return n + delta


def _classify(n: int, budget_ms: int = 2000) -> tuple[str, Counter]:
    if is_prime(n):
        return "prime", Counter({n: 1})
    result = do_factor(n, budget_ms=budget_ms)
    fs = []
    for p, e in result.factors.items():
        fs.extend([p] * e)
    fs.sort()
    cnt = Counter(fs)
    total = sum(cnt.values())
    if total == 2 and result.complete:
        return "semiprime", cnt
    return "composite", cnt


def _factor_str(cnt: Counter) -> str:
    return "*".join(
        f"{p}" if e == 1 else f"{p}^{e}"
        for p, e in sorted(cnt.items())
    )


def _open_for_append(path: str, append: bool):
    if path.endswith(".gz"):
        raw = open(path, "ab" if append else "wb")
        return io.TextIOWrapper(gzip.GzipFile(fileobj=raw, mode="ab" if append else "wb"),
                                encoding="utf-8", newline="")
    return open(path, "a" if append else "w", newline="", encoding="utf-8")


def _find_resume(path: str) -> int | None:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    try:
        if path.endswith(".gz"):
            fh = io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", newline="")
        else:
            fh = open(path, "r", encoding="utf-8", newline="")
        last = None
        with fh:
            for row in csv.DictReader(fh):
                last = row
        return int(last["n"]) if last else None
    except Exception:
        return None


_FIELDS = [
    "n", "classification", "omega_total", "omega_distinct",
    "smallest_prime_factor", "largest_prime_factor", "factorization",
]


def scan(start: int, stop: int, out_path: str,
         mode: str = "thin", diag_period: int = 1000,
         resume: bool = True, progress_every: int = 500_000,
         budget_ms: int = 2000) -> Iterator[dict]:
    """
    Scan integers from start to stop using the mod-210 wheel.
    Writes rows to out_path (CSV or .gz). Yields each row dict.

    mode: 'thin' (no diagnostics), 'full', 'sampled' (every diag_period)
    """
    append = False
    if resume:
        last_n = _find_resume(out_path)
        if last_n is not None:
            start = max(start, last_n + 1)
            append = True

    if start > stop:
        print("Nothing to do.", file=sys.stderr)
        return

    fh = _open_for_append(out_path, append)
    with fh:
        writer = csv.DictWriter(fh, fieldnames=_FIELDS)
        if not append:
            writer.writeheader()

        n = _next_wheel_candidate(start)
        samples = 0

        while n <= stop:
            for off in _OFFSETS:
                if n > stop:
                    break
                try:
                    kind, cnt = _classify(n, budget_ms=budget_ms)
                    row = {
                        "n": n,
                        "classification": kind,
                        "omega_total": sum(cnt.values()),
                        "omega_distinct": len(cnt),
                        "smallest_prime_factor": min(cnt) if cnt else n,
                        "largest_prime_factor": max(cnt) if cnt else n,
                        "factorization": _factor_str(cnt),
                    }
                    writer.writerow(row)
                    yield row
                except KeyboardInterrupt:
                    print(f"\nInterrupted — results saved to {out_path}", file=sys.stderr)
                    return
                except Exception as e:
                    print(f"[warn] n={n}: {e}", file=sys.stderr)

                samples += 1
                if samples % progress_every == 0:
                    print(f"[scan] n={n:,}", file=sys.stderr)
                    fh.flush()
                n += off

    print(f"Scan complete → {out_path}", file=sys.stderr)
