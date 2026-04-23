"""
Reusable analysis functions — scan, compare, time-series.

These return plain Python dicts/lists and have no CLI dependency,
so they can be called from notebooks or scripts without any Click involvement.
"""
from __future__ import annotations

import sys
import time
from collections import Counter


def scan_range(
    start: int,
    stop: int,
    budget: int = 2000,
    only_classification: str | None = None,
    only_structure: str | None = None,
    progress: bool = False,
) -> dict:
    """Scan [start, stop), assign structure labels, return counts + method totals.

    Returns:
        {"total": int, "counts": Counter, "methods": Counter}
    """
    from .core.factor import classify as do_classify
    from .geometry.residue import residue_profile
    from .display.json_output import structure_summary

    counts: Counter = Counter()
    method_counts: Counter = Counter()
    total = 0
    span = stop - start
    t0 = time.monotonic()
    _REPORT_EVERY = max(1, span // 100)

    for i, n in enumerate(range(start, stop)):
        classification, result = do_classify(n, budget_ms=budget)
        residue = residue_profile(n, result.factors, classification=classification)

        coil = None
        if classification.lower() == "semiprime":
            from .geometry.coil import CoilBalance
            primes = sorted(result.factors.keys())
            if len(primes) == 2:
                coil = CoilBalance(primes[0], primes[1], n)

        label = structure_summary(classification, coil=coil, residue=residue)

        if only_classification and classification.lower() != only_classification.lower():
            pass
        elif only_structure and only_structure.lower() not in label.lower():
            pass
        else:
            counts[label] += 1
            method_counts[result.method or "trial"] += 1
            total += 1

        if progress and span > 0 and (i + 1) % _REPORT_EVERY == 0:
            pct = (i + 1) / span * 100
            elapsed = time.monotonic() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            eta = (span - i - 1) / rate if rate > 0 else 0
            sys.stderr.write(f"\r  {pct:5.1f}%  n={n:,}  {rate:,.0f}/s  eta {eta:.0f}s    ")
            sys.stderr.flush()

    if progress and span > 0:
        sys.stderr.write("\r" + " " * 60 + "\r")
        sys.stderr.flush()

    return {"total": total, "counts": counts, "methods": method_counts}


def compare_summaries(summary_a: dict, summary_b: dict) -> list[dict]:
    """Build per-label comparison rows from two scan_range results.

    Each row: structure, a_count, a_percent, b_count, b_percent, delta, ratio
    """
    counts_a = summary_a["counts"]
    counts_b = summary_b["counts"]
    total_a = summary_a["total"]
    total_b = summary_b["total"]

    rows = []
    for s in set(counts_a.keys()) | set(counts_b.keys()):
        a = counts_a.get(s, 0)
        b = counts_b.get(s, 0)
        ap = (a / total_a * 100.0) if total_a else 0.0
        bp = (b / total_b * 100.0) if total_b else 0.0
        delta = bp - ap

        if ap == 0.0 and bp == 0.0:
            ratio = "1.00x"
        elif ap == 0.0:
            ratio = "new"
        else:
            ratio = f"{(bp / ap):.2f}x"

        rows.append({
            "structure": s,
            "a_count": a,
            "a_percent": ap,
            "b_count": b,
            "b_percent": bp,
            "delta": delta,
            "ratio": ratio,
        })

    return rows


def build_time_series(
    start: int,
    stop: int,
    window: int = 10000,
    step: int = 10000,
    budget: int = 2000,
    metric: str = "percent",
    top: int = 5,
    only_classification: str | None = None,
    only_structure: str | None = None,
    progress: bool = False,
) -> dict:
    """Aggregate structure distributions across sliding windows.

    Returns:
        {
            "windows": [{"start", "stop", "label", "summary"}, ...],
            "top_labels": [...],
            "series_map": {label: [float, ...]},
            "window_labels": [...],
        }
    """
    windows = []
    cursor = start
    while cursor < stop:
        win_start = cursor
        win_stop = min(cursor + window, stop)
        if win_stop <= win_start:
            break
        win_span = win_stop - win_start
        summary = scan_range(
            win_start, win_stop,
            budget=budget,
            only_classification=only_classification,
            only_structure=only_structure,
            progress=progress and win_span > 10_000,
        )
        windows.append({
            "start": win_start,
            "stop": win_stop,
            "label": f"[{win_start},{win_stop})",
            "summary": summary,
        })
        cursor += step

    aggregate_scores: dict[str, float] = {}
    for win in windows:
        counts = win["summary"]["counts"]
        total = win["summary"]["total"]
        for label, count in counts.items():
            value = (count / total * 100.0) if (metric == "percent" and total) else float(count)
            aggregate_scores[label] = aggregate_scores.get(label, 0.0) + value

    top_labels = [
        label for label, _ in sorted(
            aggregate_scores.items(),
            key=lambda kv: (-kv[1], kv[0].lower())
        )[:top]
    ]

    series_map: dict[str, list[float]] = {label: [] for label in top_labels}
    window_labels: list[str] = []

    for win in windows:
        counts = win["summary"]["counts"]
        total = win["summary"]["total"]
        window_labels.append(win["label"])
        for label in top_labels:
            count = counts.get(label, 0)
            value = (count / total * 100.0) if (metric == "percent" and total) else float(count)
            series_map[label].append(value)

    return {
        "windows": windows,
        "top_labels": top_labels,
        "series_map": series_map,
        "window_labels": window_labels,
    }
