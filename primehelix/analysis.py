"""
Reusable analysis functions — scan, compare, time-series.

These return typed schema objects (ScanResult, list[CompareRow], TimeSeriesResult)
and have no CLI dependency, so they can be called from notebooks or scripts directly.
"""
from __future__ import annotations

import sys
import time
from collections import Counter

from .schema import VALID_CLASSIFICATIONS, CompareRow, ScanResult, TimeSeriesResult, WindowSummary


def _report_progress(i: int, n: int, span: int, t0: float) -> None:
    elapsed = time.monotonic() - t0
    pct = (i + 1) / span * 100
    rate = (i + 1) / elapsed if elapsed > 0 else 0
    eta = (span - i - 1) / rate if rate > 0 else 0
    sys.stderr.write(f"\r  {pct:5.1f}%  n={n:,}  {rate:,.0f}/s  eta {eta:.0f}s    ")
    sys.stderr.flush()


def scan_range(
    start: int,
    stop: int,
    budget: int = 2000,
    only_classification: str | None = None,
    only_structure: str | None = None,
    progress: bool = False,
    detail: str = "full",
) -> ScanResult:
    """Scan [start, stop), assign structure labels, return counts + method totals.

    detail="full"           — complete labels including balance tier and residue family (default)
    detail="classification" — classification-only labels (prime/semiprime/composite/invalid);
                              skips all geometry, ~9% faster on unfiltered scans, much faster
                              when combined with only_classification filtering
    """
    if only_classification is not None and only_classification.lower() not in VALID_CLASSIFICATIONS:
        raise ValueError(
            f"unknown classification {only_classification!r}; "
            f"valid values: {sorted(VALID_CLASSIFICATIONS)}"
        )

    from .core.factor import classify as do_classify

    _fast = detail == "classification"
    if not _fast:
        from .geometry.residue import residue_profile
        from .geometry.coil import CoilBalance
        from .display.json_output import structure_summary

    counts: Counter = Counter()
    method_counts: Counter = Counter()
    total = 0
    span = stop - start
    t0 = time.monotonic()
    _REPORT_EVERY = max(1, span // 100)
    _only_cls = only_classification.lower() if only_classification else None

    for i, n in enumerate(range(start, stop)):
        classification, result = do_classify(n, budget_ms=budget)
        cls_lower = classification.lower()

        # Early exit on classification filter — avoids all geometry for non-matching integers
        if _only_cls and cls_lower != _only_cls:
            if progress and span > 0 and (i + 1) % _REPORT_EVERY == 0:
                _report_progress(i, n, span, t0)
            continue

        if _fast:
            label = classification
        else:
            residue = residue_profile(n, result.factors, classification=classification)
            coil = None
            if cls_lower == "semiprime":
                primes = sorted(result.factors.keys())
                if len(primes) == 2:
                    coil = CoilBalance(primes[0], primes[1], n)
            label = structure_summary(classification, coil=coil, residue=residue)

        if only_structure and only_structure.lower() not in label.lower():
            if progress and span > 0 and (i + 1) % _REPORT_EVERY == 0:
                _report_progress(i, n, span, t0)
            continue

        counts[label] += 1
        method_counts[result.method or "trial"] += 1
        total += 1

        if progress and span > 0 and (i + 1) % _REPORT_EVERY == 0:
            _report_progress(i, n, span, t0)

    if progress and span > 0:
        sys.stderr.write("\r" + " " * 60 + "\r")
        sys.stderr.flush()

    return ScanResult(total=total, counts=counts, methods=method_counts)


def compare_summaries(summary_a: ScanResult, summary_b: ScanResult) -> list[CompareRow]:
    """Build per-label comparison rows from two scan_range results."""
    total_a = summary_a.total
    total_b = summary_b.total

    rows = []
    for s in set(summary_a.counts.keys()) | set(summary_b.counts.keys()):
        a = summary_a.counts.get(s, 0)
        b = summary_b.counts.get(s, 0)
        ap = (a / total_a * 100.0) if total_a else 0.0
        bp = (b / total_b * 100.0) if total_b else 0.0
        delta = bp - ap

        if ap == 0.0 and bp == 0.0:
            ratio = "1.00x"
        elif ap == 0.0:
            ratio = "new"
        else:
            ratio = f"{(bp / ap):.2f}x"

        rows.append(CompareRow(
            structure=s,
            a_count=a, a_percent=ap,
            b_count=b, b_percent=bp,
            delta=delta, ratio=ratio,
        ))

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
    detail: str = "full",
) -> TimeSeriesResult:
    """Aggregate structure distributions across sliding windows."""
    windows: list[WindowSummary] = []
    cursor = start
    while cursor < stop:
        win_start = cursor
        win_stop = min(cursor + window, stop)
        if win_stop <= win_start:
            break
        win_span = win_stop - win_start
        scan = scan_range(
            win_start, win_stop,
            budget=budget,
            only_classification=only_classification,
            only_structure=only_structure,
            progress=progress and win_span > 10_000,
            detail=detail,
        )
        windows.append(WindowSummary(
            start=win_start,
            stop=win_stop,
            label=f"[{win_start},{win_stop})",
            scan=scan,
        ))
        cursor += step

    aggregate_scores: dict[str, float] = {}
    for win in windows:
        for label, count in win.scan.counts.items():
            value = (count / win.scan.total * 100.0) if (metric == "percent" and win.scan.total) else float(count)
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
        window_labels.append(win.label)
        for label in top_labels:
            count = win.scan.counts.get(label, 0)
            value = (count / win.scan.total * 100.0) if (metric == "percent" and win.scan.total) else float(count)
            series_map[label].append(value)

    return TimeSeriesResult(
        windows=windows,
        top_labels=top_labels,
        series_map=series_map,
        window_labels=window_labels,
    )
