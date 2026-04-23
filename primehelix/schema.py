"""
Public data schema for primehelix analysis results.

These dataclasses are the stable return types of analysis.scan_range,
analysis.compare_summaries, and analysis.build_time_series.
Each class carries a to_json_dict() / to_dict() method so cli.py only
needs to add command-level envelope fields (command, start, stop, …).
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

VALID_CLASSIFICATIONS: frozenset[str] = frozenset({"prime", "semiprime", "composite", "invalid"})


@dataclass
class ScanResult:
    """Result of scanning a range with scan_range()."""
    total: int
    counts: Counter
    methods: Counter

    def to_json_dict(self, *, include_methods: bool = False) -> dict:
        d: dict = {"total": self.total, "counts": dict(self.counts)}
        if include_methods:
            d["methods"] = dict(self.methods)
        return d


@dataclass
class CompareRow:
    """One label's comparison between two ranges."""
    structure: str
    a_count: int
    a_percent: float
    b_count: int
    b_percent: float
    delta: float
    ratio: str

    def to_dict(self) -> dict:
        return {
            "structure": self.structure,
            "a_count": self.a_count,
            "a_percent": self.a_percent,
            "b_count": self.b_count,
            "b_percent": self.b_percent,
            "delta": self.delta,
            "ratio": self.ratio,
        }


@dataclass
class WindowSummary:
    """One window's scan result within a time series."""
    start: int
    stop: int
    label: str
    scan: ScanResult

    def to_dict(self) -> dict:
        return {"start": self.start, "stop": self.stop, "label": self.label, "total": self.scan.total}


@dataclass
class TimeSeriesResult:
    """Result of build_time_series()."""
    windows: list[WindowSummary]
    top_labels: list[str]
    series_map: dict[str, list[float]]
    window_labels: list[str]

    def to_json_dict(self) -> dict:
        return {
            "labels": self.top_labels,
            "windows": [w.to_dict() for w in self.windows],
            "series": self.series_map,
        }
