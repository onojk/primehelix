"""
Public data schema for primehelix analysis results.

These dataclasses are the stable return types of analysis.scan_range,
analysis.compare_summaries, and analysis.build_time_series.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field


@dataclass
class ScanResult:
    """Result of scanning a range with scan_range()."""
    total: int
    counts: Counter
    methods: Counter


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


@dataclass
class TimeSeriesResult:
    """Result of build_time_series()."""
    windows: list[WindowSummary]
    top_labels: list[str]
    series_map: dict[str, list[float]]
    window_labels: list[str]
