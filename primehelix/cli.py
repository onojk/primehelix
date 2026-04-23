from __future__ import annotations
import csv

import click
from rich.console import Console
from rich.table import Table

from primehelix.display.json_output import (
    build_json_result,
    print_json,
    structure_summary,
)

console = Console()


@click.group(help="""
primehelix — structural analysis for integers

Explore how numbers are built, not just what they are.

Examples:
  primehelix classify 1300039 --coil --json
  primehelix structure-scan --start 10 --stop 100 --json
  primehelix compare-ranges --a-start 10 --a-stop 20 --b-start 20 --b-stop 30 --top-delta 5
  primehelix structure-time-series --start 1 --stop 100000 --window 10000 --step 10000 --only-classification semiprime --plot ts.png
""")
@click.version_option("0.1.2", prog_name="primehelix")
def main():
    pass


# -----------------------------
# helpers
# -----------------------------

def _print_residue_profile(residue: dict):
    table = Table(title="residue profile")
    table.add_column("field")
    table.add_column("value")

    for key, value in residue.items():
        table.add_row(str(key), str(value))

    console.print(table)


def _print_structure_summary(summary: dict, limit: int = 20):
    counts = summary["counts"]
    total = summary["total"]

    table = Table(title="structure summary")
    table.add_column("structure")
    table.add_column("count", justify="right")
    table.add_column("percent", justify="right")
    table.add_column("histogram")

    max_count = max(counts.values()) if counts else 1

    for label, count in counts.most_common(limit):
        pct = (count / total * 100.0) if total else 0.0
        bar_len = int((count / max_count) * 30)
        bar = "█" * bar_len
        table.add_row(label, str(count), f"{pct:.2f}%", bar)

    console.print(table)
    console.print(f"[dim]total numbers scanned: {total}[/dim]")


def _matches_filters(
    classification: str,
    label: str,
    only_classification: str | None,
    only_structure: str | None,
) -> bool:
    if only_classification and classification.lower() != only_classification.lower():
        return False
    if only_structure and only_structure.lower() not in label.lower():
        return False
    return True


def _summarize_filtered_range(
    start: int,
    stop: int,
    budget: int,
    only_classification: str | None = None,
    only_structure: str | None = None,
    progress: bool = False,
):
    import sys
    import time
    from collections import Counter
    from .core.factor import classify as do_classify
    from .geometry.residue import residue_profile

    counts = Counter()
    total = 0
    span = stop - start
    t0 = time.monotonic()
    _REPORT_EVERY = max(1, span // 100)  # ~1% increments

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

        if not _matches_filters(classification, label, only_classification, only_structure):
            pass
        else:
            counts[label] += 1
            total += 1

        if progress and span > 0 and (i + 1) % _REPORT_EVERY == 0:
            pct = (i + 1) / span * 100
            elapsed = time.monotonic() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            eta = (span - i - 1) / rate if rate > 0 else 0
            sys.stderr.write(
                f"\r  {pct:5.1f}%  n={n:,}  {rate:,.0f}/s  eta {eta:.0f}s    "
            )
            sys.stderr.flush()

    if progress and span > 0:
        sys.stderr.write("\r" + " " * 60 + "\r")
        sys.stderr.flush()

    return {
        "total": total,
        "counts": counts,
    }


def _compare_rows(summary_a: dict, summary_b: dict):
    counts_a = summary_a["counts"]
    counts_b = summary_b["counts"]
    total_a = summary_a["total"]
    total_b = summary_b["total"]

    labels = set(counts_a.keys()) | set(counts_b.keys())
    rows = []

    for s in labels:
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

        rows.append(
            {
                "structure": s,
                "a_count": a,
                "a_percent": ap,
                "b_count": b,
                "b_percent": bp,
                "delta": delta,
                "ratio": ratio,
            }
        )

    return rows


def _print_compare_ranges(
    summary_a: dict,
    summary_b: dict,
    label_a: str,
    label_b: str,
    limit: int = 20,
    top_delta: int | None = None,
):
    total_a = summary_a["total"]
    total_b = summary_b["total"]
    rows = _compare_rows(summary_a, summary_b)

    if top_delta:
        rows = sorted(
            rows,
            key=lambda r: (-abs(r["delta"]), r["structure"].lower())
        )[:top_delta]
    else:
        rows = sorted(
            rows,
            key=lambda r: (
                -(r["a_count"] + r["b_count"]),
                r["structure"].lower(),
            )
        )[:limit]

    title = "range comparison"
    if top_delta:
        title += f" | top delta {top_delta}"

    table = Table(title=title)
    table.add_column("structure", no_wrap=True)
    table.add_column(f"{label_a} count", justify="right")
    table.add_column(f"{label_a} %", justify="right")
    table.add_column(f"{label_b} count", justify="right")
    table.add_column(f"{label_b} %", justify="right")
    table.add_column("delta", justify="right")
    table.add_column("ratio", justify="right")

    for row in rows:
        table.add_row(
            row["structure"],
            str(row["a_count"]),
            f'{row["a_percent"]:.2f}%',
            str(row["b_count"]),
            f'{row["b_percent"]:.2f}%',
            f'{row["delta"]:+.2f}%',
            row["ratio"],
        )

    console.print(table)
    console.print(f"[dim]{label_a} total: {total_a}[/dim]")
    console.print(f"[dim]{label_b} total: {total_b}[/dim]")


# -----------------------------
# classify
# -----------------------------

@main.command()
@click.argument("n", type=str)
@click.option("--coil", is_flag=True, help="Show conical helix footprint")
@click.option("--helix", is_flag=True, help="Render ASCII double-helix visualization")
@click.option("--residue", is_flag=True, help="Show residue/arithmetic-family profile")
@click.option("--json", "as_json", is_flag=True, help="Output result as JSON")
def classify(n, coil, helix, residue, as_json):
    from .core.factor import classify as do_classify
    from .geometry.residue import residue_profile
    from .geometry.coil import coil_footprint

    N = int(n)
    classification, result = do_classify(N)
    res = residue_profile(N, result.factors, classification=classification)

    coil_data = None
    if (coil or helix) and classification == "semiprime":
        primes = sorted(result.factors.keys())
        if len(primes) == 2:
            coil_data = coil_footprint(N, primes[0], primes[1])

    if as_json:
        payload = build_json_result(
            result,
            command="classify",
            classification=classification,
            coil=coil_data,
            residue=res,
        )
        print_json(payload)
        return

    console.print(f"{N} → {classification}")
    if residue:
        _print_residue_profile(res)
    if coil and coil_data:
        from .display.output import print_coil
        print_coil(coil_data)
    if helix:
        from .display.ascii_helix import print_ascii_helix
        if coil_data:
            print_ascii_helix(coil_data)
        else:
            console.print("[yellow]Helix requires a semiprime.[/yellow]")


# -----------------------------
# factor
# -----------------------------

@main.command()
@click.argument("n", type=str)
@click.option("--verbose", is_flag=True, help="Show pipeline steps")
@click.option("--budget", default=10000, show_default=True, type=int, help="Time budget in ms")
@click.option("--json", "as_json", is_flag=True, help="Output result as JSON")
def factor(n, verbose, budget, as_json):
    from .core.factor import factor as do_factor

    N = int(n)
    result = do_factor(N, budget_ms=budget)

    if as_json:
        payload = build_json_result(result, command="factor")
        if verbose:
            payload["steps"] = list(result.steps or [])
        print_json(payload)
        return

    from .display.output import print_factor
    print_factor(result, verbose=verbose)


# -----------------------------
# structure-scan
# -----------------------------

@main.command("structure-scan")
@click.option("--start", required=True, type=int)
@click.option("--stop", required=True, type=int)
@click.option("--json", "as_json", is_flag=True)
@click.option("--only-classification", type=str)
@click.option("--only-structure", type=str)
def structure_scan(start, stop, as_json, only_classification, only_structure):
    span = stop - start
    summary = _summarize_filtered_range(
        start,
        stop,
        budget=2000,
        only_classification=only_classification,
        only_structure=only_structure,
        progress=span > 10_000 and not as_json,
    )

    if as_json:
        payload = {
            "command": "structure-scan",
            "start": start,
            "stop": stop,
            "total": summary["total"],
            "counts": dict(summary["counts"]),
        }
        if only_classification:
            payload["only_classification"] = only_classification
        if only_structure:
            payload["only_structure"] = only_structure
        print_json(payload)
        return

    _print_structure_summary(summary)


# -----------------------------
# compare-ranges
# -----------------------------

@main.command("compare-ranges")
@click.option("--a-start", required=True, type=int)
@click.option("--a-stop", required=True, type=int)
@click.option("--b-start", required=True, type=int)
@click.option("--b-stop", required=True, type=int)
@click.option("--top-delta", type=int)
@click.option("--json", "as_json", is_flag=True)
@click.option("--only-classification", type=str)
@click.option("--only-structure", type=str)
def compare_ranges(
    a_start,
    a_stop,
    b_start,
    b_stop,
    top_delta,
    as_json,
    only_classification,
    only_structure,
):
    span_a = a_stop - a_start
    span_b = b_stop - b_start
    summary_a = _summarize_filtered_range(
        a_start,
        a_stop,
        budget=2000,
        only_classification=only_classification,
        only_structure=only_structure,
        progress=span_a > 10_000 and not as_json,
    )
    summary_b = _summarize_filtered_range(
        b_start,
        b_stop,
        budget=2000,
        only_classification=only_classification,
        only_structure=only_structure,
        progress=span_b > 10_000 and not as_json,
    )

    rows = _compare_rows(summary_a, summary_b)
    if top_delta:
        rows = sorted(
            rows,
            key=lambda r: (-abs(r["delta"]), r["structure"].lower())
        )[:top_delta]

    if as_json:
        payload = {
            "command": "compare-ranges",
            "a": {
                "start": a_start,
                "stop": a_stop,
                "total": summary_a["total"],
                "counts": dict(summary_a["counts"]),
            },
            "b": {
                "start": b_start,
                "stop": b_stop,
                "total": summary_b["total"],
                "counts": dict(summary_b["counts"]),
            },
            "rows": rows,
        }

        if top_delta:
            payload["top_delta"] = top_delta
        if only_classification:
            payload["only_classification"] = only_classification
        if only_structure:
            payload["only_structure"] = only_structure

        print_json(payload)
        return

    _print_compare_ranges(
        summary_a,
        summary_b,
        label_a=f"[{a_start}, {a_stop})",
        label_b=f"[{b_start}, {b_stop})",
        top_delta=top_delta,
    )


# -----------------------------
# structure-time-series
# -----------------------------

@main.command("structure-time-series")
@click.option("--start", required=True, type=int)
@click.option("--stop", required=True, type=int)
@click.option("--window", default=10000, show_default=True, type=int)
@click.option("--step", default=10000, show_default=True, type=int)
@click.option("--metric", type=click.Choice(["percent", "count"]), default="percent", show_default=True)
@click.option("--top", default=5, show_default=True, type=int, help="Number of structure series to plot")
@click.option("--plot", "plot_path", default=None, type=str, help="Save plot to this path (PNG)")
@click.option("--json", "as_json", is_flag=True)
@click.option("--only-classification", type=str)
@click.option("--only-structure", type=str)
def structure_time_series(
    start,
    stop,
    window,
    step,
    metric,
    top,
    plot_path,
    as_json,
    only_classification,
    only_structure,
):
    from .display.plots import save_structure_time_series_plot

    if stop <= start:
        raise click.UsageError("stop must be greater than start")
    if window <= 0:
        raise click.UsageError("window must be positive")
    if step <= 0:
        raise click.UsageError("step must be positive")

    windows = []
    cursor = start
    while cursor < stop:
        win_start = cursor
        win_stop = min(cursor + window, stop)
        if win_stop <= win_start:
            break

        win_span = win_stop - win_start
        summary = _summarize_filtered_range(
            win_start,
            win_stop,
            budget=2000,
            only_classification=only_classification,
            only_structure=only_structure,
            progress=win_span > 10_000 and not as_json,
        )
        windows.append(
            {
                "start": win_start,
                "stop": win_stop,
                "label": f"[{win_start},{win_stop})",
                "summary": summary,
            }
        )
        cursor += step

    if not windows:
        raise click.UsageError("no windows generated")

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

    title = f"Structure Time Series [{start}, {stop})"
    ylabel = "Percent" if metric == "percent" else "Count"
    if only_classification:
        title += f" | class={only_classification}"
    if only_structure:
        title += f" | structure~{only_structure}"

    if plot_path:
        save_structure_time_series_plot(
            series_map=series_map,
            window_labels=window_labels,
            output_path=plot_path,
            title=title,
            ylabel=ylabel,
        )
        if not as_json:
            console.print(f"[green]Plot written to {plot_path}[/green]")
    elif not as_json:
        for label in top_labels:
            vals = series_map[label]
            console.print(f"  [cyan]{label}[/cyan]: {[f'{v:.1f}' for v in vals]}")

    if as_json:
        payload = {
            "command": "structure-time-series",
            "start": start,
            "stop": stop,
            "window": window,
            "step": step,
            "metric": metric,
            "top": top,
            "plot": plot_path,
            "labels": top_labels,
            "windows": [
                {
                    "start": win["start"],
                    "stop": win["stop"],
                    "label": win["label"],
                    "total": win["summary"]["total"],
                }
                for win in windows
            ],
            "series": series_map,
        }
        if only_classification:
            payload["only_classification"] = only_classification
        if only_structure:
            payload["only_structure"] = only_structure
        print_json(payload)


if __name__ == "__main__":
    main()
