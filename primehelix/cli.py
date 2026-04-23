from __future__ import annotations
import csv
import json as _json
import sys

import click
from rich.console import Console
from rich.table import Table

from primehelix import __version__
from primehelix.analysis import scan_range, compare_summaries, build_time_series
from primehelix.schema import VALID_CLASSIFICATIONS
from primehelix.display.json_output import (
    build_json_result,
    label_entropy,
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
@click.version_option(__version__, prog_name="primehelix")
def main():
    pass


# -----------------------------
# helpers
# -----------------------------

def _parse_n(raw: str) -> int:
    try:
        return int(raw)
    except ValueError:
        raise click.BadParameter(f"expected an integer, got {raw!r}")


def _print_residue_profile(residue: dict):
    table = Table(title="residue profile")
    table.add_column("field")
    table.add_column("value")
    for key, value in residue.items():
        table.add_row(str(key), str(value))
    console.print(table)


def _print_structure_summary(summary, limit: int = 20):
    counts = summary.counts
    total = summary.total
    table = Table(title="structure summary")
    table.add_column("structure")
    table.add_column("count", justify="right")
    table.add_column("percent", justify="right")
    table.add_column("histogram")
    max_count = max(counts.values()) if counts else 1
    for label, count in counts.most_common(limit):
        pct = (count / total * 100.0) if total else 0.0
        bar = "█" * int((count / max_count) * 30)
        table.add_row(label, str(count), f"{pct:.2f}%", bar)
    console.print(table)
    console.print(f"[dim]total numbers scanned: {total}[/dim]")


def _print_compare_ranges(
    summary_a,
    summary_b,
    label_a: str,
    label_b: str,
    limit: int = 20,
    top_delta: int | None = None,
):
    rows = compare_summaries(summary_a, summary_b)
    if top_delta:
        rows = sorted(rows, key=lambda r: (-abs(r.delta), r.structure.lower()))[:top_delta]
    else:
        rows = sorted(rows, key=lambda r: (-(r.a_count + r.b_count), r.structure.lower()))[:limit]

    title = "range comparison" + (f" | top delta {top_delta}" if top_delta else "")
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
            row.structure,
            str(row.a_count),
            f'{row.a_percent:.2f}%',
            str(row.b_count),
            f'{row.b_percent:.2f}%',
            f'{row.delta:+.2f}%',
            row.ratio,
        )
    console.print(table)
    console.print(f"[dim]{label_a} total: {summary_a.total}[/dim]")
    console.print(f"[dim]{label_b} total: {summary_b.total}[/dim]")


def _write_csv(path: str, fieldnames: list[str], rows: list[dict]):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: str, payload: dict):
    with open(path, "w") as f:
        _json.dump(payload, f, indent=2)


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
    from .geometry.coil import CoilBalance, coil_footprint

    N = _parse_n(n)
    classification, result = do_classify(N)
    res = residue_profile(N, result.factors, classification=classification)

    coil_balance = None
    coil_data = None
    if classification == "semiprime":
        primes = sorted(result.factors.keys())
        if len(primes) == 2:
            coil_balance = CoilBalance(primes[0], primes[1], N)
            if coil or helix:
                coil_data = coil_footprint(N, primes[0], primes[1])

    if as_json:
        payload = build_json_result(
            result,
            command="classify",
            classification=classification,
            coil=coil_data,
            residue=res,
        )
        if coil_balance and not coil_data:
            payload["structure"] = structure_summary(classification, coil=coil_balance, residue=res)
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

    N = _parse_n(n)
    result = do_factor(N, budget_ms=budget)

    if as_json:
        payload = build_json_result(result, command="factor")
        payload.pop("classification", None)
        payload["steps"] = list(result.steps or []) if verbose else []
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
@click.option("--profile", is_flag=True, help="Include factorization method distribution")
@click.option("--only-classification",
              type=click.Choice(sorted(VALID_CLASSIFICATIONS), case_sensitive=False))
@click.option("--only-structure", type=str)
@click.option("--export-csv", "export_csv", default=None, type=str,
              help="Write label counts to CSV at this path")
@click.option("--export-json", "export_json_path", default=None, type=str,
              help="Write full JSON result to this file")
@click.option("--fast", "fast_mode", is_flag=True,
              help="Skip geometry — return classification-only labels (faster for large ranges)")
def structure_scan(start, stop, as_json, profile, only_classification, only_structure,
                   export_csv, export_json_path, fast_mode):
    span = stop - start
    summary = scan_range(
        start, stop,
        budget=2000,
        only_classification=only_classification,
        only_structure=only_structure,
        progress=span > 10_000 and not as_json,
        detail="classification" if fast_mode else "full",
    )

    payload = {
        "command": "structure-scan",
        "start": start,
        "stop": stop,
        "entropy": label_entropy(summary.counts, summary.total),
        **summary.to_json_dict(include_methods=profile),
    }
    if only_classification:
        payload["only_classification"] = only_classification
    if only_structure:
        payload["only_structure"] = only_structure

    if export_csv:
        rows = [
            {"label": label, "count": count,
             "percent": f"{count/summary.total*100:.4f}" if summary.total else "0"}
            for label, count in summary.counts.most_common()
        ]
        _write_csv(export_csv, ["label", "count", "percent"], rows)
        if not as_json:
            console.print(f"[green]CSV written to {export_csv}[/green]")

    if export_json_path:
        _write_json(export_json_path, payload)
        if not as_json:
            console.print(f"[green]JSON written to {export_json_path}[/green]")

    if as_json:
        print_json(payload)
        return

    _print_structure_summary(summary)

    if profile and summary.methods:
        methods = summary.methods
        mtotal = sum(methods.values())
        table = Table(title="factorization method profile")
        table.add_column("method")
        table.add_column("count", justify="right")
        table.add_column("percent", justify="right")
        for m, c in methods.most_common():
            table.add_row(m, str(c), f"{c/mtotal*100:.2f}%")
        console.print(table)


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
@click.option("--only-classification",
              type=click.Choice(sorted(VALID_CLASSIFICATIONS), case_sensitive=False))
@click.option("--only-structure", type=str)
@click.option("--export-csv", "export_csv", default=None, type=str,
              help="Write comparison rows to CSV at this path")
@click.option("--export-json", "export_json_path", default=None, type=str,
              help="Write full JSON result to this file")
@click.option("--fast", "fast_mode", is_flag=True,
              help="Skip geometry — return classification-only labels (faster for large ranges)")
def compare_ranges(
    a_start, a_stop, b_start, b_stop,
    top_delta, as_json, only_classification, only_structure,
    export_csv, export_json_path, fast_mode,
):
    span_a = a_stop - a_start
    span_b = b_stop - b_start
    detail = "classification" if fast_mode else "full"
    summary_a = scan_range(
        a_start, a_stop, budget=2000,
        only_classification=only_classification, only_structure=only_structure,
        progress=span_a > 10_000 and not as_json,
        detail=detail,
    )
    summary_b = scan_range(
        b_start, b_stop, budget=2000,
        only_classification=only_classification, only_structure=only_structure,
        progress=span_b > 10_000 and not as_json,
        detail=detail,
    )

    rows = compare_summaries(summary_a, summary_b)
    if top_delta:
        rows = sorted(rows, key=lambda r: (-abs(r.delta), r.structure.lower()))[:top_delta]

    ea = label_entropy(summary_a.counts, summary_a.total)
    eb = label_entropy(summary_b.counts, summary_b.total)
    payload = {
        "command": "compare-ranges",
        "a": {"start": a_start, "stop": a_stop, "entropy": ea, **summary_a.to_json_dict()},
        "b": {"start": b_start, "stop": b_stop, "entropy": eb, **summary_b.to_json_dict()},
        "entropy_delta": round(eb - ea, 4) if ea is not None and eb is not None else None,
        "rows": [r.to_dict() for r in rows],
    }
    if top_delta:
        payload["top_delta"] = top_delta
    if only_classification:
        payload["only_classification"] = only_classification
    if only_structure:
        payload["only_structure"] = only_structure

    if export_csv:
        _write_csv(
            export_csv,
            ["structure", "a_count", "a_percent", "b_count", "b_percent", "delta", "ratio"],
            [r.to_dict() for r in rows],
        )
        if not as_json:
            console.print(f"[green]CSV written to {export_csv}[/green]")

    if export_json_path:
        _write_json(export_json_path, payload)
        if not as_json:
            console.print(f"[green]JSON written to {export_json_path}[/green]")

    if as_json:
        print_json(payload)
        return

    _print_compare_ranges(
        summary_a, summary_b,
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
@click.option("--only-classification",
              type=click.Choice(sorted(VALID_CLASSIFICATIONS), case_sensitive=False))
@click.option("--only-structure", type=str)
@click.option("--export-csv", "export_csv", default=None, type=str,
              help="Write per-window series data to CSV at this path")
@click.option("--export-json", "export_json_path", default=None, type=str,
              help="Write full JSON result to this file")
@click.option("--fast", "fast_mode", is_flag=True,
              help="Skip geometry — return classification-only labels (faster for large ranges)")
def structure_time_series(
    start, stop, window, step, metric, top,
    plot_path, as_json, only_classification, only_structure,
    export_csv, export_json_path, fast_mode,
):
    if stop <= start:
        raise click.UsageError("stop must be greater than start")
    if window <= 0:
        raise click.UsageError("window must be positive")
    if step <= 0:
        raise click.UsageError("step must be positive")

    ts = build_time_series(
        start, stop,
        window=window, step=step,
        budget=2000, metric=metric, top=top,
        only_classification=only_classification,
        only_structure=only_structure,
        progress=not as_json,
        detail="classification" if fast_mode else "full",
    )

    if not ts.windows:
        raise click.UsageError("no windows generated")

    title = f"Structure Time Series [{start}, {stop})"
    ylabel = "Percent" if metric == "percent" else "Count"
    if only_classification:
        title += f" | class={only_classification}"
    if only_structure:
        title += f" | structure~{only_structure}"

    payload = {
        "command": "structure-time-series",
        "start": start,
        "stop": stop,
        "window": window,
        "step": step,
        "metric": metric,
        "top": top,
        "plot": plot_path,
        **ts.to_json_dict(),
    }
    if only_classification:
        payload["only_classification"] = only_classification
    if only_structure:
        payload["only_structure"] = only_structure

    if export_csv:
        csv_rows = []
        for i, win in enumerate(ts.windows):
            row = {"window": ts.window_labels[i], "start": win.start, "stop": win.stop}
            for label in ts.top_labels:
                row[label] = f"{ts.series_map[label][i]:.4f}"
            csv_rows.append(row)
        _write_csv(export_csv, ["window", "start", "stop"] + ts.top_labels, csv_rows)
        if not as_json:
            console.print(f"[green]CSV written to {export_csv}[/green]")

    if export_json_path:
        _write_json(export_json_path, payload)
        if not as_json:
            console.print(f"[green]JSON written to {export_json_path}[/green]")

    if plot_path:
        try:
            from .display.plots import save_structure_time_series_plot
        except ImportError as e:
            raise click.UsageError(str(e))
        save_structure_time_series_plot(
            series_map=ts.series_map,
            window_labels=ts.window_labels,
            output_path=plot_path,
            title=title,
            ylabel=ylabel,
        )
        if not as_json:
            console.print(f"[green]Plot written to {plot_path}[/green]")
    elif not as_json:
        for label in ts.top_labels:
            vals = ts.series_map[label]
            console.print(f"  [cyan]{label}[/cyan]: {[f'{v:.1f}' for v in vals]}")

    if as_json:
        print_json(payload)


if __name__ == "__main__":
    main()
