"""
Rich terminal output for all primehelix commands.
"""
from __future__ import annotations
import math

from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.text import Text
from rich import box

console = Console()


def _factor_str(factors: dict[int, int]) -> str:
    parts = []
    for p in sorted(factors):
        e = factors[p]
        parts.append(str(p) if e == 1 else f"{p}^{e}")
    return " × ".join(parts)


def print_classify(n: int, classification: str, result) -> None:
    color = {"prime": "green", "semiprime": "yellow", "composite": "red"}.get(classification, "white")
    label = Text(classification.upper(), style=f"bold {color}")

    table = Table(box=box.ROUNDED, show_header=False, padding=(0, 1))
    table.add_column("Field", style="dim")
    table.add_column("Value")

    table.add_row("n", str(n))
    table.add_row("bits", str(n.bit_length()))
    table.add_row("digits", str(len(str(n))))
    table.add_row("classification", label)
    table.add_row("factorization", _factor_str(result.factors))
    table.add_row("method", result.method)
    table.add_row("time", f"{result.elapsed_ms:.1f} ms")

    if not result.complete:
        table.add_row("note", Text("⚠ factorization incomplete (timeout)", style="yellow"))

    console.print(Panel(table, title="[bold]classify[/bold]", border_style=color))

    if result.steps:
        console.print("[dim]Pipeline steps:[/dim]")
        for s in result.steps:
            console.print(f"  [dim]·[/dim] {s}")


def print_factor(result) -> None:
    color = "green" if result.complete else "yellow"
    table = Table(box=box.ROUNDED, show_header=False, padding=(0, 1))
    table.add_column("Field", style="dim")
    table.add_column("Value")

    table.add_row("n", str(result.n))
    table.add_row("bits", str(result.n.bit_length()))
    table.add_row("factorization", _factor_str(result.factors))
    table.add_row("method", result.method)
    table.add_row("time", f"{result.elapsed_ms:.1f} ms")
    table.add_row("complete", "yes" if result.complete else Text("no (timeout)", style="yellow"))

    console.print(Panel(table, title="[bold]factor[/bold]", border_style=color))

    if result.steps:
        console.print("[dim]Pipeline steps:[/dim]")
        for s in result.steps:
            console.print(f"  [dim]·[/dim] {s}")


def print_coil(fp, show_signature: bool = False) -> None:
    table = Table(box=box.ROUNDED, show_header=False, padding=(0, 1))
    table.add_column("Field", style="dim")
    table.add_column("Value")

    table.add_row("n", str(fp.n))
    table.add_row("p (small factor)", str(fp.p))
    table.add_row("q (large factor)", str(fp.q))
    table.add_row("bit_gap", str(fp.bit_gap))
    table.add_row("balance |p-q|/√n", f"{fp.balance:.6f}")
    table.add_row("", "")
    table.add_row("d(n → q)", f"{fp.d_n_to_q:.6f}")
    table.add_row("d(q → p)", f"{fp.d_q_to_p:.6f}")
    table.add_row("d(p → 1)", f"{fp.d_p_to_1:.6f}")
    table.add_row("", "")
    table.add_row("normalized (f1,f2,f3)", f"({fp.f1:.4f}, {fp.f2:.4f}, {fp.f3:.4f})")
    table.add_row("slope s1", f"{fp.s1:.4e}")
    table.add_row("slope s2", f"{fp.s2:.4e}")
    table.add_row("slope s3", f"{fp.s3:.4e}")

    if fp.bit_gap == 0:
        table.add_row("type", Text("balanced semiprime (RSA-like)", style="green"))
    else:
        table.add_row("type", Text("lopsided semiprime", style="yellow"))

    if show_signature:
        table.add_row("", "")
        table.add_row("sig (geom)", fp.sig_geom[:32] + "…")
        table.add_row("sig (invariant)", fp.sig_invariant[:32] + "…")

    console.print(Panel(table, title="[bold]coil footprint[/bold]", border_style="cyan"))


def print_bitbucket(bb, density_table=None) -> None:
    table = Table(box=box.ROUNDED, show_header=False, padding=(0, 1))
    table.add_column("Field", style="dim")
    table.add_column("Value")

    table.add_row("n", str(bb.n))
    table.add_row("bit-length k", str(bb.k))
    table.add_row("bucket range", f"[{bb.lo}, {bb.hi}]")
    table.add_row("bucket width", str(bb.width))
    table.add_row("offset", str(bb.offset))
    table.add_row("offset (normalized)", f"{bb.offset_norm:.6f}  (0=bucket start, 1=bucket end)")

    console.print(Panel(table, title="[bold]bit-bucket[/bold]", border_style="blue"))

    if density_table:
        dt = Table(title="Prime Density by Bit-Length", box=box.SIMPLE_HEAVY)
        dt.add_column("k", justify="right")
        dt.add_column("primes", justify="right")
        dt.add_column("density", justify="right")
        dt.add_column("asymptotic 1/(k·ln2)", justify="right")
        dt.add_column("ratio", justify="right")

        for row in density_table:
            style = "bold cyan" if row.k == bb.k else ""
            marker = " ← here" if row.k == bb.k else ""
            dt.add_row(
                str(row.k) + marker,
                f"{row.prime_count:,}",
                f"{row.density:.6f}",
                f"{row.asymptotic:.6f}",
                f"{row.ratio:.3f}",
                style=style,
            )
        console.print(dt)


def print_tangent(eq, ts, ideal=None) -> None:
    table = Table(box=box.ROUNDED, show_header=False, padding=(0, 1))
    table.add_column("Field", style="dim")
    table.add_column("Value")

    table.add_row("[bold]Equal split[/bold] L = 2√n", "")
    table.add_row("  L", eq.L[:40])
    table.add_row("  half (√n)", eq.half[:40])
    table.add_row("  remainder (half² - n)", eq.remainder[:40])
    table.add_row("", "")
    table.add_row("[bold]Tangent split[/bold] L = n+1", "")
    table.add_row("  L", str(ts.L))
    table.add_row("  discriminant (n-1)²", str(ts.discriminant))
    table.add_row("  √disc exact?", "yes" if ts.sqrt_disc_exact else "no")

    if ideal:
        table.add_row("", "")
        table.add_row("[bold]Ideal split[/bold] from factors", "")
        table.add_row(f"  L = p+q", str(ideal.L))
        table.add_row(f"  discriminant (p-q)²", str(ideal.discriminant))
        table.add_row("  roots (p, q)", f"{ideal.roots[0]}, {ideal.roots[1]}")

    console.print(Panel(table, title="[bold]tangent diagnostics[/bold]", border_style="magenta"))


def print_ecm(n: int, factor: int | None, elapsed_ms: float, B1: int, curves: int) -> None:
    if factor:
        cofactor = n // factor
        table = Table(box=box.ROUNDED, show_header=False, padding=(0, 1))
        table.add_column("Field", style="dim")
        table.add_column("Value")
        table.add_row("n", str(n))
        table.add_row("factor", str(min(factor, cofactor)))
        table.add_row("cofactor", str(max(factor, cofactor)))
        table.add_row("B1", f"{B1:,}")
        table.add_row("curves", str(curves))
        table.add_row("time", f"{elapsed_ms:.1f} ms")
        console.print(Panel(table, title="[bold]ECM[/bold]", border_style="green"))
    else:
        console.print(Panel(
            f"[yellow]No factor found[/yellow]\n"
            f"n={n}\nB1={B1:,}  curves={curves}  time={elapsed_ms:.1f}ms",
            title="[bold]ECM[/bold]", border_style="yellow"
        ))


def print_qs(n: int, factor: int | None, elapsed_ms: float) -> None:
    if factor:
        cofactor = n // factor
        table = Table(box=box.ROUNDED, show_header=False, padding=(0, 1))
        table.add_column("Field", style="dim")
        table.add_column("Value")
        table.add_row("n", str(n))
        table.add_row("factor", str(min(factor, cofactor)))
        table.add_row("cofactor", str(max(factor, cofactor)))
        table.add_row("time", f"{elapsed_ms:.1f} ms")
        console.print(Panel(table, title="[bold]Quadratic Sieve[/bold]", border_style="green"))
    else:
        console.print(Panel(
            f"[yellow]No factor found[/yellow]\nn={n}  time={elapsed_ms:.1f}ms",
            title="[bold]Quadratic Sieve[/bold]", border_style="yellow"
        ))
