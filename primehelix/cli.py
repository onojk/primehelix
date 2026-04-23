"""
primehelix — unified console toolkit for prime number theory,
integer factorization, and geometric number exploration.
"""
from __future__ import annotations
import time
import sys

import click
from rich.console import Console
from rich.table import Table
from primehelix.display.json_output import build_json_result, print_json

console = Console()


@click.group()
@click.version_option("0.1.1", prog_name="primehelix")
def main():
    """
    primehelix: prime theory, factorization, and geometric number exploration.

    \b
    Commands:
      classify   Classify n as prime / semiprime / composite
      factor     Full factoring pipeline (trial→p-1→p+1→rho→ecm→qs)
      coil       Conical helix footprint for semiprimes
      bitbucket  Bit-bucket placement and prime density analysis
      scan       Wheel-accelerated range scanner to CSV
      ecm        Lenstra ECM — elliptic curve factoring
      qs         Quadratic Sieve factoring

    \b
    Output modes:
      default    Rich human-readable console output
      --json     Machine-readable JSON for scripting and automation
    """


def _print_residue_profile(residue: dict):
    table = Table(title="residue profile")
    table.add_column("field")
    table.add_column("value")

    for key, value in residue.items():
        table.add_row(str(key), str(value))

    console.print(table)


@main.command()
@click.argument("n", type=str)
@click.option("--coil", is_flag=True, help="Show conical helix footprint (semiprimes only)")
@click.option("--tangent", is_flag=True, help="Show tangent-split diagnostics")
@click.option("--residue", "show_residue", is_flag=True, help="Show arithmetic residue-family profile")
@click.option("--budget", default=10000, show_default=True, help="Time budget in ms")
@click.option("--json", "as_json", is_flag=True, help="Output machine-readable JSON")
def classify(n: str, coil: bool, tangent: bool, show_residue: bool, budget: int, as_json: bool):
    """
    Classify N as prime, semiprime, or composite.

    Use --json for machine-readable output.
    Use --coil with semiprimes to include geometric structure.
    Use --residue to include arithmetic family profiling.
    """
    from .core.factor import classify as do_classify
    from .display.output import print_classify, print_coil, print_tangent
    from .geometry.residue import residue_profile

    try:
        N = int(n.strip().replace(",", "").replace("_", ""))
    except ValueError:
        console.print("[red]Error:[/red] N must be an integer.")
        sys.exit(1)

    if N < 0:
        console.print("[red]Error:[/red] N must be non-negative.")
        sys.exit(1)

    classification, result = do_classify(N, budget_ms=budget)

    residue = residue_profile(N, result.factors, classification=classification)

    if as_json:
        coil_data = None

        if coil and classification.lower() == "semiprime":
            from .geometry.coil import coil_footprint
            primes = sorted(result.factors.keys())
            if len(primes) == 2:
                coil_data = coil_footprint(N, primes[0], primes[1])
            elif len(primes) == 1:
                p = primes[0]
                coil_data = coil_footprint(N, p, p)

        payload = build_json_result(
            result,
            command="classify",
            classification=classification,
            coil=coil_data,
            residue=residue,
        )
        print_json(payload)
        return

    print_classify(N, classification, result)

    if show_residue:
        _print_residue_profile(residue)

    if coil and classification.lower() == "semiprime":
        from .geometry.coil import coil_footprint
        primes = sorted(result.factors.keys())
        if len(primes) == 2:
            fp = coil_footprint(N, primes[0], primes[1])
            print_coil(fp)
        elif len(primes) == 1:
            p = primes[0]
            fp = coil_footprint(N, p, p)
            print_coil(fp)

    if tangent:
        from .geometry.tangent import equal_split, tangent_split, ideal_split
        eq = equal_split(N)
        ts = tangent_split(N)
        ideal = None
        if classification.lower() == "semiprime":
            flat = []
            for p, e in result.factors.items():
                flat.extend([p] * e)
            if len(flat) == 2:
                ideal = ideal_split(N, flat[0], flat[1])
        print_tangent(eq, ts, ideal)


@main.command()
@click.argument("n", type=str)
@click.option("--budget", default=10000, show_default=True, help="Time budget in ms")
@click.option(
    "--method",
    type=click.Choice(["auto", "trial", "rho", "pm1", "ecm", "qs"]),
    default="auto",
    show_default=True,
)
@click.option("--verbose", is_flag=True, help="Show pipeline steps")
@click.option("--json", "as_json", is_flag=True, help="Output machine-readable JSON")
def factor(n: str, budget: int, method: str, verbose: bool, as_json: bool):
    """
    Factor N using the full pipeline or a specific method.

    Use --verbose to show pipeline steps.
    Use --json for machine-readable output.
    """
    from .display.output import print_factor

    try:
        N = int(n.strip().replace(",", "").replace("_", ""))
    except ValueError:
        console.print("[red]Error:[/red] N must be an integer.")
        sys.exit(1)

    t0 = time.monotonic()

    if method == "auto":
        from .core.factor import factor as do_factor
        result = do_factor(N, budget_ms=budget)

    elif method == "trial":
        from .core.primes import _SMALL_PRIMES
        from .core.factor import FactorResult

        factors = {}
        tmp = N
        for p in _SMALL_PRIMES:
            while tmp % p == 0:
                factors[p] = factors.get(p, 0) + 1
                tmp //= p
        if tmp > 1:
            factors[tmp] = factors.get(tmp, 0) + 1

        result = FactorResult(
            n=N,
            factors=factors,
            method="trial",
            elapsed_ms=(time.monotonic() - t0) * 1000,
            complete=True,
        )

    elif method == "rho":
        from .core.rho import pollard_rho
        from .core.factor import FactorResult

        f = pollard_rho(N, timeout_ms=budget)
        if f:
            result = FactorResult(
                n=N,
                factors={f: 1, N // f: 1},
                method="rho",
                elapsed_ms=(time.monotonic() - t0) * 1000,
                complete=True,
            )
        else:
            result = FactorResult(
                n=N,
                factors={N: 1},
                method="rho",
                elapsed_ms=(time.monotonic() - t0) * 1000,
                complete=False,
            )

    elif method == "pm1":
        from .core.pm1 import pollard_pm1
        from .core.factor import FactorResult

        f = pollard_pm1(N)
        if f:
            result = FactorResult(
                n=N,
                factors={f: 1, N // f: 1},
                method="p-1",
                elapsed_ms=(time.monotonic() - t0) * 1000,
                complete=True,
            )
        else:
            result = FactorResult(
                n=N,
                factors={N: 1},
                method="p-1",
                elapsed_ms=(time.monotonic() - t0) * 1000,
                complete=False,
            )

    elif method == "ecm":
        from .core.ecm import ecm as do_ecm
        from .core.factor import FactorResult

        f = do_ecm(N, timeout_ms=budget)
        if f:
            result = FactorResult(
                n=N,
                factors={f: 1, N // f: 1},
                method="ecm",
                elapsed_ms=(time.monotonic() - t0) * 1000,
                complete=True,
            )
        else:
            result = FactorResult(
                n=N,
                factors={N: 1},
                method="ecm",
                elapsed_ms=(time.monotonic() - t0) * 1000,
                complete=False,
            )

    elif method == "qs":
        from .core.qs import quadratic_sieve
        from .core.factor import FactorResult

        f = quadratic_sieve(N)
        if f:
            result = FactorResult(
                n=N,
                factors={f: 1, N // f: 1},
                method="qs",
                elapsed_ms=(time.monotonic() - t0) * 1000,
                complete=True,
            )
        else:
            result = FactorResult(
                n=N,
                factors={N: 1},
                method="qs",
                elapsed_ms=(time.monotonic() - t0) * 1000,
                complete=False,
            )

    if as_json:
        payload = build_json_result(result, command="factor")
        print_json(payload)
        return

    print_factor(result, verbose=verbose)


@main.command()
@click.argument("n", type=str)
@click.option("--signature", is_flag=True, help="Show SHA-256 signatures")
@click.option("--ascii", "ascii_mode", is_flag=True, help="Show ASCII helix visualization")
@click.option("--r0", default=1.0, show_default=True)
@click.option("--alpha", default=0.0125, show_default=True)
@click.option("--beta", default=0.005, show_default=True)
@click.option("--L", default=360.0, show_default=True)
def coil(
    n: str,
    signature: bool,
    ascii_mode: bool,
    r0: float,
    alpha: float,
    beta: float,
    l: float,
):
    """Show conical helix footprint for N (must be semiprime or provide factors)."""
    from .core.factor import classify as do_classify
    from .geometry.coil import coil_footprint
    from .display.output import print_coil

    try:
        N = int(n.strip().replace(",", "").replace("_", ""))
    except ValueError:
        console.print("[red]Error:[/red] N must be an integer.")
        sys.exit(1)

    classification, result = do_classify(N)

    if classification.lower() == "prime":
        console.print(f"[green]{N} is prime[/green] — no coil footprint (requires semiprime).")
        return

    primes_flat = []
    for p, e in result.factors.items():
        primes_flat.extend([p] * e)
    primes_flat.sort()

    if len(primes_flat) < 2:
        console.print("[yellow]Could not fully factor N within budget.[/yellow]")
        return

    p, q = primes_flat[0], primes_flat[-1]
    fp = coil_footprint(N, p, q, r0=r0, alpha=alpha, beta=beta, L=l)

    if ascii_mode:
        from .display.ascii_helix import print_ascii_helix
        print_ascii_helix(fp)
        return

    print_coil(fp, show_signature=signature)


@main.command()
@click.argument("n", type=str)
@click.option("--density", is_flag=True, help="Show prime density table for nearby buckets")
@click.option("--max-k", default=20, show_default=True, help="Max bit-length for density table")
def bitbucket(n: str, density: bool, max_k: int):
    """Show bit-bucket placement and prime density analysis for N."""
    from .geometry.bitbucket import bit_bucket, prime_density_table
    from .display.output import print_bitbucket

    try:
        N = int(n.strip().replace(",", "").replace("_", ""))
    except ValueError:
        console.print("[red]Error:[/red] N must be an integer.")
        sys.exit(1)

    bb = bit_bucket(N)
    dt = None
    if density:
        effective_max = min(max_k, bb.k + 2)
        dt = prime_density_table(max_k=effective_max)
    print_bitbucket(bb, density_table=dt)


@main.command()
@click.option("--start", default=2, show_default=True, type=int)
@click.option("--stop", required=True, type=int)
@click.option("--out", default="scan.csv", show_default=True)
@click.option("--mode", type=click.Choice(["thin", "full", "sampled"]), default="thin")
@click.option("--no-resume", is_flag=True, help="Start fresh, don't resume")
@click.option("--budget", default=2000, show_default=True, help="Per-number budget in ms")
def scan(start: int, stop: int, out: str, mode: str, no_resume: bool, budget: int):
    """Wheel-accelerated range scan from START to STOP, writing to CSV."""
    from .scan.wheel import scan as do_scan

    count = 0
    for row in do_scan(start, stop, out, mode=mode, resume=not no_resume, budget_ms=budget):
        count += 1
        if count % 10000 == 0:
            console.print(f"[dim]  scanned {count:,} numbers, last n={row['n']:,}[/dim]")

    console.print(f"[green]Done.[/green] {count:,} numbers written to {out}")


@main.command()
@click.argument("n", type=str)
@click.option("--B1", default=50000, show_default=True, type=int)
@click.option("--curves", default=100, show_default=True, type=int)
@click.option("--budget", default=10000, show_default=True, type=int, help="Time budget in ms")
def ecm(n: str, b1: int, curves: int, budget: int):
    """Lenstra ECM — elliptic curve factoring."""
    from .core.ecm import ecm as do_ecm
    from .display.output import print_ecm

    try:
        N = int(n.strip().replace(",", "").replace("_", ""))
    except ValueError:
        console.print("[red]Error:[/red] N must be an integer.")
        sys.exit(1)

    t0 = time.monotonic()
    f = do_ecm(N, B1=b1, curves=curves, timeout_ms=budget)
    elapsed = (time.monotonic() - t0) * 1000
    print_ecm(N, f, elapsed, b1, curves)


@main.command()
@click.argument("n", type=str)
@click.option("--B-scale", default=5.0, show_default=True, type=float, help="Factor-base size multiplier")
def qs(n: str, b_scale: float):
    """Quadratic Sieve — for hard semiprimes with no small factors."""
    from .core.qs import quadratic_sieve
    from .display.output import print_qs

    try:
        N = int(n.strip().replace(",", "").replace("_", ""))
    except ValueError:
        console.print("[red]Error:[/red] N must be an integer.")
        sys.exit(1)

    console.print(f"[dim]Running Quadratic Sieve on {len(str(N))}-digit number…[/dim]")
    t0 = time.monotonic()
    f = quadratic_sieve(N, B_scale=b_scale)
    elapsed = (time.monotonic() - t0) * 1000
    print_qs(N, f, elapsed)


if __name__ == "__main__":
    main()
