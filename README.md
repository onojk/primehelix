# primehelix

[![PyPI version](https://img.shields.io/pypi/v/primehelix.svg)](https://pypi.org/project/primehelix/)
[![Python versions](https://img.shields.io/pypi/pyversions/primehelix.svg)](https://pypi.org/project/primehelix/)
[![CI](https://github.com/onojk/primehelix/actions/workflows/ci.yml/badge.svg)](https://github.com/onojk/primehelix/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Structural analysis for integers — explore how numbers are built, not just what they are.**

Most number theory tools answer *what*: is this prime, what are the factors. primehelix also answers *how*: what arithmetic family does each factor belong to, how balanced is the factor pair, where does the number sit on a conical helix, and how do these structural patterns shift across large ranges.

Every integer gets a compact **structure label** — `semiprime | lopsided | mod4_1x3`, `prime | gaussian`, `composite` — that encodes classification, geometric balance, and residue-family membership into one scannable token. Those labels are the spine of the tool: classify one number, scan a million, compare two ranges, plot trends over time.

---

## Install

```bash
pip install primehelix                # core: classify, factor, scan, compare
pip install 'primehelix[plot]'        # add matplotlib for --plot
```

On Linux, GMP is required for full performance (gmpy2):
```bash
sudo apt install libgmp-dev libmpfr-dev libmpc-dev
pip install primehelix
```

```bash
primehelix classify 1300039 --helix
primehelix classify 1300039 --json
primehelix structure-scan --start 1 --stop 100000
primehelix compare-ranges --a-start 1 --a-stop 50000 --b-start 50000 --b-stop 100000 --top-delta 6
```

**From source:**
```bash
git clone https://github.com/onojk/primehelix.git
cd primehelix
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

---

## Core concepts

### Structure labels

Every integer is assigned a structure label composed of up to three parts joined by ` | `:

```
semiprime | lopsided | mod4_1x3
prime | gaussian
composite
invalid
```

| Part | What it encodes |
|------|----------------|
| Classification | `prime`, `semiprime`, `composite`, `invalid` |
| Balance | `balanced`, `moderate`, `lopsided` — based on bit-length gap between factors |
| Residue family | `mod4_1x3`, `mod4_3x3`, `pythagorean`, `gaussian`, etc. |

These labels are stable strings — safe to grep, aggregate, diff between ranges, and track over time.

### Residue families

Odd primes split into two families by their residue mod 4:

- **Pythagorean primes** (p ≡ 1 mod 4) — expressible as a sum of two squares
- **Gaussian primes** (p ≡ 3 mod 4) — remain prime in the Gaussian integers

Semiprimes inherit a **mod4 pair** from their two factors: `1x1`, `1x3`, or `3x3`. This pair is stable under the prime number theorem — its distribution converges predictably as range grows, and shifts between ranges are measurable.

### Conical helix geometry

Integers are mapped to points on a conical helix in 3D:

```
r(n) = r₀ + α·n      radius grows with n
θ(n) = 2π·n / L      angular position
z(n) = β·n           vertical position
```

For a semiprime `n = p·q`, the arc distances between `n`, `q`, `p`, and `1` form a geometric footprint. The **bit-gap** between the factors controls how the helix spreads — balanced RSA-like primes produce a symmetric shape; lopsided pairs compress one strand. The `--helix` flag renders this as an ASCII double-helix in the terminal.

---

## Commands

### `classify` — classify and inspect one integer

```bash
primehelix classify 1300039
primehelix classify 1300039 --helix       # ASCII helix visualization
primehelix classify 1300039 --coil        # geometric footprint metrics
primehelix classify 1300039 --residue     # full residue profile
primehelix classify 1300039 --json        # machine-readable output
```

**`--helix` output** (1300039 = 13 × 100003, bit_gap=13):

```
1300039 → semiprime

Helix (p=13, q=100003)
balance=87.696, bit_gap=13

                      +-------------------*
                     +                     *
                     *---------------------+
                        *               +
                            +~~~~~~~*
                                +
                            +~~~~~~~*
                          +           *
                       *-----------------+
                    *                       +
                +-------------------------------*
              +                                   *
               *---------------------------------+
```

The spread and compression of the helix reflect the actual factor structure — a balanced semiprime like 110000479000513 (= 10000019 × 11000027, bit_gap=0) produces a tight symmetric pattern; a lopsided one like 1300039 produces a wide expanding cone.

**`--json` output:**

```json
{
  "command": "classify",
  "n": 1300039,
  "classification": "semiprime",
  "factors": {"13": 1, "100003": 1},
  "factorization": "13 * 100003",
  "method": "trial",
  "complete": true,
  "structure": "semiprime | lopsided | mod4_1x3",
  "residue": {
    "semiprime_mod4_pair": "1x3",
    "semiprime_mod4_note": "mixed 1 mod 4 and 3 mod 4 factor families",
    "factor_families_mod4": ["pythagorean", "gaussian"]
  }
}
```

---

### `factor` — full factoring pipeline

```bash
primehelix factor 2147483646
primehelix factor 2147483646 --verbose    # show pipeline steps
primehelix factor 2147483646 --json
primehelix factor 2147483646 --json --verbose
```

**Pipeline:** trial division → Pollard p−1 → Williams p+1 → Pollard Rho (Brent) → Lenstra ECM → Quadratic Sieve

**Output:**
```
  n              │ 2147483646
  factorization  │ 2 × 3^2 × 7 × 11 × 31 × 151 × 331
  method         │ rho
  complete       │ yes

Pipeline steps (--verbose):
  · trial: 2  · trial: 3  · trial: 3  · trial: 7
  · trial: 11  · trial: 31  · rho: 151
```

Primality testing uses **Baillie–PSW** (Miller–Rabin base-2 + strong Lucas PRP) — deterministic for all 64-bit integers. Prime cofactors are always proven before the factorization is marked complete.

---

### `structure-scan` — count structure labels across a range

```bash
primehelix structure-scan --start 1 --stop 1000000
primehelix structure-scan --start 1 --stop 1000000 --only-classification semiprime
primehelix structure-scan --start 1 --stop 1000000 --json
```

Scans every integer in `[start, stop)`, assigns a structure label, and returns counts with a histogram. Progress is shown on stderr for ranges over 10,000 numbers.

**Output (abridged):**

```
              structure summary
 ┌────────────────────────────────┬────────┬────────┬──────────────────────┐
 │ structure                      │  count │percent │ histogram            │
 ├────────────────────────────────┼────────┼────────┼──────────────────────┤
 │ composite                      │ 531820 │ 53.18% │ ██████████████████████│
 │ prime | gaussian               │  87432 │  8.74% │ ████████             │
 │ prime | pythagorean            │  80251 │  8.03% │ ███████              │
 │ semiprime | moderate | mod4_1x3│  93114 │  9.31% │ ████████             │
 │ semiprime | lopsided | mod4_1x3│  ...   │  ...   │ ...                  │
 └────────────────────────────────┴────────┴────────┴──────────────────────┘
```

---

### `compare-ranges` — diff structure distributions between two ranges

```bash
primehelix compare-ranges \
  --a-start 1 --a-stop 50000 \
  --b-start 50000 --b-stop 100000 \
  --top-delta 6
```

Shows which structure types grew or shrank most between two ranges, with counts, percentages, delta, and ratio.

**Output:**

```
             range comparison | top delta 6
 ┌─────────────────────────────────┬───────┬───────┬───────┬───────┬───────┬───────┐
 │ structure                       │ [1,50k│  [1,  │[50k,  │ [50k, │ delta │ ratio │
 │                                 │ count │  50k% │100k)  │100k)% │       │       │
 ├─────────────────────────────────┼───────┼───────┼───────┼───────┼───────┼───────┤
 │ composite                       │ 32755 │ 65.5% │ 34273 │ 68.5% │ +3.0% │ 1.05x │
 │ semiprime | moderate | mod4_1x3 │  2203 │  4.4% │  1711 │  3.4% │ -0.9% │ 0.78x │
 │ prime | gaussian                │  2583 │  5.2% │  2225 │  4.5% │ -0.7% │ 0.86x │
 │ semiprime | lopsided | mod4_1x3 │  2376 │  4.8% │  2684 │  5.4% │ +0.6% │ 1.13x │
 └─────────────────────────────────┴───────┴───────┴───────┴───────┴───────┴───────┘
```

Use `--only-classification semiprime` to isolate one class. Use `--json` to pipe results downstream.

---

### `structure-time-series` — track structural trends over sliding windows

```bash
primehelix structure-time-series \
  --start 1 --stop 1000000 \
  --window 100000 --step 100000 \
  --only-classification semiprime \
  --top 5 \
  --plot semiprime_ts.png
```

Divides `[start, stop)` into overlapping or non-overlapping windows, computes structure distributions in each, selects the top-N series by aggregate weight, and plots them as a line chart. Omit `--plot` for a compact text summary instead.

---

## JSON output

`classify` and `factor` both support `--json`. The schema is stable across patch versions:

| Field | Present in | Notes |
|-------|-----------|-------|
| `command` | both | `"classify"` or `"factor"` |
| `n` | both | integer |
| `classification` | classify | `"prime"`, `"semiprime"`, `"composite"`, `"invalid"` |
| `factors` | both | `{"p": exponent, ...}` |
| `prime_factors` | both | flat list, e.g. `[3, 3, 7]` for 3²×7 |
| `factorization` | both | `"2 * 3^2 * 7"` (ASCII) |
| `method` | both | last algorithm used |
| `elapsed_ms` | both | wall time in milliseconds |
| `complete` | both | `true` if all factors proven prime |
| `structure` | classify | compact label — `"semiprime \| lopsided \| mod4_1x3"` |
| `steps` | factor with `--verbose` | pipeline step trail; empty list otherwise |
| `coil` | classify with `--coil` | geometric footprint + insight string |
| `residue` | classify | mod4/mod6/mod30 profile |

The `structure-scan` and `compare-ranges` commands also include:

| Field | Command | Notes |
|-------|---------|-------|
| `entropy` | structure-scan | Shannon entropy (bits) of the label distribution. 0 = single label, log₂(k) = uniform over k labels |
| `a.entropy`, `b.entropy` | compare-ranges | Entropy of each range independently |
| `entropy_delta` | compare-ranges | `b.entropy − a.entropy`; positive = B more diverse |

Breaking changes to this schema will be documented in release notes and accompanied by a minor version bump.

---

## Guarantees and limits

**Deterministic:**
- Structure labels and residue families are computed from factorization alone — identical input always produces identical output.
- Primality testing uses Baillie–PSW (Miller-Rabin base-2 + strong Lucas PRP), which is deterministic for all integers up to 2⁶⁴. No known counterexamples exist.
- `complete: true` means every factor has been proven prime. The factorization is exact.

**May time out:**
- The factoring pipeline has a configurable budget (`--budget`, default 10 000 ms). For numbers with large prime factors that resist trial division and Pollard Rho, the pipeline may exhaust its budget and return `complete: false` with a partial factorization.
- For most integers up to ~15 digits, factorization completes in milliseconds. Harder numbers (e.g. RSA-like products of two large primes) may time out.

**Stable and scriptable:**
- `classify`, `structure-scan`, `compare-ranges`, `structure-time-series` with `--json` produce stable, machine-readable output safe to pipe, grep, and aggregate.
- Structure labels are stable strings — they are designed to be safe keys for counting and comparison across runs.

**Experimental:**
- `--coil` and `--helix` output (geometric footprint, ASCII visualization) reflects a model under active development. The coordinate values and balance thresholds may change between minor versions.
- The insight strings in `coil.insight` are heuristic and human-readable only — do not parse them programmatically.

---

## Empirical findings

All measurements below were produced by running primehelix against [1, 1 000 000). The commands are fully reproducible.

### Overall structure distribution

```bash
primehelix structure-scan --start 1 --stop 1000000 --json
```

| Classification | Count | Share |
|----------------|------:|------:|
| composite | 711,465 | 71.15% |
| semiprime | 210,035 | 21.00% |
| prime | 78,498 | 7.85% |
| invalid (n ≤ 1) | 1 | — |

Primes split almost exactly evenly between the two residue families — 50.09% gaussian (p ≡ 3 mod 4) and 49.91% pythagorean (p ≡ 1 mod 4) — consistent with Dirichlet's theorem on primes in arithmetic progressions.

### Semiprime balance distribution

Among the 210,035 semiprimes in [1, 1M):

| Balance tier | Count | Share |
|--------------|------:|------:|
| lopsided (bit_gap > 8 or balance ≥ 10) | 153,718 | **73.2%** |
| moderate | 54,427 | 25.9% |
| balanced (bit_gap ≤ 1 and balance < 0.15) | 1,677 | 0.80% |

Lopsided pairs dominate by a wide margin. Balanced semiprimes — the RSA-like products of two primes of nearly equal bit-length — are extremely rare below 1M: under 1 in 125.

### Mod4 pair distribution: all semiprimes vs lopsided-only

| Mod4 pair | All semiprimes | Lopsided only | Shift |
|-----------|---------------:|-------------:|------:|
| mod4_1x3 (mixed families) | 40.0% | 36.4% | −3.6 pp |
| mod4_3x3 (both gaussian) | 23.7% | 22.9% | −0.9 pp |
| mod4_1x1 (both pythagorean) | 16.4% | 13.7% | −2.7 pp |
| even-involved (factor of 2) | 19.8% | **27.0%** | **+7.2 pp** |

The lopsided constraint systematically shifts the distribution toward even-involved pairs. The explanation is structural: any semiprime of the form 2×p is always lopsided (p is at least 2 bits larger than 2 for p ≥ 5), so the entire even semiprime population is absorbed into the lopsided bucket. Mixed (1x3) and symmetric (1x1, 3x3) pairs are all proportionally reduced.

### Lopsidedness grows with range

```bash
primehelix compare-ranges \
  --a-start 1 --a-stop 500000 \
  --b-start 500000 --b-stop 1000000 \
  --only-classification semiprime --top-delta 6 --json
```

| Structure | delta | ratio |
|-----------|------:|------:|
| semiprime \| lopsided \| mod4_1x3 | +2.69% | 1.11× |
| semiprime \| moderate \| mod4_1x3 | −2.05% | 0.85× |
| semiprime \| lopsided \| mod4_3x3 | +1.51% | 1.09× |
| semiprime \| moderate \| mod4_3x3 | −1.39% | 0.81× |
| semiprime \| lopsided \| mod4_1x1 | +1.29% | 1.14× |
| semiprime \| moderate \| mod4_1x1 | −0.66% | 0.90× |

As the range shifts from [1, 500k) to [500k, 1M), lopsided semiprimes gain share and moderate ones shrink — uniformly across all three odd mod4 families. The mechanism: small primes (2, 3, 5, 7, …) are repeatedly reused as the smaller factor in semiprimes that reach into higher ranges, producing an ever-wider bit-gap between the two factors.

---

## Develop and test

```bash
git clone https://github.com/onojk/primehelix.git
cd primehelix
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pytest tests/ -v
```

**Sync after updates:**
```bash
cd ~/primehelix && git pull && source .venv/bin/activate && pytest tests/ -v
```

---

## Architecture

```
primehelix/
├── cli.py                  — 5 Click commands + scan helpers
├── core/
│   ├── primes.py           — Baillie-PSW (Miller-Rabin + strong Lucas PRP)
│   ├── factor.py           — Pipeline orchestration
│   ├── rho.py              — Pollard Rho (Brent, batch-GCD)
│   ├── pm1.py              — Pollard p−1 / Williams p+1
│   ├── ecm.py              — Lenstra ECM (pure Python + gmpy2)
│   └── qs.py               — Quadratic Sieve (GF(2) left nullspace)
├── geometry/
│   ├── coil.py             — Conical helix model, CoilFootprint, CoilBalance
│   ├── residue.py          — Mod4/mod6/mod30 residue profiling
│   ├── bitbucket.py        — Bit-bucket placement and density
│   └── tangent.py          — Equal/tangent/ideal split diagnostics
├── display/
│   ├── output.py           — Rich terminal panels and tables
│   ├── json_output.py      — JSON schema, structure_summary label builder
│   ├── plots.py            — Matplotlib time-series line charts
│   └── ascii_helix.py      — ASCII double-helix renderer
└── scan/
    └── wheel.py            — Mod-210 wheel scanner, resumable gzip CSV
```

---

## Origins

primehelix consolidates five research repositories:

| Repo | Contribution |
|------|-------------|
| `geom_factor` | Quadratic Sieve, bit-bucket theory, geometric model |
| `rsacrack` | Factoring pipeline, coil classifier |
| `ECC-Tools` | ECM reference (C + libecm) |
| `Cprime` | GMP-backed C CLI (trial + p−1 + Rho) |
| `onojk123` | Wheel scanner, tangent prime test |

---

## Author

Jonathan Kendall
https://github.com/onojk

---

## License

MIT
