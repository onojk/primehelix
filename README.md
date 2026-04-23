# primehelix

**Structural analysis for integers — explore how numbers are built, not just what they are.**

Most number theory tools answer *what*: is this prime, what are the factors. primehelix also answers *how*: what arithmetic family does each factor belong to, how balanced is the factor pair, where does the number sit on a conical helix, and how do these structural patterns shift across large ranges.

Every integer gets a compact **structure label** — `semiprime | lopsided | mod4_1x3`, `prime | gaussian`, `composite` — that encodes classification, geometric balance, and residue-family membership into one scannable token. Those labels are the spine of the tool: classify one number, scan a million, compare two ranges, plot trends over time.

---

## Quick start

```bash
git clone https://github.com/onojk/primehelix.git
cd primehelix
python3 -m venv .venv && source .venv/bin/activate
pip install -e .
```

```bash
primehelix classify 1300039 --helix
primehelix classify 1300039 --json
primehelix structure-scan --start 1 --stop 100000
primehelix compare-ranges --a-start 1 --a-stop 50000 --b-start 50000 --b-stop 100000 --top-delta 6
```

On Linux, install GMP before `pip install`:
```bash
sudo apt install libgmp-dev libmpfr-dev libmpc-dev
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

The following measurements come from running primehelix against the integers [1, 100 000).

### Structure distribution baseline

```
primehelix structure-scan --start 1 --stop 100000 --json
```

| Structure | Count | Share |
|-----------|------:|------:|
| composite | 67 028 | 67.03% |
| semiprime \| lopsided \| mod4_1x3 | 5 060 | 5.06% |
| prime \| gaussian | 4 808 | 4.81% |
| prime \| pythagorean | 4 783 | 4.78% |
| semiprime \| moderate \| mod4_1x3 | 3 914 | 3.91% |
| semiprime \| lopsided \| mod4_3x3 | 3 349 | 3.35% |
| semiprime \| lopsided \| mod4_2x3 | 2 559 | 2.56% |
| semiprime \| lopsided \| mod4_2x1 | 2 528 | 2.53% |
| semiprime \| moderate \| mod4_3x3 | 2 112 | 2.11% |
| semiprime \| moderate \| mod4_1x1 | 1 798 | 1.80% |
| semiprime \| balanced (all pairs) | 217 | 0.22% |

**Notes:** Primes split nearly evenly between the two residue families (gaussian ≈ pythagorean), consistent with Dirichlet's theorem. Among semiprimes, lopsided pairs (larger factor at least 8 bits wider) are more frequent than moderate pairs by roughly 2.5×. Balanced semiprimes — pairs where both factors have nearly equal bit-length — are rare: under 0.25% of all integers in this range.

### Lopsidedness grows with range

```
primehelix compare-ranges \
  --a-start 1 --a-stop 50000 \
  --b-start 50000 --b-stop 100000 \
  --only-classification semiprime --top-delta 8 --json
```

| Structure | delta | ratio |
|-----------|------:|------:|
| semiprime \| lopsided \| mod4_1x3 | +4.20% | 1.21× |
| semiprime \| moderate \| mod4_1x3 | −3.01% | 0.83× |
| semiprime \| lopsided \| mod4_3x3 | +2.38% | 1.18× |
| semiprime \| moderate \| mod4_3x3 | −2.31% | 0.77× |

As the integer range shifts from [1, 50k) to [50k, 100k), lopsided semiprimes gain share at the expense of moderate ones — consistently across all mod4 families. This reflects the growing gap between small primes (reused as the smaller factor) and the larger prime cofactors needed to reach higher products.

---

## Install and test

```bash
git clone https://github.com/onojk/primehelix.git
cd primehelix
python3 -m venv .venv
source .venv/bin/activate
pip install -e .           # core: classify, factor, scan, compare
pip install -e ".[plot]"   # add matplotlib for --plot
pip install -e ".[dev]"    # everything including tests
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
