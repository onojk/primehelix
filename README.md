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
  "structure": "semiprime | mod4_1x3",
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

| Field | Always present | Notes |
|-------|---------------|-------|
| `command` | yes | `"classify"` or `"factor"` |
| `n` | yes | integer |
| `classification` | classify only | `"prime"`, `"semiprime"`, `"composite"`, `"invalid"` |
| `factors` | yes | `{"p": exponent, ...}` |
| `factorization` | yes | `"2 * 3^2 * 7"` (ASCII) |
| `method` | yes | last algorithm used |
| `complete` | yes | `true` if all factors proven prime |
| `structure` | classify only | compact label string |
| `steps` | with `--verbose` | pipeline step trail |
| `coil` | with `--coil` | geometric footprint + insight |
| `residue` | classify only | mod4/mod6/mod30 profile |

Breaking changes to this schema will be documented in release notes and accompanied by a minor version bump.

---

## Install and test

```bash
git clone https://github.com/onojk/primehelix.git
cd primehelix
python3 -m venv .venv
source .venv/bin/activate
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
