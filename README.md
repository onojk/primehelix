# primehelix

**A unified console toolkit for prime number theory, integer factorization, and geometric number exploration.**

`primehelix` consolidates five research projects into one coherent program — combining fast factoring algorithms, cryptographic primality testing, a geometric/conical-helix model of number structure, bit-bucket analysis, and a wheel-accelerated range scanner.

---

## Quick look

```
$ primehelix classify 2147483646

╭──────────────────────────── classify ────────────────────────────╮
│  n               │ 2147483646                                     │
│  bits            │ 31                                             │
│  digits          │ 10                                             │
│  classification  │ COMPOSITE                                      │
│  factorization   │ 2 × 3^2 × 7 × 11 × 31 × 151 × 331            │
│  method          │ rho                                            │
│  time            │ 0.4 ms                                         │
╰────────────────────────────────────────────────────────────────────╯
```

```
$ primehelix factor 2147483646 --verbose

╭──────────────────────────── factor ──────────────────────────────╮
│  n              │ 2147483646                                      │
│  bits           │ 31                                              │
│  factorization  │ 2 × 3^2 × 7 × 11 × 31 × 151 × 331             │
│  method         │ rho                                             │
│  time           │ 0.4 ms                                          │
│  complete       │ yes                                             │
╰────────────────────────────────────────────────────────────────────╯
Pipeline steps:
  · trial: 2
  · trial: 3
  · trial: 3
  · trial: 7
  · trial: 11
  · trial: 31
  · rho: 151
```

```
$ primehelix classify 65537

  classification  │ PRIME
  method          │ bpsw
```

```
$ primehelix factor 2147483646 --json
{"n": 2147483646, "factors": {"2": 1, "3": 2, "7": 1, "11": 1, "31": 1, "151": 1, "331": 1}, "method": "rho", "elapsed_ms": 0.363, "complete": true, "steps": []}
```

---

## Factoring pipeline

Numbers are attacked in escalating order of cost — each stage is only reached if the previous one doesn't return a factor:

| Stage | Algorithm | Good for |
|-------|-----------|----------|
| 1 | **Trial division** | Any factor ≤ 149 |
| 2 | **Pollard p−1** | Factors with smooth p−1 |
| 3 | **Williams p+1** | Factors with smooth p+1 |
| 4 | **Pollard Rho (Brent)** | Mid-range composites |
| 5 | **Lenstra ECM** | Factors up to ~60 digits |
| 6 | **Quadratic Sieve** | Hard semiprimes with large factors |

Each stage respects a time budget; the pipeline escalates automatically.

### Primality testing

Uses **Baillie–PSW** (Miller–Rabin base-2 + strong Lucas PRP) — deterministic for all 64-bit integers, no known counterexamples beyond. Any unresolved cofactor is primality-tested before the factorization is marked complete.

---

## Geometric model

### Conical helix
Integers are mapped to points on a **conical helix** in 3D:

```
r(n) = r₀ + α·n      (radius grows with n)
θ(n) = 2π·n / L      (angular position)
z(n) = β·n           (vertical position)
```

For a semiprime `n = p·q`, the arc distances `d(n→q)`, `d(q→p)`, `d(p→1)` form a normalized "footprint" that characterizes the factor pair geometrically — balanced RSA-like primes produce a distinct signature from lopsided ones.

### Bit-bucket analysis
Integers are grouped by bit-length into buckets `[2^(k-1), 2^k − 1]`. Within each bucket, the offset `n − 2^(k-1)` normalizes position for cross-scale comparison. Prime density per bucket decays as `≈ 1/(k·ln2)`, matching the Prime Number Theorem in bit-length terms.

---

## Install

```bash
git clone https://github.com/onojk/primehelix.git
cd primehelix
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

**Requirements:** Python 3.10+, gmpy2

On Linux, install GMP first:
```bash
sudo apt install libgmp-dev libmpfr-dev libmpc-dev
```

### Sync and test

```bash
cd ~/primehelix && git pull && source .venv/bin/activate && pytest tests/ -v
```

---

## Commands

### classify

```bash
primehelix classify 91                          # semiprime: 7 × 13
primehelix classify 2147483647                  # prime (Mersenne M31)
primehelix classify 110000479000513 --coil      # + conical helix footprint
primehelix classify 65535 --tangent             # + tangent-split diagnostics
primehelix classify 2147483646 --json           # machine-readable output
```

### factor

```bash
primehelix factor 2147483646                    # full pipeline, clean output
primehelix factor 2147483646 --verbose          # show pipeline steps
primehelix factor 2147483646 --json             # JSON output
primehelix factor 2147483646 --json --verbose   # JSON with steps array
primehelix factor 8633 --method qs              # force Quadratic Sieve
primehelix factor <n> --budget 30000            # extend time budget (ms)
```

### bitbucket

```bash
primehelix bitbucket 97
primehelix bitbucket 65537 --density
```

### coil

```bash
primehelix coil 91
primehelix coil 110000479000513 --signature
```

### scan

```bash
primehelix scan --stop 100000 --out primes.csv
primehelix scan --start 1000000 --stop 2000000 --mode sampled --out scan.csv.gz
```

### ecm

```bash
primehelix ecm 104729314187 --B1 50000 --curves 100
```

### qs

```bash
primehelix qs 8633
primehelix qs 9804659461513846513
```

---

## Architecture

```
primehelix/
├── cli.py              — Click entry point (7 commands)
├── core/
│   ├── primes.py       — BPSW (Miller-Rabin + strong Lucas PRP)
│   ├── factor.py       — Pipeline orchestration and recursion
│   ├── rho.py          — Pollard Rho (Brent variant, batch-GCD)
│   ├── pm1.py          — Pollard p−1 / Williams p+1
│   ├── ecm.py          — Lenstra ECM (pure Python + gmpy2)
│   └── qs.py           — Quadratic Sieve (GF(2) left nullspace)
├── geometry/
│   ├── coil.py         — Conical helix footprint + SHA-256 signatures
│   ├── bitbucket.py    — Bit-bucket placement and density tables
│   └── tangent.py      — Equal-split / tangent-split / ideal-split
├── scan/
│   └── wheel.py        — Mod-210 wheel scanner, resumable gzip CSV
└── display/
    └── output.py       — Rich terminal panels and tables
```

---

## Origins

`primehelix` consolidates five research repositories:

| Repo | Contribution |
|------|-------------|
| `geom_factor` | Quadratic Sieve, bit-bucket theory, geometric model |
| `rsacrack` | Factoring pipeline, coil classifier |
| `ECC-Tools` | ECM reference (C + libecm) |
| `Cprime` | GMP-backed C CLI (trial + p−1 + Rho) |
| `onojk123` | Wheel scanner, tangent prime test |

---

## License

MIT
