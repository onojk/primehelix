# primehelix

**A unified console toolkit for prime number theory, integer factorization, and geometric number exploration.**

`primehelix` consolidates five research projects into one coherent program — combining fast factoring algorithms, cryptographic primality testing, a geometric/conical-helix model of number structure, bit-bucket analysis, and a wheel-accelerated range scanner.

---

## Concepts

### Factoring Pipeline
Numbers are attacked in escalating order of cost:

1. **Trial division** — instant for small factors
2. **Pollard p−1** — fast when p−1 is smooth
3. **Williams p+1** — fast when p+1 is smooth
4. **Pollard Rho (Brent)** — probabilistic, good for mid-range composites
5. **ECM (Lenstra)** — elliptic curve method, scales to ~60-digit factors
6. **Quadratic Sieve** — for hard semiprimes where both factors are large

### Primality Testing
Uses **Baillie–PSW** (Miller–Rabin base-2 + strong Lucas PRP) — deterministic for all 64-bit integers, no known counterexamples beyond.

### Geometric / Coil Model
Integers are mapped to points on a **conical helix** in 3D space:
- `r(n) = r₀ + α·n`
- `θ(n) = 2π·n / L`
- `z(n) = β·n`

For a semiprime `n = p·q`, the distances `d(n→q)`, `d(q→p)`, `d(p→1)` along the helix form a normalized "footprint" that characterizes the factor pair geometrically.

### Bit-Bucket Analysis
Integers are grouped by bit-length into buckets `[2^(k-1), 2^k − 1]`. Within each bucket, the **offset** `n − 2^(k-1)` normalizes position for cross-scale comparison. Prime density per bucket decays as `≈ 1/(k·ln2)`, matching the Prime Number Theorem in bit-length terms.

### Wheel-Accelerated Scanner
A mod-210 wheel (coprime residues mod 2·3·5·7) skips ~77% of candidates. The scanner is resumable, supports gzip output, and adapts its stride based on local density.

---

## Install

```bash
git clone https://github.com/onojk/primehelix.git
cd primehelix
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

**Requirements:** Python 3.10+, gmpy2 (requires libgmp-dev on Linux)

```bash
sudo apt install libgmp-dev libmpfr-dev libmpc-dev
pip install gmpy2
```

---

## Usage

```
primehelix <command> [options]
```

### classify
```bash
primehelix classify 91
primehelix classify 2147483647
primehelix classify 110000479000513 --coil
```

### factor
```bash
primehelix factor 15
primehelix factor 110000479000513
primehelix factor 8633 --method qs
primehelix factor <big-number> --budget 10000
```

### scan
```bash
primehelix scan --stop 100000 --out primes.csv
primehelix scan --start 1000000 --stop 2000000 --mode sampled --out scan.csv.gz
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

### qs
```bash
primehelix qs 9804659461513846513
```

### ecm
```bash
primehelix ecm 104729314187 --B1 50000 --curves 100
```

---

## Architecture

```
primehelix/
├── cli.py              — Click-based entry point
├── core/
│   ├── primes.py       — BPSW (Miller-Rabin + Lucas), wheel sieve
│   ├── factor.py       — Full pipeline orchestration
│   ├── rho.py          — Pollard Rho (Brent, batch-GCD)
│   ├── pm1.py          — Pollard p−1 / Williams p+1
│   └── qs.py           — Quadratic Sieve with GF(2) nullspace
│   └── ecm.py          — Lenstra ECM (pure Python + gmpy2)
├── geometry/
│   ├── coil.py         — Conical helix model and footprint
│   ├── bitbucket.py    — Bit-bucket grouping and density
│   └── tangent.py      — Tangent-split diagnostics
├── scan/
│   └── wheel.py        — Mod-210 wheel scanner, resumable CSV
└── display/
    └── output.py       — Rich terminal formatting
```

---

## Origins

`primehelix` consolidates five research repositories:

| Repo | Contribution |
|------|-------------|
| `geom_factor` | Quadratic sieve, bit-bucket theory, geometric visualizations |
| `rsacrack` | Factoring pipeline, coil classifier, Flask API |
| `ECC-Tools` | ECM via libecm (C reference implementation) |
| `Cprime` | GMP-backed C CLI (trial + P-1 + Rho) |
| `onojk123` | Wheel scanner, tangent prime test |

---

## License

MIT
