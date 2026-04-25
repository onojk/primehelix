# PrimeHelix

[![PyPI version](https://img.shields.io/pypi/v/primehelix.svg)](https://pypi.org/project/primehelix/)
[![Python versions](https://img.shields.io/pypi/pyversions/primehelix.svg)](https://pypi.org/project/primehelix/)
[![CI](https://github.com/onojk/primehelix/actions/workflows/ci.yml/badge.svg)](https://github.com/onojk/primehelix/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## What is PrimeHelix

PrimeHelix is a Python CLI tool and research platform for studying semiprime structure through the normalized logarithmic coordinate θ = log p / log n, where n = p·q is a semiprime with p ≤ q.

Every integer is classified by factor balance and residue family and assigned a compact structural label:

```
semiprime | lopsided | mod4_1x3
semiprime | balanced | mod4_3x3
prime | gaussian
composite
```

Five CLI commands:

| Command | What it does |
|---------|-------------|
| `classify` | Inspect one integer — label, factors, residue family, ASCII helix |
| `scan` | Count structure labels across a range |
| `compare` | Diff structure distributions between two ranges |
| `trend` | Track structural shifts over sliding windows |
| `structure-scan` | Full label breakdown with entropy and method profile |

---

## Paper

**"On the Distribution of Smaller Factors in Semiprimes"**
Jonathan Kendall, 2026

[paper/paper.pdf](paper/paper.pdf)

An empirical study of the coordinate θ = log p / log n across semiprimes n ≤ 10⁸,
computed by exact enumeration using a smallest-prime-factor sieve. The proportion
P(N, θ₀) of semiprimes satisfying θ ≤ θ₀ increases monotonically with N and is
consistent with convergence to a limiting value, but the rate is extremely slow —
consistent with O(1/log log N) decay as predicted by the log log N accumulation of
small prime contributions (Montgomery–Vaughan, §7.4). Over N ∈ [10⁴, 10⁸],
log log N varies only from 2.22 to 2.91, making a local linear fit in log N appear
accurate within this window. Geometric and spectral embeddings (additive phase
wrapping, multiplicative phase maps, spectral fingerprint, factor-ratio spectrum)
all reduce to log-density effects or small-prime artifacts — θ is not recoverable
from these natural embeddings.

Math Stack Exchange discussion: https://math.stackexchange.com/questions/5134337

---

## Key Findings

- **Lopsided dominance**: ~73% of semiprimes at N = 10⁶, rising to ~82% at N = 10⁸.
  The bias strengthens with scale.

- **Convergence rate**: O(1/log log N), consistent with Montgomery–Vaughan §7.4.
  The deficit δ(N, θ₀) = 1 − P(N, θ₀) decreases by roughly 5–7% over four
  orders of magnitude in N.

- **Natural coordinate**: θ = log p / log n is scale-invariant and directly encodes
  factor asymmetry. It is the most stable descriptor of semiprime structure.

- **Exact identity**: log(p/q) = (2θ − 1) log n holds for every semiprime. Given θ
  and log n, the factor ratio is completely determined with no residual freedom.

- **Spectral null results**: Additive phase wrapping, multiplicative phase maps n^(it),
  spectral fingerprint C(t) = |Σ n^(it)| / N, and factor-ratio spectrum
  C_ratio(t) = |Σ e^(it·log(p/q))| / N all reduce to log-density effects or discrete
  artifacts from small primes. θ is not recoverable from any of these embeddings.

- **b(θ) quantification**: The empirical convergence rate b(θ₀) varies smoothly with θ₀
  and is characterized across four thresholds (θ₀ ∈ {0.20, 0.25, 0.30, 0.35}) by the
  regression b(θ₀) ≈ k(1 − θ₀)^α, k ≈ 0.015, α ≈ 3.1 (R² = 0.97).

---

## Visualizations

Two native 3D visualizers using pyglet and OpenGL:

**`primehelix_viz.py`** — Helix visualization with balanced semiprimes on the central axis.

**`primehelix_lattice.py`** — Lattice graph with family strands, axis spokes, and cross chords.

```bash
python3 primehelix_lattice.py --n 2000 --radius 80
```

Controls:

| Input | Action |
|-------|--------|
| Drag | Rotate |
| Scroll wheel | Zoom (proportional) |
| Space | Toggle auto-rotate |
| `1` / `2` / `3` | Toggle strand / spoke / chord edges |

Color coding: purple = 1×1 family, teal = 1×3 family, orange-red = 3×3 family,
gray = even-involved, amber = balanced semiprimes (axis).

---

## Datasets

Pre-computed scans are in [`datasets/`](datasets/):

| File | Range | Records |
|------|-------|---------|
| `semiprime_1e6_labels.csv` | [1, 10⁶) | 210,035 semiprimes |
| `semiprime_1e7_labels.csv` | [1, 10⁷) | 1,904,324 semiprimes |
| `semiprime_1e8_labels.csv` | [1, 10⁸) | 17,427,258 semiprimes |
| `compare_semiprime_1e6_halves.csv` | [1, 500k) vs [500k, 1M) | per-label deltas |
| `compare_semiprime_1e7_halves.csv` | [1, 5M) vs [5M, 10M) | per-label deltas |

Each file includes the CLI command used to produce it. 1B scan in progress.

---

## Installation

```bash
pip install primehelix
```

On Linux, install GMP first for full performance:

```bash
sudo apt install libgmp-dev libmpfr-dev libmpc-dev
pip install primehelix
```

---

## Research Scripts

Scripts in the repository root, produced during the empirical study:

| Script | Purpose |
|--------|---------|
| `verify_loglogN_rate.py` | Fits Model A (linear in log N) and Model B (c/log log N) to the deficit; confirms O(1/log log N) theoretical rate |
| `exact_lopsided_counts.py` | Exact enumeration of lopsided semiprimes across N ∈ {10⁴ … 10⁸} |
| `dickman_compare.py` | Compares empirical convergence to Dickman-function predictions |
| `helix_factor_prediction_test.py` | Three experiments testing whether helix angle carries factor information — null result |
| `multiplicative_phase_scan.py` | Phase clustering by small factor p at generic and zeta-zero t values |
| `spectral_fingerprint.py` | C(t) = \|Σ n^(it)\| / N for semiprimes, primes, and random integers; marks first 10 Riemann zeta zeros |
| `factor_ratio_spectral.py` | Three spectral functions over the log(p/q) domain; autocorrelation and power spectrum |
| `fit_b_theta_model.py` | Regression fits for b(θ₀) scaling law |

---

## Geometry

The cylinder model maps semiprime structure to 3D geometry:

- **Lopsided semiprimes** wind around the cylinder as the helix. Moving up the cylinder corresponds to increasing N.
- **Balanced semiprimes** cluster near the central axis — the spine the helix winds around.
- **Family strands** connect semiprimes in the same mod4 residue family (1×1, 1×3, 3×3, even-involved).
- **Spokes** connect semiprimes at the same N value across families.
- **Cross chords** link semiprimes sharing a small prime factor.

The geometry is a visualization tool. It does not constitute a factoring method or encode
any structure beyond what is already in θ and the residue family labels.

---

## License

MIT

---

## Author

Jonathan Kendall — https://github.com/onojk
