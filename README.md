# On the Distribution of Smaller Factors in Semiprimes

Research code and reproducibility scripts for the paper of the same name.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## Requirements

- Python ≥ 3.9
- numpy, matplotlib, scipy

```bash
pip install numpy matplotlib scipy
```

RAM requirements:
- **~400 MB** for the SPF sieve to N = 10⁸ (`verify_semiprime_tables.py`)
- **~1 GB** for the pair-counting extension to N = 10⁹ (`extend_conjecture1.py`)

---

## Reproduce all results

Run the canonical verification script (~8 minutes, ~400 MB RAM):

```bash
python3 verify_semiprime_tables.py
```

This runs all 6 tasks and 3 checkpoint validations. Every task should print `PASS`.

To extend Conjecture 1 to N = 10⁹ and fit subleading correction models
(~1 GB RAM, ~30–60 minutes depending on hardware):

```bash
python3 extend_conjecture1.py
```

---

## Outputs

### `verify_semiprime_tables.py`

| File | Description |
|------|-------------|
| `results/task1_table1.csv` | Table 1: P(N, θ₀) values across θ₀ ∈ {0.20, 0.25, 0.30, 0.35} |
| `results/task2_table2.csv` | Table 2: δ(N, θ₀) and b(θ₀) linear regression |
| `results/task3_table3.csv` | Table 3: b_eff / f(θ₀) stability ratio |
| `results/task5_extended_table2.csv` | Extended Table 2 (denser θ₀ grid) |
| `results/task6_stabilisation.csv` | Stabilisation data (b_eff vs N) |
| `figures/activation_plot.png` | Activation threshold figure |

### `extend_conjecture1.py`

| File | Description |
|------|-------------|
| `results/conjecture1_extended.csv` | b_eff at N ∈ {10⁸, 2×10⁸, 5×10⁸, 10⁹} for θ₀ ∈ [0.10, 0.45] |
| `results/subleading_models.csv` | Six subleading model fits (S1–S6): A, R², max abs residual |
| `stabilisation_extended.png` | b_eff(N, θ₀) extended to N = 10⁹ |

### Checkpoint validation

`verify_semiprime_tables.py` validates three exact semiprime counts before any computation:

| N | Expected count |
|---|---------------|
| 10⁶ | 210,035 |
| 10⁷ | 1,904,324 |
| 10⁸ | 17,427,258 |

All three must match exactly. A mismatch aborts with an error.

---

## Paper

**"On the Distribution of Smaller Factors in Semiprimes"**
Jonathan Kendall, 2026

[paper/paper.pdf](paper/paper.pdf) · [paper/paper.tex](paper/paper.tex)

Math Stack Exchange discussion: https://math.stackexchange.com/questions/5134337

An empirical study of the coordinate θ = log p / log n across semiprimes n ≤ 10⁹,
computed by exact enumeration. The proportion P(N, θ₀) of semiprimes satisfying
θ ≤ θ₀ increases monotonically with N, consistent with convergence at O(1/log log N)
as predicted by Montgomery–Vaughan §7.4. The effective convergence rate
b_eff = δ · log log N approaches the Chebyshev-corrected heuristic f(θ₀) = log(1/θ₀) − θ₀ − log 2 + ½
to within a factor of 1.21 across all tested θ₀. Subleading structure (six models tested,
best R² = 0.25) remains open.

---

## Reproducibility note

`verify_paper_tables.py` is the original reference script used during paper development.
`verify_semiprime_tables.py` is the canonical reproducibility script and should be preferred
for all verification. Both produce identical checkpoint counts and PASS/FAIL outcomes.

All outputs are fully deterministic: no random seeds, no sampling. Results depend only on
exact semiprime enumeration via smallest-prime-factor sieve.

---

## PrimeHelix package

This repository also contains the `primehelix` Python package — a CLI tool for studying
semiprime structure via the θ coordinate and residue families.

```bash
pip install primehelix
```

Five commands: `classify`, `scan`, `compare`, `trend`, `structure-scan`.

See `pyproject.toml` for full package details.

---

## Author

Jonathan Kendall — https://github.com/onojk
