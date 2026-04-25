# On the Distribution of Smaller Factors in Semiprimes
*Jonathan Kendall, 2026*

An empirical study of the coordinate θ = log p / log n for semiprimes n = pq (p ≤ q).
Exact enumeration up to N ≤ 10⁸ shows that the proportion P(N, θ₀) of semiprimes
with θ ≤ θ₀ increases monotonically with N, consistent with O(1/log log N) decay toward
a limiting value — the classical rate for related semiprime counts (Montgomery–Vaughan §7.4).
A novel discrete activation phenomenon is documented: each prime p contributes a visible
step in P(N, θ₀) when N crosses p^(1/θ₀). A refined heuristic prefactor
f(θ₀) = log(1/θ₀) − θ₀ − log 2 + ½ reduces the maximum relative prediction error
for δ · log log N from 4.3× (naive Mertens) to under 1.21× across θ₀ ∈ [0.10, 0.45].
Several geometric and spectral embeddings (additive phase wrapping, multiplicative phase
maps, spectral fingerprint, factor-ratio spectrum) all fail to recover the structure
carried by θ, consistent with a separation principle: factor asymmetry is not recoverable
from any function of n alone that is multiplicative or additive in log-space.

## Requirements
- Python 3.8+, numpy, scipy, matplotlib
- ~400 MB RAM for SPF sieve to N = 1e8
- Runtime: ~8 minutes

## Reproduce all results
```
python3 verify_semiprime_tables.py
```

`verify_semiprime_tables.py` is the canonical script (commit fc43ef0).
`verify_paper_tables.py` is an earlier version kept for reference only.

## Outputs
- results/table1.csv, table2.csv, table3.csv
- activation_theta_0_10.png
- stabilisation.png

## Paper
```
cd paper && pdflatex paper.tex && pdflatex paper.tex
```

## Verified checksums (exact, no sampling)
- 210,035 semiprimes ≤ 1e6
- 1,904,324 semiprimes ≤ 1e7
- 17,427,258 semiprimes ≤ 1e8

## Repository structure
```
paper/          — paper.tex, paper.pdf (canonical)
paper/drafts/   — earlier draft versions
results/        — CSV outputs from verification
verify_semiprime_tables.py  — canonical verification script
cover_letter.tex            — Integers journal submission
```
