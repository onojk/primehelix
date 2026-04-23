# PrimeHelix

PrimeHelix is a computational and visualization toolkit for exploring the structural patterns of integers—especially primes, semiprimes, and composite numbers—through classification, residue analysis, and time-series behavior.

The goal is to move beyond simple counting and uncover emergent structure in number distributions.

---

## Features

- Fast integer classification:
  - Prime
  - Semiprime
  - Composite
- Structural labeling:
  - Balanced vs lopsided factorization
- Residue analysis:
  - mod 4 pair structures (e.g., 1×3, 3×3, etc.)
- Time-series analysis:
  - Sliding window structure distribution
- Visualization:
  - Plotting of structural trends over ranges
  - ASCII helix representation (experimental)

---

## Installation

Clone the repo and install dependencies:

git clone https://github.com/onojk/primehelix.git
cd primehelix
python3 -m venv .venv
source .venv/bin/activate
pip install -e .

---

## Example Usage

### Structure Time Series (Semiprimes)

primehelix structure-time-series \
  --start 1 --stop 100000 \
  --window 10000 --step 10000 \
  --only-classification semiprime \
  --only-structure lopsided \
  --plot lopsided_ts.png

---

### Compare Two Ranges

primehelix compare-ranges \
  --a-start 1 --a-stop 50000 \
  --b-start 50000 --b-stop 100000 \
  --only-classification semiprime \
  --top-delta 15

---

## Key Observations

- mod4_1×3 becomes dominant (~35%) as range increases  
- mod4_3×3 stabilizes (~22%)  
- Structures involving even factors decline over scale  
- mod4_1×1 increases gradually  

This suggests convergence toward a stable residue distribution in semiprimes.

---

## Development

Run tests:

pytest tests/ -v

---

## Author

Jonathan Kendall  
https://github.com/onojk
