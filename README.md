# primehelix

**primehelix** is a command-line toolkit for exploring the *structure of integers*.

Most tools stop at factorization.  
**primehelix goes further: it reveals how numbers are built.**

---

## ✨ Features

- Prime / semiprime / composite classification
- Multi-stage factoring pipeline  
  *(trial → p-1 → p+1 → rho → ECM → QS)*
- Residue-based structure labeling (modular patterns)
- Geometric interpretation (conical helix / “coil”)
- Structural comparison across integer ranges
- Top-change detection (`--top-delta`)
- JSON, CSV, JSONL outputs for automation
- Plot generation (matplotlib)
- Fully regression-tested core

---

## 🚀 Installation

```bash
git clone https://github.com/onojk/primehelix.git
cd primehelix
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
⚡ Quick Start
Classify a number
primehelix classify 1300039 --coil --residue
Scan a range
primehelix structure-scan --start 1 --stop 100
Generate a plot
primehelix structure-scan --start 1 --stop 100 --plot structures.png
Compare two ranges
primehelix compare-ranges \
  --a-start 1 --a-stop 100 \
  --b-start 100 --b-stop 200
Focus on biggest changes
primehelix compare-ranges \
  --a-start 1 --a-stop 100 \
  --b-start 100 --b-stop 200 \
  --top-delta 5
🧠 Core Concepts
1. Classification

Each integer is classified as:

prime
semiprime
composite
2. Residue Structure

primehelix labels numbers using modular patterns.

Examples:

prime | gaussian
prime | pythagorean
semiprime | moderate | mod4_1x3
composite
3. Geometric Structure (Coil)

For semiprimes:

n = p × q

Mapped to a conical helix:

balance → |p − q| / √n
bit-gap → difference in bit-length
structural symmetry
4. Structural Comparison

Compare how number structure evolves across ranges:

primehelix compare-ranges --a-start 10 --a-stop 20 --b-start 20 --b-stop 30

Output includes:

counts
percentages
delta (percentage-point change)
ratio (relative change)
emergence / disappearance of structures
5. Top Structural Changes
--top-delta N

Shows only the most significant structural shifts.

📊 Output Modes
Mode	Description
default	rich CLI tables
--json	structured output
--jsonl	streaming records
--csv	spreadsheet export
--plot	PNG visualizations
🧪 Testing
pytest tests/ -v

Includes:

factoring correctness
classification
residue structure
JSON schema validation
structure scanning
compare-ranges analytics
📁 Project Structure
primehelix/
  core/        factoring + primality
  geometry/    coil + residue systems
  display/     CLI + plotting
  scan/        range analysis
  cli.py       entrypoint
📌 Version

v0.1.2

Includes:

compare-ranges analytics
ratio + delta metrics
top-delta filtering
plotting improvements
36 passing tests
📜 License

MIT
