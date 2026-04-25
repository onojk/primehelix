#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Usage:
#   ./glean_tenenbaum_snippets.sh [optional_pdf_path]
#
# If no argument is given, defaults to Downloads location.
# ------------------------------------------------------------

PDF="${1:-$HOME/Downloads/gsm163-endmatter.pdf}"

OUTDIR="$HOME/primehelix/tenenbaum_snippets"
TXT="$OUTDIR/tenenbaum_full.txt"
SNIPS="$OUTDIR/relevant_snippets_for_claude.txt"

mkdir -p "$OUTDIR"

# ------------------------------------------------------------
# Check PDF exists
# ------------------------------------------------------------
if [ ! -f "$PDF" ]; then
    echo "ERROR: PDF not found at:"
    echo "  $PDF"
    echo
    echo "Run like this if needed:"
    echo "  $0 ~/Downloads/gsm163-endmatter.pdf"
    exit 1
fi

# ------------------------------------------------------------
# Check pdftotext exists
# ------------------------------------------------------------
if ! command -v pdftotext >/dev/null 2>&1; then
    echo "ERROR: pdftotext not installed."
    echo "Install with:"
    echo "  sudo apt install poppler-utils"
    exit 1
fi

echo "Using PDF:"
echo "  $PDF"
echo

echo "Extracting PDF text..."
pdftotext -layout "$PDF" "$TXT"

echo "Building relevant snippets..."

cat > "$SNIPS" <<'HEADER'
Relevant excerpts from Gérald Tenenbaum,
Introduction to Analytic and Probabilistic Number Theory

Purpose:
Connect semiprime lopsidedness

    n = p q,  p <= n^theta0

to known machinery around smooth numbers, rough numbers,
Dickman function, Buchstab function, sieve theory, and prime-factor models.

Empirical observation:

    delta(N, theta0) ~ f(theta0) / log log N

but experimentally:

    f(theta0) != log(1/theta0)

and finite-N behavior shows discrete prime activation thresholds.

============================================================
HEADER

# ------------------------------------------------------------
# Extraction function
# ------------------------------------------------------------
extract_context () {
    local pattern="$1"
    local label="$2"

    echo "" >> "$SNIPS"
    echo "============================================================" >> "$SNIPS"
    echo "$label" >> "$SNIPS"
    echo "Search pattern: $pattern" >> "$SNIPS"
    echo "============================================================" >> "$SNIPS"
    echo "" >> "$SNIPS"

    grep -n -i -C 18 "$pattern" "$TXT" >> "$SNIPS" || true
}

# ------------------------------------------------------------
# Key extractions
# ------------------------------------------------------------

extract_context "Dickman" "DICKMAN FUNCTION"
extract_context "Buchstab" "BUCHSTAB FUNCTION"
extract_context "friable" "FRIABLE INTEGERS"
extract_context "smooth" "SMOOTH NUMBERS"
extract_context "saddle-point" "SADDLE POINT METHOD"
extract_context "Rankin" "RANKIN METHOD"
extract_context "Kubilius" "KUBILIUS MODEL"
extract_context "Psi(x, y)" "PSI(x,y)"
extract_context "Ψ(x, y)" "PSI UNICODE"
extract_context "Phi(x, y)" "PHI(x,y)"
extract_context "Φ(x, y)" "PHI UNICODE"
extract_context "largest prime factor" "LARGEST PRIME FACTOR"
extract_context "smallest prime factor" "SMALLEST PRIME FACTOR"
extract_context "P +(n)" "P PLUS NOTATION"
extract_context "P -(n)" "P MINUS NOTATION"
extract_context "prime factors" "PRIME FACTORS"
extract_context "Selberg" "SELBERG SIEVE"
extract_context "Brun" "BRUN SIEVE"
extract_context "sieve" "SIEVE METHODS"
extract_context "Mertens" "MERTENS THEOREM"
extract_context "Hardy" "HARDY-RAMANUJAN"
extract_context "Ramanujan" "RAMANUJAN"
extract_context "Erdos" "ERDOS"
extract_context "Erdős" "ERDOS (ACCENT)"
extract_context "Turan" "TURAN-KUBILIUS"
extract_context "Turán" "TURAN-KUBILIUS (ACCENT)"
extract_context "omega" "OMEGA FUNCTIONS"
extract_context "Ω" "BIG OMEGA"
extract_context "ω" "SMALL OMEGA"

# ------------------------------------------------------------
# Claude prompt appended
# ------------------------------------------------------------

cat >> "$SNIPS" <<'FOOTER'

============================================================
PROMPT TO FEED CLAUDE
============================================================

Using the excerpts above, analyze:

Let n = p q be a semiprime with p <= q, and define:

    p <= n^theta0

Empirical result:

    delta(N, theta0) ~ f(theta0) / log log N

but:

    f(theta0) != log(1/theta0)

Tasks:

1. Derive the continuous heuristic giving log(1/theta0)
2. Explain why discrete primes distort the prefactor
3. Analyze sum:

       sum_{p <= N^theta0} pi(N/p)

4. Explain finite-N step behavior
5. Suggest corrected prefactor model

FOOTER

echo
echo "Done."
echo
echo "Full text:"
echo "  $TXT"
echo
echo "Claude-ready snippets:"
echo "  $SNIPS"
echo
echo "Preview:"
echo "------------------------------------------------------------"
head -80 "$SNIPS"
