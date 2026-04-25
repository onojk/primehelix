# Submission checklist — Integers journal

## Compile
- [ ] cd paper && pdflatex paper.tex: zero errors
- [ ] cd paper && pdflatex paper.tex (twice): zero undefined references
- [ ] All figures render (no missing file warnings)
- [ ] cover_letter.tex compiles clean

## Code verification
- [ ] python3 verify_semiprime_tables.py runs to completion
- [ ] Table 1: all 20 values match paper to < 0.005 pp
- [ ] Table 2: R2 > 0.97 for theta0 in {0.20, 0.25, 0.30, 0.35}
- [ ] Table 3: max beff/f <= 1.21 (verified: 1.2015)
- [ ] Checksums pass: 210035 / 1904324 / 17427258

## Claims correctly scoped
- [ ] P_inf(theta0) = 1 stated as unproven assumption
- [ ] Proposition (jump size) marked heuristic
- [ ] Boundary-layer section in appendix, marked speculative
- [ ] f(theta0) validity range [0.10, 0.45] stated explicitly
- [ ] Indicator decomposition paragraph present in Section 3.3
- [ ] Subleading correction addressed in Section 5

## Repository
- [ ] README.md present and accurate
- [ ] paper/drafts/ contains old versions (not paper.tex)
- [ ] github.com/onojk/primehelix is public
- [ ] All scripts committed and pushed
- [ ] results/ CSVs committed

## Still pending (do not claim in paper)
- [ ] N=1e9 stabilisation data (sieve still running)
- [ ] Extended Conjecture 1 plot with N=1e9 points
