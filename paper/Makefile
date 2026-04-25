# Makefile for compiling paper2.tex (semiprime factor asymmetry paper)
# Requires: pdflatex, python (for generating figures if needed)

.PHONY: all clean

# Primary target
all: paper2.pdf

# Generate the figures from the data (if plots.py or semiprime_data.csv change)
# This rule assumes plots.py exists and uses semiprime_data.csv to regenerate all PNGs
activation_theta_0.30.png theta_vs_logratio.png stabilisation_test.png: plots.py semiprime_data.csv
	python plots.py

# Compile the LaTeX document (run twice for cross-references and TOC)
paper2.pdf: paper2.tex activation_theta_0.30.png theta_vs_logratio.png stabilisation_test.png
	pdflatex paper2.tex
	pdflatex paper2.tex

# Remove all auxiliary files and the generated PDF
clean:
	rm -f *.aux *.log *.out *.toc *.pdf
