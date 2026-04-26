import numpy as np
import pandas as pd
from pathlib import Path

# Load data
df = pd.read_csv("semiprime_data.csv")

N_vals = df["N"].to_numpy()
loglogN = np.log(np.log(N_vals))

theta_main = [0.20, 0.25, 0.30, 0.35]
theta_pref = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45]

# 🔴 THIS WAS MISSING
out = []

# =========================
# Table 1
# =========================
out.append(r"""
\begin{table}[htbp]
\centering
\caption{Proportion \(P(N,\theta_0)\) of semiprimes \(n \leq N\) satisfying
\(\theta \leq \theta_0\). All values are percentages.}
\label{tab:convergence}
\begin{tabular}{rcccc}
\toprule
\(N\) & \(\theta_0 = 0.20\) & \(\theta_0 = 0.25\) & \(\theta_0 = 0.30\) & \(\theta_0 = 0.35\) \\
\midrule
""")

for _, row in df.iterrows():
    out.append(
        rf"\(10^{int(np.log10(row['N']))}\) & "
        rf"{row['P_0.20']:.2f} & {row['P_0.25']:.2f} & "
        rf"{row['P_0.30']:.2f} & {row['P_0.35']:.2f} \\"
    )

out.append(r"""
\bottomrule
\end{tabular}
\end{table}
""")

# =========================
# Table 2 (NumPy regression)
# =========================
out.append(r"""
\begin{table}[htbp]
\centering
\caption{Parameters of the local linear fit
\(\delta \approx a(\theta_0)-b(\theta_0)\log N\) over \(N\in[10^4,10^8]\).}
\label{tab:fit}
\begin{tabular}{cccc}
\toprule
\(\theta_0\) & \(a(\theta_0)\) & \(b(\theta_0)\) & \(R^2\) \\
\midrule
""")

X = np.log(N_vals)

for theta in theta_main:
    P = df[f"P_{theta:.2f}"].to_numpy() / 100.0
    delta = 1.0 - P

    A = np.vstack([X, np.ones(len(X))]).T
    slope, a = np.linalg.lstsq(A, delta, rcond=None)[0]

    b = -slope

    delta_pred = a + slope * X
    ss_res = np.sum((delta - delta_pred) ** 2)
    ss_tot = np.sum((delta - np.mean(delta)) ** 2)
    r2 = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0

    out.append(rf"{theta:.2f} & {a:.4f} & {b:.6f} & {r2:.4f} \\")

out.append(r"""
\bottomrule
\end{tabular}
\end{table}
""")

# =========================
# Table 3 (prefactor)
# =========================
out.append(r"""
\begin{table}[htbp]
\centering
\caption{Empirical \(b_{\mathrm{eff}}\) vs. heuristic prefactors.}
\label{tab:prefactor}
\resizebox{\textwidth}{!}{%
\begin{tabular}{cccccc}
\toprule
\(\theta_0\) & \(\log(1/\theta_0)\) & \(f(\theta_0)\) & \(b_{\mathrm{eff}}\)
& \(b_{\mathrm{eff}}/\log(1/\theta_0)\) & \(b_{\mathrm{eff}}/f(\theta_0)\) \\
\midrule
""")

for theta in theta_pref:
    col = f"P_{theta:.2f}"
    if col not in df.columns:
        continue

    P = df[col].to_numpy() / 100.0
    delta = 1.0 - P

    beff = np.mean(delta * loglogN)

    naive = np.log(1.0 / theta)
    f = naive - theta - np.log(2.0) + 0.5

    out.append(
        rf"{theta:.2f} & {naive:.3f} & {f:.3f} & {beff:.3f} & "
        rf"{beff/naive:.3f} & {beff/f:.3f} \\"
    )

out.append(r"""
\bottomrule
\end{tabular}%
}
\end{table}
""")

# Save file
Path("generated_tables.tex").write_text("\n".join(out))

print("✅ generated_tables.tex written successfully")
