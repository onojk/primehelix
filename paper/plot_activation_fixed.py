import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("semiprime_data.csv")
N_vals = df["N"]

theta0 = 0.30
col_name = f"P_{theta0:.2f}"

if col_name not in df.columns:
    raise ValueError(f"Column {col_name} not found. Available: {df.columns.tolist()}")

P_vals = df[col_name]

fig, ax = plt.subplots(figsize=(8, 5))

ax.semilogx(
    N_vals,
    P_vals,
    "o-",
    color="blue",
    linewidth=2,
    markersize=8,
    label=fr"$\theta_0 = {theta0:.2f}$",
)

ax.set_xlabel(r"$N$", fontsize=12)
ax.set_ylabel(fr"$P(N, {theta0:.2f})$ (%)", fontsize=12)
ax.set_title(fr"Activation thresholds for $\theta_0 = {theta0:.2f}$", fontsize=14)
ax.grid(True, which="both", linestyle=":", alpha=0.5)

xmin, xmax = ax.get_xlim()

# Keep only readable thresholds
primes = [17, 23, 31, 41]

for p in primes:
    x_thresh = p ** (1.0 / theta0)

    if xmin <= x_thresh <= xmax:
        ax.axvline(
            x=x_thresh,
            linestyle="--",
            color="gray",
            alpha=0.65,
            linewidth=1,
        )

        ax.text(
            x_thresh,
            -0.13,
            f"p={p}",
            transform=ax.get_xaxis_transform(),
            rotation=90,
            ha="center",
            va="top",
            fontsize=8,
            color="darkred",
            clip_on=False,
        )

ax.text(
    0.02,
    0.96,
    r"Dashed lines: $N = p^{1/\theta_0}$",
    transform=ax.transAxes,
    va="top",
    fontsize=9,
)

ax.legend(loc="lower right")

fig.subplots_adjust(bottom=0.20)
fig.tight_layout()

outname = f"activation_theta_{theta0:.2f}_fixed.png"
fig.savefig(outname, dpi=300, bbox_inches="tight")
print(f"Figure saved as {outname}")
