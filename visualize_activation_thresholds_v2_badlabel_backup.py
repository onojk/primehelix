#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

theta0 = 0.30

N_vals = np.array([1e4, 1e5, 1e6, 1e7, 1e8])
P_vals = np.array([67.73, 69.40, 70.51, 71.83, 72.74])

# ONLY show a few clean primes (this is key)
PRIMES = [17, 19, 23, 29]

def threshold(p):
    return p ** (1.0 / theta0)

fig, ax = plt.subplots(figsize=(9, 5))

# main curve
ax.plot(
    N_vals,
    P_vals,
    marker='o',
    linewidth=2,
    color='blue',
    label=r'$\theta_0 = 0.30$'
)

ax.set_xscale("log")
ax.set_xlabel(r"$N$")
ax.set_ylabel(r"$P(N,0.30)$ (%)")

# better spacing
ymin = min(P_vals) - 0.3
ymax = max(P_vals) + 0.6
ax.set_ylim(ymin, ymax)

# clean grid
ax.grid(True, which="major", alpha=0.3)
ax.grid(True, which="minor", linestyle=":", alpha=0.2)

# --- KEY FIX: clean labeling ---
for i, p in enumerate(PRIMES):
    t = threshold(p)

    if t < N_vals[0] or t > N_vals[-1]:
        continue

    ax.axvline(t, linestyle='--', color='gray', alpha=0.5)

    # place label ABOVE plot, not inside clutter
    ax.annotate(
        f"$p={p}$",
        xy=(t, ymax),
        xytext=(0, 6),
        textcoords="offset points",
        ha='center',
        va='bottom',
        fontsize=9,
        color='darkred'
    )

# title with padding (prevents collision)
ax.set_title(
    r"Activation thresholds for $\theta_0 = 0.30$",
    pad=20
)

ax.legend(loc="lower right")

plt.tight_layout()

plt.savefig(
    "figures/activation_theta_0.30.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()
