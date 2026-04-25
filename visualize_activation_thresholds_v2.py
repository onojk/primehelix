#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

theta0 = 0.30

N_vals = np.array([1e4, 1e5, 1e6, 1e7, 1e8])
P_vals = np.array([67.73, 69.40, 70.51, 71.83, 72.74])

primes = [17, 19, 23, 29]
label_y = [72.9, 72.55, 72.25, 71.9]


def threshold(p):
    return p ** (1.0 / theta0)


fig, ax = plt.subplots(figsize=(9.0, 5.2))

ax.plot(
    N_vals,
    P_vals,
    marker="o",
    linewidth=2.4,
    markersize=6.5,
    color="blue",
    label=r"$\theta_0 = 0.30$",
    zorder=3,
)

ax.set_xscale("log")
ax.set_xlabel(r"$N$", fontsize=12)
ax.set_ylabel(r"$P(N,0.30)$ (%)", fontsize=12)
ax.set_title(r"Activation thresholds for $\theta_0 = 0.30$", fontsize=14, pad=12)

ax.set_ylim(67.35, 73.10)

ax.grid(True, which="major", alpha=0.28)
ax.grid(True, which="minor", linestyle=":", alpha=0.18)

for i, p in enumerate(primes):
    x = threshold(p)

    ax.axvline(
        x,
        color="gray",
        linestyle="--",
        linewidth=1.1,
        alpha=0.55,
        zorder=1,
    )

    ax.text(
        x,
        label_y[i],
        rf"$p={p}$",
        fontsize=9,
        color="darkred",
        ha="center",
        va="center",
        zorder=4,
    )

ax.legend(loc="lower right", fontsize=10)

fig.tight_layout()

fig.savefig("figures/activation_theta_0.30.png", dpi=300, bbox_inches="tight")
fig.savefig("figures/activation_threshold_theta030.png", dpi=300, bbox_inches="tight")

print("Saved:")
print("  figures/activation_theta_0.30.png")
print("  figures/activation_threshold_theta030.png")
