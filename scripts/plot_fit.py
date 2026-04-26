import numpy as np
import matplotlib.pyplot as plt


# Data
N = np.array([1e4, 1e5, 1e6, 1e7])
P = np.array([59.28, 60.77, 62.47, 63.12]) / 100.0


# Your fitted parameters
C = 0.6133
alpha = 0.1844


def model(N):
    return 1 - C / (np.log(N) ** alpha)


# Smooth curve
N_smooth = np.logspace(4, 8, 200)
P_smooth = model(N_smooth)


plt.figure(figsize=(8, 5))

plt.plot(N_smooth, 100 * P_smooth, label="Fitted model")
plt.scatter(N, 100 * P, label="Data points", zorder=5)

plt.xscale("log")
plt.xlabel("N (log scale)")
plt.ylabel("Lopsided (%)")
plt.title("Fit vs Data")

plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("fit_vs_data.png", dpi=200)

print("Saved: fit_vs_data.png")
