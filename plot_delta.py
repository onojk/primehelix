import numpy as np
import matplotlib.pyplot as plt

# Data
N = np.array([1e4, 1e5, 1e6, 1e7])
P = np.array([59.28, 60.77, 62.47, 63.12]) / 100.0

delta = 1 - P

plt.figure(figsize=(8, 5))

plt.plot(N, delta, marker="o")

plt.xscale("log")
plt.yscale("log")

plt.xlabel("N (log scale)")
plt.ylabel("Delta = 1 - P(N)")
plt.title("Decay of non-lopsided semiprimes")

plt.grid(True, which="both", alpha=0.3)

plt.tight_layout()
plt.savefig("delta_plot.png", dpi=200)

print("Saved: delta_plot.png")
