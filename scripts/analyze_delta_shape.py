import numpy as np
import matplotlib.pyplot as plt

N = np.array([1e4, 1e5, 1e6, 1e7, 1e8])
P = np.array([59.28, 60.77, 62.47, 63.74, 64.95]) / 100.0

delta = 1 - P
logN = np.log(N)

plt.figure(figsize=(8, 5))

# raw delta
plt.plot(logN, delta, marker="o", label="Δ(N) vs log(N)")

# try linear fit in log-space
coeffs = np.polyfit(logN, delta, 1)
fit = np.polyval(coeffs, logN)

plt.plot(logN, fit, linestyle="--", label="Linear fit in log(N)")

plt.xlabel("log(N)")
plt.ylabel("Δ(N)")
plt.title("Shape of decay vs log(N)")
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("delta_vs_logN.png", dpi=200)

print("Saved: delta_vs_logN.png")
print(f"Slope ≈ {coeffs[0]:.6f}, intercept ≈ {coeffs[1]:.6f}")
