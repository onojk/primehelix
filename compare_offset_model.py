import numpy as np
import matplotlib.pyplot as plt

N = np.array([1e4, 1e5, 1e6, 1e7])
P = np.array([59.28, 60.77, 62.47, 63.12]) / 100.0
delta = 1 - P

# Fit simple offset model: delta ≈ a - b / log(N)
logN = np.log(N)
coeffs = np.polyfit(1 / logN, delta, 1)

a = coeffs[1]
b = -coeffs[0]

print(f"a ≈ {a:.4f}, b ≈ {b:.4f}")

model = a - b / logN

plt.figure(figsize=(8, 5))
plt.plot(N, delta, "o-", label="Observed Δ(N)")
plt.plot(N, model, "--", label="a - b/log(N) fit")

plt.xscale("log")
plt.yscale("log")

plt.xlabel("N")
plt.ylabel("Δ(N)")
plt.title("Offset decay model")

plt.legend()
plt.grid(True, which="both", alpha=0.3)

plt.tight_layout()
plt.savefig("offset_model.png", dpi=200)

print("Saved: offset_model.png")
