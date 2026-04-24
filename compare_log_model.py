import numpy as np
import matplotlib.pyplot as plt

N = np.array([1e4, 1e5, 1e6, 1e7])
P = np.array([59.28, 60.77, 62.47, 63.12]) / 100.0
delta = 1 - P

log_model = 1 / np.log(N)

plt.figure(figsize=(8, 5))

plt.plot(N, delta, marker="o", label="Observed Δ(N)")
plt.plot(N, log_model, linestyle="--", label="1 / log(N)")

plt.xscale("log")
plt.yscale("log")

plt.xlabel("N")
plt.ylabel("Δ(N)")
plt.title("Compare decay to 1/log(N)")

plt.legend()
plt.grid(True, which="both", alpha=0.3)

plt.tight_layout()
plt.savefig("compare_log_model.png", dpi=200)

print("Saved: compare_log_model.png")
