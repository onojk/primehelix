import numpy as np
import matplotlib.pyplot as plt
import os
import random

# Deterministic behavior
os.environ["MPLBACKEND"] = "Agg"
random.seed(0)
np.random.seed(0)


# Your measured data
theta = np.array([0.20, 0.25, 0.30, 0.35])
b = np.array([0.008712, 0.007253, 0.005672, 0.005250])


def r_squared(y_true, y_pred):
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    return 1.0 - ss_res / ss_tot


# =========================
# MODEL 1: k(1 - θ)^α
# =========================
x1 = np.log(1.0 - theta)
y = np.log(b)

alpha, log_k = np.polyfit(x1, y, 1)
k = np.exp(log_k)

b_fit_1 = k * (1.0 - theta) ** alpha
r2_1 = r_squared(b, b_fit_1)


# =========================
# MODEL 2: c / θ^γ
# =========================
x2 = np.log(theta)
m2, log_c = np.polyfit(x2, y, 1)
gamma = -m2
c = np.exp(log_c)

b_fit_2 = c / (theta ** gamma)
r2_2 = r_squared(b, b_fit_2)


# =========================
# PRINT RESULTS
# =========================
print("\n=== FIT b(θ) MODELS ===\n")

print("Model 1: b(θ) = k * (1 - θ)^α")
print(f"k     = {k:.10f}")
print(f"alpha = {alpha:.6f}")
print(f"R²    = {r2_1:.6f}")
print()

print("Model 2: b(θ) = c / θ^γ")
print(f"c     = {c:.10f}")
print(f"gamma = {gamma:.6f}")
print(f"R²    = {r2_2:.6f}")
print()


# =========================
# PLOT
# =========================
theta_smooth = np.linspace(0.18, 0.38, 300)

b_smooth_1 = k * (1.0 - theta_smooth) ** alpha
b_smooth_2 = c / (theta_smooth ** gamma)

plt.figure(figsize=(8, 5))

plt.scatter(theta, b, label="Observed b(θ)", zorder=5)
plt.plot(theta_smooth, b_smooth_1, label="k(1-θ)^α fit")
plt.plot(theta_smooth, b_smooth_2, linestyle="--", label="c/θ^γ fit")

plt.xlabel("θ")
plt.ylabel("b(θ)")
plt.title("Model fits for decay-rate parameter b(θ)")
plt.grid(True, alpha=0.3)
plt.legend()

plt.tight_layout()

output_file = "b_theta_model_fit.png"
plt.savefig(output_file, dpi=200)

print(f"Saved: {output_file}")
