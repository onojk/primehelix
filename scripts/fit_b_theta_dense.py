import numpy as np
import matplotlib.pyplot as plt
import os
import random

# Deterministic
os.environ["MPLBACKEND"] = "Agg"
random.seed(0)
np.random.seed(0)


# ============================================
# INPUT: extended theta grid + your measured b
# (you will update these if needed later)
# ============================================

theta = np.array([
    0.15, 0.20, 0.25, 0.30, 0.35, 0.40
])

# You MUST fill this with real values once computed
# For now we interpolate your known data (temporary)
theta_known = np.array([0.20, 0.25, 0.30, 0.35])
b_known = np.array([0.008712, 0.007253, 0.005672, 0.005250])

b = np.interp(theta, theta_known, b_known)


def r_squared(y_true, y_pred):
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    return 1.0 - ss_res / ss_tot


# ============================================
# MODEL FITS
# ============================================

# Model A: b = k (1 - θ)^α
x1 = np.log(1.0 - theta)
y = np.log(b)

alpha, log_k = np.polyfit(x1, y, 1)
k = np.exp(log_k)

b_fit_A = k * (1.0 - theta) ** alpha
r2_A = r_squared(b, b_fit_A)


# Model B: b = c / θ^γ
x2 = np.log(theta)
m2, log_c = np.polyfit(x2, y, 1)
gamma = -m2
c = np.exp(log_c)

b_fit_B = c / (theta ** gamma)
r2_B = r_squared(b, b_fit_B)


# ============================================
# DIAGNOSTIC TESTS (CRUCIAL)
# ============================================

theta_b = theta * b


# ============================================
# PRINT RESULTS
# ============================================

print("\n=== DENSE θ FIT ===\n")

print("Model A: b = k(1-θ)^α")
print(f"k     = {k:.8f}")
print(f"alpha = {alpha:.6f}")
print(f"R²    = {r2_A:.6f}")
print()

print("Model B: b = c / θ^γ")
print(f"c     = {c:.8f}")
print(f"gamma = {gamma:.6f}")
print(f"R²    = {r2_B:.6f}")
print()

print("θ * b(θ) values:")
for t, val in zip(theta, theta_b):
    print(f"θ={t:.2f} → θ*b={val:.6f}")


# ============================================
# PLOTS
# ============================================

theta_smooth = np.linspace(0.14, 0.42, 400)

b_smooth_A = k * (1.0 - theta_smooth) ** alpha
b_smooth_B = c / (theta_smooth ** gamma)


# --- Fit comparison
plt.figure(figsize=(8, 5))
plt.scatter(theta, b, label="Observed b(θ)", zorder=5)
plt.plot(theta_smooth, b_smooth_A, label="k(1-θ)^α")
plt.plot(theta_smooth, b_smooth_B, linestyle="--", label="c/θ^γ")
plt.xlabel("θ")
plt.ylabel("b(θ)")
plt.title("Dense θ fit for b(θ)")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("b_theta_dense_fit.png", dpi=200)


# --- θ * b(θ) test
plt.figure(figsize=(8, 5))
plt.plot(theta, theta_b, marker="o")
plt.xlabel("θ")
plt.ylabel("θ * b(θ)")
plt.title("Test: θ * b(θ) (should be flat if ~1/θ law)")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("theta_times_b.png", dpi=200)


# --- log-log test
plt.figure(figsize=(8, 5))
plt.plot(np.log(theta), np.log(b), marker="o")
plt.xlabel("log(θ)")
plt.ylabel("log(b)")
plt.title("Log-log test: b vs θ")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("loglog_b_theta.png", dpi=200)


print("\nSaved:")
print(" - b_theta_dense_fit.png")
print(" - theta_times_b.png")
print(" - loglog_b_theta.png")
