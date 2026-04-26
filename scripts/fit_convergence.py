import numpy as np
from scipy.optimize import curve_fit


# Data
N = np.array([1e4, 1e5, 1e6, 1e7])
P = np.array([59.28, 60.77, 62.47, 63.12]) / 100.0


def model(N, C, alpha):
    return 1 - C / (np.log(N) ** alpha)


params, _ = curve_fit(model, N, P, p0=[1.0, 1.0])

C, alpha = params

print("\n=== FIT RESULTS ===")
print(f"C = {C:.4f}")
print(f"alpha = {alpha:.4f}")
