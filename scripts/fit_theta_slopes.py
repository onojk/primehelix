import numpy as np
import matplotlib.pyplot as plt

from plot_multi_theta import N, compute_curve


THETAS = [0.20, 0.25, 0.30, 0.35]


def main():
    log_n = np.log(np.array(N, dtype=float))

    theta_values = []
    slopes = []
    intercepts = []

    print("\n=== FIT Δ(N, θ) = a(θ) - b(θ) log(N) ===")
    print(f"{'theta':>8} {'a(theta)':>12} {'b(theta)':>12} {'slope':>12}")
    print("-" * 48)

    for theta in THETAS:
        p_values = np.array(compute_curve(theta), dtype=float) / 100.0
        delta = 1.0 - p_values

        # Linear fit: delta = m*log(N) + c
        m, c = np.polyfit(log_n, delta, 1)

        # Rewrite as delta = a - b*log(N)
        a = c
        b = -m

        theta_values.append(theta)
        slopes.append(b)
        intercepts.append(a)

        print(f"{theta:>8.2f} {a:>12.6f} {b:>12.6f} {m:>12.6f}")

    plt.figure(figsize=(8, 5))
    plt.plot(theta_values, slopes, marker="o", linewidth=2)

    plt.xlabel("θ")
    plt.ylabel("b(θ)")
    plt.title("Decay-rate parameter b(θ)")

    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    output_file = "theta_slope_fit.png"
    plt.savefig(output_file, dpi=200)

    print(f"\nSaved: {output_file}")


if __name__ == "__main__":
    main()
