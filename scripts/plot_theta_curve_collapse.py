import matplotlib.pyplot as plt

from plot_multi_theta import N, compute_curve


THETAS = [0.20, 0.25, 0.30, 0.35]


def main():
    plt.figure(figsize=(9, 5))

    for theta in THETAS:
        p_values = compute_curve(theta)
        delta_values = [100.0 - p for p in p_values]

        # Collapse by removing the theta-dependent baseline
        baseline = delta_values[0]
        collapsed = [d - baseline for d in delta_values]

        plt.plot(
            N,
            collapsed,
            marker="o",
            linewidth=2,
            label=f"θ = {theta:.2f}",
        )

    plt.xscale("log")

    plt.xlabel("N cutoff (log scale)")
    plt.ylabel("Δ(N,θ) - Δ(1e4,θ)")
    plt.title("Curve Collapse Test for Δ(N,θ)")

    plt.axhline(0, linewidth=1, alpha=0.4)
    plt.grid(True, alpha=0.3)
    plt.legend()

    plt.tight_layout()

    output_file = "theta_curve_collapse.png"
    plt.savefig(output_file, dpi=200)

    print(f"Saved: {output_file}")


if __name__ == "__main__":
    main()
