import matplotlib.pyplot as plt

from plot_multi_theta import N, compute_curve


THETAS = [0.20, 0.25, 0.30, 0.35]


def main():
    plt.figure(figsize=(9, 5))

    for theta in THETAS:
        p_values = compute_curve(theta)
        delta_values = [100.0 - p for p in p_values]

        plt.plot(
            N,
            delta_values,
            marker="o",
            linewidth=2,
            label=f"θ = {theta:.2f}",
        )

    plt.xscale("log")

    plt.xlabel("N cutoff (log scale)")
    plt.ylabel("Balanced / non-lopsided semiprimes (%)")
    plt.title("Normalized Δ(N, θ) = 1 - P(N, θ)")

    plt.grid(True, alpha=0.3)
    plt.legend()

    plt.tight_layout()

    output_file = "multi_theta_delta.png"
    plt.savefig(output_file, dpi=200)

    print(f"Saved: {output_file}")


if __name__ == "__main__":
    main()
