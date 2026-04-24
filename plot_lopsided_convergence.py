import matplotlib.pyplot as plt


# Exact results (from your computation)
N = [10_000, 100_000, 1_000_000, 10_000_000, 100_000_000]
P = [59.28, 60.77, 62.47, 63.74, 64.95]


def main():
    plt.figure(figsize=(8, 5))

    plt.plot(
        N,
        P,
        marker="o",
        linewidth=2,
        label="Exact enumeration"
    )

    plt.xscale("log")

    plt.title("Exact Lopsided Semiprime Convergence (θ = 0.25)")
    plt.xlabel("N cutoff (log scale)")
    plt.ylabel("Lopsided semiprimes (%)")

    plt.grid(True, alpha=0.3)
    plt.legend()

    plt.tight_layout()

    output_file = "exact_lopsided_convergence.png"
    plt.savefig(output_file, dpi=200)

    print(f"Saved plot: {output_file}")


if __name__ == "__main__":
    main()
