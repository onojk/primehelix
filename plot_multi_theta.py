import matplotlib.pyplot as plt

# Exact results you already computed (θ = 0.25)
N = [10_000, 100_000, 1_000_000, 10_000_000, 100_000_000]
P_025 = [59.28, 60.77, 62.47, 63.74, 64.95]

# Recompute for other θ using your exact counter
# (reuses your exact_lopsided_counts logic via a small helper)
from exact_lopsided_counts import sieve
from bisect import bisect_left, bisect_right
from math import isqrt

def count_ratio(limit, theta, primes):
    total = 0
    lopsided = 0

    for p in primes:
        if p * p > limit:
            break

        q_min_total = p
        q_max = limit // p

        total += bisect_right(primes, q_max) - bisect_left(primes, q_min_total)

        # General θ condition: p <= (pq)^θ  ⇔  p^(1/θ - 1) <= q
        power = int(round(1/theta - 1))
        q_min_lop = max(p, p ** power)

        if q_min_lop <= q_max:
            lopsided += bisect_right(primes, q_max) - bisect_left(primes, q_min_lop)

    return 100.0 * lopsided / total

def compute_curve(theta):
    maxN = max(N)
    primes = sieve(maxN // 2)
    return [count_ratio(n, theta, primes) for n in N]

def main():
    plt.figure(figsize=(9, 5))

    # θ set
    thetas = [0.20, 0.25, 0.30, 0.35]

    for t in thetas:
        if t == 0.25:
            y = P_025
        else:
            y = compute_curve(t)

        plt.plot(N, y, marker="o", label=f"θ = {t:.2f}")

    plt.xscale("log")
    plt.xlabel("N cutoff (log scale)")
    plt.ylabel("Lopsided semiprimes (%)")
    plt.title("Lopsided Semiprimes vs N for multiple θ")
    plt.grid(True, alpha=0.3)
    plt.legend()

    plt.tight_layout()
    plt.savefig("multi_theta_convergence.png", dpi=200)
    print("Saved: multi_theta_convergence.png")

if __name__ == "__main__":
    main()
