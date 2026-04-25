#!/usr/bin/env python3

"""
activation_steps.py

Detect and visualize discrete activation steps in P(N, θ0)
for semiprimes n = p*q (p <= q).

Outputs:
- CSV of P(N, θ0)
- Detected jump points
- Plot with predicted activation thresholds N = p^(1/θ0)
"""

import math
import csv
import argparse
import matplotlib.pyplot as plt


# ----------------------------
# Prime + SPF sieve
# ----------------------------
def sieve_spf(limit):
    spf = list(range(limit + 1))
    for i in range(2, int(limit ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, limit + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def is_semiprime(n, spf):
    p = spf[n]
    q = n // p
    return p * q == n and spf[q] == q


# ----------------------------
# Compute P(N, θ0)
# ----------------------------
def compute_P_values(max_n, step, theta0):
    spf = sieve_spf(max_n)

    Ns = []
    Pvals = []

    total = 0
    good = 0

    for n in range(2, max_n + 1):
        if is_semiprime(n, spf):
            total += 1
            p = spf[n]

            theta = math.log(p) / math.log(n)
            if theta <= theta0:
                good += 1

        if n % step == 0:
            if total > 0:
                Ns.append(n)
                Pvals.append(good / total)

    return Ns, Pvals


# ----------------------------
# Detect jumps (step increases)
# ----------------------------
def detect_jumps(Ns, Pvals, threshold=0.002):
    jumps = []
    for i in range(1, len(Pvals)):
        delta = Pvals[i] - Pvals[i - 1]
        if delta > threshold:
            jumps.append((Ns[i], delta))
    return jumps


# ----------------------------
# Predict activation thresholds
# ----------------------------
def predicted_thresholds(theta0, max_n):
    primes = simple_primes(1000)
    thresholds = []

    for p in primes:
        N = p ** (1 / theta0)
        if N <= max_n:
            thresholds.append((p, int(N)))

    return thresholds


def simple_primes(limit):
    sieve = [True] * (limit + 1)
    sieve[0:2] = [False, False]
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            for j in range(i * i, limit + 1, i):
                sieve[j] = False
    return [i for i in range(2, limit + 1) if sieve[i]]


# ----------------------------
# Main
# ----------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=1_000_000)
    parser.add_argument("--step", type=int, default=1000)
    parser.add_argument("--theta", type=float, default=0.25)
    parser.add_argument("--out", default="activation_steps")

    args = parser.parse_args()

    print(f"Computing P(N, θ0={args.theta}) up to {args.max_n}...")

    Ns, Pvals = compute_P_values(args.max_n, args.step, args.theta)

    jumps = detect_jumps(Ns, Pvals)
    thresholds = predicted_thresholds(args.theta, args.max_n)

    # ----------------------------
    # Save CSV
    # ----------------------------
    with open(f"{args.out}.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["N", "P(N, theta0)"])
        for n, p in zip(Ns, Pvals):
            writer.writerow([n, p])

    # ----------------------------
    # Plot
    # ----------------------------
    plt.figure(figsize=(10, 6))
    plt.plot(Ns, Pvals, label="P(N, θ0)", linewidth=2)

    # Plot jumps
    for n, delta in jumps:
        plt.axvline(n, linestyle="--", alpha=0.3)

    # Plot predicted thresholds
    for p, n in thresholds:
        plt.axvline(n, color="red", alpha=0.2)
        plt.text(n, max(Pvals)*0.95, str(p), rotation=90, fontsize=6)

    plt.title(f"Activation Steps in P(N, θ0={args.theta})")
    plt.xlabel("N")
    plt.ylabel("P(N, θ0)")
    plt.legend()

    plt.tight_layout()
    plt.savefig(f"{args.out}.png", dpi=150)
    plt.show()

    # ----------------------------
    # Console summary
    # ----------------------------
    print("\nDetected jumps:")
    for n, d in jumps[:10]:
        print(f"N={n}, Δ={d:.4f}")

    print("\nPredicted thresholds (first 10):")
    for p, n in thresholds[:10]:
        print(f"p={p}, N≈{n}")


if __name__ == "__main__":
    main()
