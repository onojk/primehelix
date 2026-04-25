#!/usr/bin/env python3
"""
Visualize discrete activation thresholds for semiprime lopsidedness.
Based on paper Section 3.5: primes activate at N = p^(1/θ₀),
creating step-like jumps in P(N,θ₀).

Author: Based on Kendall (2026) - On the Distribution of Smaller Factors in Semiprimes
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple
import sympy

def activation_threshold(prime: int, theta: float) -> float:
    """N = p^(1/θ) where prime p begins contributing"""
    return prime ** (1.0 / theta)

def find_activation_primes(theta: float, N_max: float) -> List[Tuple[int, float]]:
    """Find all primes p where p^(1/θ) ≤ N_max"""
    max_prime = int(N_max ** theta) + 100
    primes = list(sympy.primerange(2, max_prime))
    
    activations = []
    for p in primes:
        N_activate = activation_threshold(p, theta)
        if N_activate <= N_max:
            activations.append((p, N_activate))
    
    return sorted(activations, key=lambda x: x[1])

def compute_P_N_theta_with_steps(N_values: List[int], theta: float, 
                                  primes_activations: List[Tuple[int, float]]) -> List[float]:
    """
    Compute P(N,θ) with explicit step detection.
    This is a heuristic model showing the discrete nature.
    """
    results = []
    total_semiprimes_est = 0
    
    for N in N_values:
        # Estimate cumulative contributions from primes that have activated
        cumulative = 0
        for p, N_activate in primes_activations:
            if N >= N_activate:
                # Each activated prime contributes approximately π(N/p) - π(p-1)
                # Simplified: roughly (N/p) / log(N/p) for large N
                if N > p:
                    approx_contrib = (N / p) / max(1, math.log(N / p))
                    cumulative += approx_contrib
        
        # Estimate total semiprimes up to N: ~ N log log N / log N
        if N > 2:
            total_semiprimes_est = N * math.log(max(2, math.log(N))) / max(1, math.log(N))
        
        proportion = min(1.0, cumulative / max(1, total_semiprimes_est))
        results.append(proportion)
    
    return results

def plot_activation_steps(theta: float, N_max: int = 10**8):
    """Generate Figure showing discrete activation thresholds"""
    
    # Get activation primes
    activations = find_activation_primes(theta, N_max)
    
    print(f"\nθ = {theta}")
    print(f"Activation thresholds for first few primes:")
    for p, N_act in activations[:10]:
        print(f"  p = {p:3d}  →  N_activate = {N_act:.2e}")
    
    # Create N values on log scale
    N_values = np.logspace(4, 8, 500).astype(int)
    
    # Compute P(N,θ)
    P_values = compute_P_N_theta_with_steps(N_values, theta, activations)
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: P(N,θ) with visible steps
    ax1.semilogx(N_values, P_values, 'b-', linewidth=1.5, alpha=0.7)
    
    # Mark activation thresholds
    for p, N_act in activations[:15]:  # First 15 primes
        y_at_activation = np.interp(np.log10(N_act), np.log10(N_values), P_values)
        ax1.axvline(x=N_act, color='red', linestyle=':', alpha=0.3, linewidth=0.8)
        if p <= 31:  # Label only small primes
            ax1.text(N_act, y_at_activation + 0.02, f'p={p}', 
                    rotation=90, fontsize=7, alpha=0.7)
    
    ax1.set_xlabel('N (log scale)', fontsize=12)
    ax1.set_ylabel(f'P(N, θ={theta})', fontsize=12)
    ax1.set_title(f'Discrete Activation Steps at N = p$^{{1/θ}}$\nθ = {theta}', fontsize=11)
    ax1.grid(True, alpha=0.3, which='both')
    ax1.set_ylim(0, 1)
    
    # Plot 2: Derivative approximation (shows step locations more clearly)
    dP = np.diff(P_values)
    dN = np.diff(np.log10(N_values))
    derivative = dP / dN
    
    ax2.plot(np.log10(N_values[1:]), derivative, 'g-', alpha=0.7, linewidth=1)
    ax2.set_xlabel('log₁₀(N)', fontsize=12)
    ax2.set_ylabel('dP/d(log N) (approx)', fontsize=12)
    ax2.set_title('Step Detection via Derivative', fontsize=11)
    ax2.grid(True, alpha=0.3)
    
    # Mark steps on derivative plot
    for p, N_act in activations[:15]:
        if N_act >= 10**4 and N_act <= N_max:
            ax2.axvline(x=math.log10(N_act), color='red', linestyle=':', alpha=0.3)
            if p <= 13:
                ax2.text(math.log10(N_act), max(derivative)*0.8, f'p={p}', 
                        rotation=90, fontsize=7, alpha=0.7)
    
    plt.tight_layout()
    
    # Save figure
    filename = f'activation_steps_theta_{theta:.2f}.png'
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {filename}")
    
    return fig

def multi_theta_comparison(thetas: List[float], N_max: int = 10**8):
    """Compare activation patterns across different θ values"""
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(thetas)))
    
    for idx, theta in enumerate(thetas):
        ax = axes[idx]
        
        activations = find_activation_primes(theta, N_max)
        N_values = np.logspace(4, 8, 300).astype(int)
        P_values = compute_P_N_theta_with_steps(N_values, theta, activations)
        
        ax.semilogx(N_values, P_values, color=colors[idx], linewidth=1.5, 
                   label=f'θ={theta}')
        
        # Mark only the most significant steps
        for p, N_act in activations[:8]:
            if p in [2, 3, 5, 7]:  # Only smallest primes for clarity
                ax.axvline(x=N_act, color='red', linestyle=':', alpha=0.2, linewidth=0.5)
        
        ax.set_xlabel('N (log scale)', fontsize=10)
        ax.set_ylabel(f'P(N,θ)', fontsize=10)
        ax.set_title(f'θ = {theta} (first steps at p=2,3,5,7)', fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9)
        ax.set_ylim(0, 1)
    
    plt.suptitle('Activation Thresholds Across Different θ Values\n' +
                f'N_max = {N_max:.0e}', fontsize=12)
    plt.tight_layout()
    
    filename = 'activation_steps_comparison.png'
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {filename}")

def table_of_activation_thresholds(theta: float = 0.10, max_N: int = 10**9):
    """Generate table showing prime activation thresholds (like paper's Section 3.5)"""
    
    print("\n" + "="*60)
    print(f"Activation Thresholds for θ = {theta}")
    exponent = 1.0 / theta
    print(f"N = p^(1/θ) = p^{exponent:.1f}")
    print("="*60)
    print(f"{'Prime p':<10} {'N_activate':<15} {'log10(N)':<12} {'In Range?'}")
    print("-"*60)
    
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    
    for p in primes:
        N_act = p ** (1.0/theta)
        in_range = N_act <= max_N
        status = "✓" if in_range else "✗ (beyond range)"
        
        # Match paper's exact numbers for θ=0.10
        if theta == 0.10 and p == 2:
            print(f"{p:<10} {N_act:<15.0f} {math.log10(N_act):<12.2f} {status} (paper: 1024)")
        elif theta == 0.10 and p == 3:
            print(f"{p:<10} {N_act:<15.0f} {math.log10(N_act):<12.2f} {status} (paper: ~59000)")
        elif theta == 0.10 and p == 5:
            print(f"{p:<10} {N_act:<15.2e} {math.log10(N_act):<12.2f} {status} (paper: ~9.8×10⁶)")
        elif theta == 0.10 and p == 7:
            print(f"{p:<10} {N_act:<15.2e} {math.log10(N_act):<12.2f} {status} (paper: ~2.8×10⁸)")
        else:
            print(f"{p:<10} {N_act:<15.2e} {math.log10(N_act):<12.2f} {status}")

if __name__ == "__main__":
    print("="*60)
    print("Semiprime Activation Threshold Visualizer")
    print("Based on Kendall (2026) Section 3.5")
    print("Discrete steps at N = p^(1/θ)")
    print("="*60)
    
    # Generate paper's example for θ=0.10 (mentioned in paper)
    print("\n[1] Reproducing paper example: θ = 0.10")
    table_of_activation_thresholds(theta=0.10, max_N=10**9)
    
    # Generate visualization for θ=0.25 (from Table 1)
    print("\n[2] Generating activation step plots...")
    fig1 = plot_activation_steps(theta=0.25, N_max=10**8)
    
    # Generate for θ=0.20 (from paper's Table 1)
    fig2 = plot_activation_steps(theta=0.20, N_max=10**8)
    
    # Multi-theta comparison
    print("\n[3] Generating multi-theta comparison...")
    thetas_to_compare = [0.10, 0.20, 0.30, 0.40]
    multi_theta_comparison(thetas_to_compare, N_max=10**8)
    
    print("\n✅ All visualizations saved as PNG files.")
    print("\nKey observation from paper (Section 3.5):")
    print("  The convergence of P(N,θ₀) is not smooth but proceeds")
    print("  in steps indexed by the prime activation sequence {p^(1/θ₀)}.")
    print("  This explains why smooth models perform less well at small θ₀.")
    
    plt.show()
