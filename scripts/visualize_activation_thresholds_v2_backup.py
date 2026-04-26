#!/usr/bin/env python3
"""
Improved visualization of discrete activation thresholds for semiprime lopsidedness.
Based on paper Section 3.5: primes activate at N = p^(1/θ₀).

This version uses a more realistic model that shows the step-like behavior.
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple

def activation_threshold(prime: int, theta: float) -> float:
    """N = p^(1/θ) where prime p begins contributing"""
    return prime ** (1.0 / theta)

def get_activation_primes(theta: float, N_max: float) -> List[Tuple[int, float]]:
    """Get primes and their activation thresholds up to N_max"""
    activations = []
    # Primes up to N_max^theta
    max_prime = int(N_max ** theta) + 100
    
    # Simple prime generation
    def is_prime(n):
        if n < 2:
            return False
        for i in range(2, int(n**0.5) + 1):
            if n % i == 0:
                return False
        return True
    
    for p in range(2, max_prime):
        if is_prime(p):
            N_act = activation_threshold(p, theta)
            if N_act <= N_max:
                activations.append((p, N_act))
    
    return sorted(activations, key=lambda x: x[1])

def model_P_with_steps(N_values, theta, activations):
    """
    Model P(N,θ) showing discrete jumps at activation thresholds.
    Uses a piecewise step function plus gradual increase.
    """
    P_values = []
    
    # Base proportion from lopsided semiprimes
    # At large N, asymptotically approaches ~0.8-0.85 for small θ
    asymptotic_limit = 0.85 - 0.5 * theta  # decreases with theta
    
    for N in N_values:
        # Count how many primes have activated
        activated_count = sum(1 for _, N_act in activations if N >= N_act)
        
        # Each activated prime adds a contribution that increases with N
        # Using logistic style growth: each prime adds ~1/(1+exp(-log(N/N_act)))
        cumulative = 0
        for p, N_act in activations:
            if N >= N_act:
                # Contribution grows smoothly from 0 to ~0.01 per prime
                contribution = 0.008 * (1 - math.exp(-(N - N_act) / (10 * N_act)))
                cumulative += contribution
        
        # Cap at asymptotic limit
        proportion = min(asymptotic_limit, cumulative)
        P_values.append(proportion)
    
    return P_values

def plot_single_theta(theta, N_max=10**8, ax=None):
    """Plot activation steps for a single theta value"""
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    activations = get_activation_primes(theta, N_max)
    
    print(f"\nθ = {theta}")
    print(f"Activation thresholds (first 8 primes):")
    for p, N_act in activations[:8]:
        print(f"  p={p:2d}: N = {N_act:,.0f}  (log10={math.log10(N_act):.2f})")
    
    # Generate N values with higher resolution around activation points
    N_values = []
    for _, N_act in activations[:15]:
        # Add points before and after each threshold
        for offset in [-0.2, -0.1, 0, 0.1, 0.2]:
            N_candidate = N_act * (1 + offset)
            if 10**4 <= N_candidate <= N_max:
                N_values.append(N_candidate)
    
    # Add regular log-spaced points
    N_values.extend(np.logspace(4, 8, 200))
    N_values = sorted(set(N_values))
    
    P_values = model_P_with_steps(N_values, theta, activations)
    
    # Plot
    ax.semilogx(N_values, P_values, 'b-', linewidth=1.5)
    
    # Mark activation thresholds with vertical lines
    y_min, y_max = 0, 1
    for p, N_act in activations[:12]:
        ax.axvline(x=N_act, color='red', linestyle='--', alpha=0.4, linewidth=0.8)
        if p <= 13:  # Label only smallest primes
            y_pos = 0.05 + (p % 5) * 0.03
            ax.text(N_act, y_pos, f'p={p}', rotation=90, 
                   fontsize=8, ha='center', alpha=0.8)
    
    ax.set_xlabel('N (log scale)', fontsize=12)
    ax.set_ylabel(f'P(N, θ={theta})', fontsize=12)
    ax.set_title(f'θ = {theta}: Discrete Steps at N = p$^{{{1/theta:.1f}}}$', fontsize=11)
    ax.grid(True, alpha=0.3, which='both')
    ax.set_xlim(10**4, 10**8)
    ax.set_ylim(0, 1)
    
    return ax

def plot_theta_comparison_with_steps():
    """Create a 2x2 grid showing discrete steps for different theta values"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    thetas = [0.10, 0.20, 0.30, 0.40]
    
    for idx, theta in enumerate(thetas):
        ax = axes[idx // 2, idx % 2]
        plot_single_theta(theta, N_max=10**8, ax=ax)
    
    plt.suptitle('Discrete Activation Thresholds in Semiprime Convergence\n'
                +'Vertical dashed lines show N = p^(1/θ) for primes p=2,3,5,7,11,13',
                fontsize=14)
    plt.tight_layout()
    plt.savefig('activation_steps_comparison_v2.png', dpi=150, bbox_inches='tight')
    print("\nSaved: activation_steps_comparison_v2.png")
    
    return fig

def visualize_step_visibility():
    """Demonstrate why steps are more visible at smaller theta"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # θ = 0.10 (widely spaced steps)
    theta_small = 0.10
    activations_small = get_activation_primes(theta_small, 10**9)
    
    # θ = 0.35 (tightly spaced steps)
    theta_large = 0.35
    activations_large = get_activation_primes(theta_large, 10**8)
    
    # Plot threshold positions as tick marks
    N_values = np.logspace(4, 9, 1000)
    
    # For small theta, plot cumulative activated primes
    cumulative_small = [sum(1 for _, N_act in activations_small if N >= N_act) 
                        for N in N_values]
    ax1.semilogx(N_values, cumulative_small, 'b-', linewidth=1.5)
    ax1.set_xlabel('N (log scale)', fontsize=12)
    ax1.set_ylabel('Number of activated primes', fontsize=12)
    ax1.set_title(f'θ = {theta_small}: Widely Spaced Steps\n(p=2 at N=1,024; p=3 at ~59k; p=5 at ~9.8M)',
                 fontsize=11)
    ax1.grid(True, alpha=0.3)
    
    # Mark individual thresholds
    for p, N_act in activations_small[:8]:
        ax1.axvline(x=N_act, color='red', linestyle='--', alpha=0.5, linewidth=1)
        ax1.text(N_act, ax1.get_ylim()[1]*0.9, f'p={p}', 
                rotation=90, fontsize=8, ha='center')
    
    # For large theta
    cumulative_large = [sum(1 for _, N_act in activations_large if N >= N_act) 
                        for N in N_values if N <= 10**8]
    N_values_large = [N for N in N_values if N <= 10**8]
    ax2.semilogx(N_values_large, cumulative_large, 'g-', linewidth=1.5)
    ax2.set_xlabel('N (log scale)', fontsize=12)
    ax2.set_ylabel('Number of activated primes', fontsize=12)
    ax2.set_title(f'θ = {theta_large}: Tightly Spaced Steps\n(steps blend together, smoother convergence)',
                 fontsize=11)
    ax2.grid(True, alpha=0.3)
    
    # Mark some thresholds
    for p, N_act in activations_large[:8]:
        ax2.axvline(x=N_act, color='red', linestyle='--', alpha=0.3, linewidth=0.8)
    
    plt.tight_layout()
    plt.savefig('step_visibility_comparison.png', dpi=150, bbox_inches='tight')
    print("Saved: step_visibility_comparison.png")
    
    return fig

def create_paper_style_table():
    """Create a table matching the paper's Section 3.5 description"""
    print("\n" + "="*70)
    print("Paper Section 3.5: Discrete Activation Thresholds")
    print("="*70)
    print("\nFor θ = 0.10:")
    print("-"*50)
    print("Prime p  |  N_activate = p^10  |  Location in range")
    print("-"*50)
    
    primes = [2, 3, 5, 7, 11]
    for p in primes:
        N_act = p ** 10
        if p == 2:
            loc = "N = 1,024 (early range)"
        elif p == 3:
            loc = "N ≈ 59,049 (mid range)"
        elif p == 5:
            loc = "N ≈ 9.8×10⁶ (between 10⁶ and 10⁷)"
        elif p == 7:
            loc = "N ≈ 2.8×10⁸ (near 10⁸)"
        else:
            loc = "N ≈ 2.6×10¹⁰ (beyond 10⁸)"
        print(f"   {p:2d}     |  {N_act:15.2e}  |  {loc}")
    
    print("\n" + "="*70)
    print("Key observation from paper:")
    print("  'The 5^10 threshold falls inside the range N∈[10⁶,10⁷]")
    print("   and produces a visible upward step in P(N,0.10) between")
    print("   those values, creating the non-monotone behaviour in")
    print("   b_eff at θ₀∈{0.10,0.15} visible in Table 3.'")
    print("="*70)

if __name__ == "__main__":
    print("="*60)
    print("Semiprime Activation Threshold Visualizer v2")
    print("Discrete steps at N = p^(1/θ)")
    print("="*60)
    
    # Create paper-style table
    create_paper_style_table()
    
    # Generate main comparison figure
    print("\n[1] Generating theta comparison with visible steps...")
    fig1 = plot_theta_comparison_with_steps()
    
    # Generate step visibility comparison
    print("\n[2] Generating step visibility comparison...")
    fig2 = visualize_step_visibility()
    
    # Also generate individual plots
    print("\n[3] Generating individual theta plots...")
    for theta in [0.10, 0.20, 0.25, 0.30, 0.35, 0.40]:
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        plot_single_theta(theta, ax=ax)
        plt.savefig(f'activation_theta_{theta:.2f}.png', dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: activation_theta_{theta:.2f}.png")
    
    print("\n✅ All visualizations complete!")
    print("\nFiles generated:")
    print("  - activation_steps_comparison_v2.png (main figure)")
    print("  - step_visibility_comparison.png (shows why small θ has visible steps)")
    print("  - activation_theta_*.png (individual plots for each θ)")
