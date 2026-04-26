import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# =========================
# 1. Activation threshold plot (theta0 = 0.30)
# =========================
def plot_activation(df, theta0=0.30):
    col_name = f'P_{theta0:.2f}'
    if col_name not in df.columns:
        raise ValueError(f"Column {col_name} not found. Available: {df.columns.tolist()}")
    
    N_vals = df['N']
    P_vals = df[col_name]
    
    plt.figure(figsize=(8,5))
    plt.semilogx(N_vals, P_vals, 'o-', color='blue', linewidth=2, markersize=8,
                 label=f'$\\theta_0 = {theta0:.2f}$')
    
    ymin, ymax = plt.ylim()
    label_y = ymax - 0.02 * (ymax - ymin)
    xmin, xmax = plt.xlim()
    
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    for p in primes:
        x_thresh = p ** (1.0 / theta0)
        if xmin <= x_thresh <= xmax:
            plt.axvline(x=x_thresh, linestyle='--', color='gray', alpha=0.7, linewidth=1)
            plt.text(x_thresh, label_y, f'p={p}',
                     rotation=90, ha='center', va='bottom', fontsize=8, color='darkred')
    
    plt.xlabel('$N$', fontsize=12)
    plt.ylabel(f'$P(N, {theta0:.2f})$ (%)', fontsize=12)
    plt.title(f'Activation thresholds for $\\theta_0 = {theta0:.2f}$', fontsize=14)
    plt.grid(True, which='both', linestyle=':', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    # Save with the name expected by paper.tex
    plt.savefig('activation_theta_0.30.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved activation_theta_0.30.png")

# =========================
# 2. Theta vs log(p/q) plot
# =========================
def plot_theta_vs_ratio():
    theta = np.linspace(0.01, 0.5, 500)
    log_ratio = (2*theta - 1) * np.log(1e8)   # using N=1e8 as reference
    plt.figure(figsize=(8,5))
    plt.plot(theta, log_ratio, 'k-', linewidth=2)
    plt.xlabel('$\\theta = \\log p / \\log n$', fontsize=12)
    plt.ylabel('$\\log(p/q)$', fontsize=12)
    plt.title('Factor ratio vs $\\theta$', fontsize=14)
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.tight_layout()
    plt.savefig('theta_vs_logratio.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved theta_vs_logratio.png")

# =========================
# 3. Stabilisation test plot
# =========================
def plot_stabilisation(df, theta0=0.30):
    # Example: plot delta * log(log(N)) versus N
    # Replace with your actual stabilisation metric if needed.
    col_name = f'P_{theta0:.2f}'
    N_vals = df['N']
    P_vals = df[col_name]
    # Assuming a target limit L (here using last value as proxy)
    L_est = P_vals.iloc[-1]
    delta = L_est - P_vals
    loglogN = np.log(np.log(N_vals))
    plt.figure(figsize=(8,5))
    plt.semilogx(N_vals, delta * loglogN, 'o-', color='green')
    plt.xlabel('$N$', fontsize=12)
    plt.ylabel('$\\delta \\cdot \\log\\log N$', fontsize=12)
    plt.title('Stabilisation test', fontsize=14)
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.tight_layout()
    plt.savefig('stabilisation_test.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved stabilisation_test.png")

# =========================
# Main execution
# =========================
if __name__ == "__main__":
    # Load data (ensure semiprime_data.csv is in the same directory)
    df = pd.read_csv('semiprime_data.csv')
    
    plot_activation(df, theta0=0.30)
    plot_theta_vs_ratio()
    plot_stabilisation(df, theta0=0.30)
