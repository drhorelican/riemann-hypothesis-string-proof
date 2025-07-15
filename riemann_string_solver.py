# =============================================================================
# String-Theoretic Solver for the Riemann Hypothesis
# Author: Anonymous Collaborative Team
# Date: July 15, 2025
#
# Description:
# This script performs a numerical verification of the hypothesis that the
# non-trivial zeros of the Riemann zeta function are eigenvalues of a
# string-theoretic Hamiltonian.
#
# It performs the following steps:
# 1. Computes the first 100,000 non-trivial zeros of the Riemann zeta function.
# 2. Constructs a discretized string Hamiltonian with Calabi-Yau corrections.
# 3. Computes the first 100,000 eigenvalues of the Hamiltonian.
# 4. Compares the spectra and calculates the error.
# 5. Performs a statistical analysis (GUE) of the level spacings.
# 6. Generates plots and saves all data to CSV files in the /results/ directory.
# =============================================================================

import numpy as np
from mpmath import mp, zetazero
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
import pandas as pd
import os
from datetime import datetime

# --- Configuration ---
NUM_ZEROS = 100000          # Number of zeros/eigenvalues to compute
GRID_SIZE = NUM_ZEROS * 2   # Discretization grid size for the Hamiltonian (N)
MP_DPS = 25                 # Decimal precision for mpmath
RESULTS_DIR = "results"     # Directory to save outputs

# Hamiltonian Parameters (dimensionless, from calibration)
R_PARAM = 1.0213
LAMBDA_CY = 0.0427

def ensure_dir(directory):
    """Ensures that the specified directory exists."""
    if not os.path.exists(directory):
        os.makedirs(directory)

def compute_riemann_zeros(n_zeros):
    """Computes the imaginary parts of the first n non-trivial Riemann zeros."""
    print(f"[{datetime.now()}] Step 1: Computing {n_zeros} Riemann zeros...")
    mp.dps = MP_DPS
    zeros = [float(zetazero(i).imag) for i in range(1, n_zeros + 1)]
    print(f"[{datetime.now()}] Computation of zeros complete.")
    return np.array(zeros)

def build_string_hamiltonian(N, R, lambda_cy):
    """Constructs the sparse matrix for the string Hamiltonian."""
    print(f"[{datetime.now()}] Step 2: Building the string Hamiltonian (Grid Size N={N})...")
    # Kinetic term (finite difference approximation of Laplacian)
    main_diag = 2 * np.ones(N)
    off_diag = -1 * np.ones(N - 1)
    kinetic_term = diags([main_diag, off_diag, off_diag], [0, -1, 1], shape=(N, N))

    # Potential term from AdS curvature
    x = np.arange(N) - N // 2
    potential_term = diags((x / R)**2, 0, shape=(N, N))

    # Calabi-Yau correction term (simplified diagonal form)
    indices = np.arange(N)
    cy_correction_diag = lambda_cy * np.sin(np.pi * indices / N)**2
    cy_correction_term = diags(cy_correction_diag, 0, shape=(N, N))

    H = kinetic_term + potential_term + cy_correction_term
    print(f"[{datetime.now()}] Hamiltonian matrix built.")
    return H

def compute_eigenvalues(H, n_eigenvalues):
    """Computes the lowest eigenvalues of the Hamiltonian matrix."""
    print(f"[{datetime.now()}] Step 3: Diagonalizing the Hamiltonian to find eigenvalues...")
    # Using 'SM' for smallest magnitude eigenvalues, as they correspond to lowest energy states
    eigenvalues = eigsh(H, k=n_eigenvalues, which='SM', return_eigenvectors=False)
    print(f"[{datetime.now()}] Diagonalization complete.")
    return np.sort(eigenvalues)

def analyze_and_plot_results(zeros, eigenvalues, results_dir):
    """Analyzes the results and generates all plots and data files."""
    print(f"[{datetime.now()}] Step 4: Analyzing results and generating outputs...")

    # --- Error Analysis ---
    n_points = len(zeros)
    errors = np.abs(eigenvalues - zeros) / zeros * 100
    mean_error = np.mean(errors)
    print(f"--> Average Spectral Error: {mean_error:.4f}%")
    
    # --- DataFrames and CSV Export ---
    df_zeros = pd.DataFrame({'n': range(1, n_points + 1), 'gamma_n': zeros})
    df_eigen = pd.DataFrame({'n': range(1, n_points + 1), 'eigenvalue': eigenvalues})
    df_errors = pd.DataFrame({'n': range(1, n_points + 1), 'error_percent': errors})
    
    df_zeros.to_csv(f"{results_dir}/riemann_zeros_{n_points}.csv", index=False)
    df_eigen.to_csv(f"{results_dir}/eigenvalues_{n_points}.csv", index=False)
    df_errors.to_csv(f"{results_dir}/errors_{n_points}.csv", index=False)
    
    # --- Plot 1: Spectrum Comparison ---
    plt.figure(figsize=(12, 6))
    plt.plot(df_zeros['n'], df_zeros['gamma_n'], 'b-', label='Riemann Zeros $\\gamma_n$')
    plt.plot(df_eigen['n'], df_eigen['eigenvalue'], 'r--', label='Hamiltonian Eigenvalues $E_n$')
    plt.title(f'Comparison of Riemann Zeros and Hamiltonian Eigenvalues (N={n_points})')
    plt.xlabel('Index (n)')
    plt.ylabel('Value')
    plt.legend()
    plt.grid(True)
    plt.savefig(f"{results_dir}/spectrum_comparison_{n_points}.png")
    plt.close()

    # --- Plot 2: GUE Statistics ---
    def gue_distribution(s):
        return (32 / np.pi**2) * s**2 * np.exp(-4 * s**2 / np.pi)

    spacings = np.diff(zeros)
    normalized_spacings = spacings / np.mean(spacings)
    
    plt.figure(figsize=(10, 6))
    plt.hist(normalized_spacings, bins=50, density=True, label='Normalized Spacings of Zeros')
    s_vals = np.linspace(0, 3, 100)
    plt.plot(s_vals, gue_distribution(s_vals), 'r-', label='GUE Prediction')
    plt.title(f'GUE Statistics of Level Spacings (N={n_points})')
    plt.xlabel('Normalized Spacing (s)')
    plt.ylabel('Probability Density P(s)')
    plt.legend()
    plt.grid(True)
    plt.savefig(f"{results_dir}/gue_comparison_{n_points}.png")
    plt.close()

    # --- Plot 3: Error Plot ---
    plt.figure(figsize=(10, 6))
    plt.plot(df_errors['n'], df_errors['error_percent'])
    plt.xscale('log')
    plt.title(f'Percentage Error vs. Zero Index (N={n_points})')
    plt.xlabel('Index (n) [log scale]')
    plt.ylabel('Relative Error (%)')
    plt.grid(True)
    plt.savefig(f"{results_dir}/error_vs_n_{n_points}.png")
    plt.close()

    print(f"[{datetime.now()}] All plots and data files saved to '{results_dir}/'.")


if __name__ == "__main__":
    ensure_dir(RESULTS_DIR)
    
    # Execute the full workflow
    riemann_zeros = compute_riemann_zeros(NUM_ZEROS)
    hamiltonian_matrix = build_string_hamiltonian(GRID_SIZE, R_PARAM, LAMBDA_CY)
    hamiltonian_eigenvalues = compute_eigenvalues(hamiltonian_matrix, NUM_ZEROS)
    analyze_and_plot_results(riemann_zeros, hamiltonian_eigenvalues, RESULTS_DIR)
    
    print(f"[{datetime.now()}] Workflow complete.")
