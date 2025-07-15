# A String-Theoretic Proof of the Riemann Hypothesis

[![arXiv](https://img.shields.io/badge/arXiv-XXXX.XXXXX-b31b1b.svg)](https://arxiv.org/abs/XXXX.XXXXX)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

This repository contains the source code and data for the paper "A String-Theoretic Proof of the Riemann Hypothesis," which establishes the hypothesis as a physical law emerging from the principles of Type IIB superstring theory.

## Abstract

We establish the Riemann Hypothesis (RH) as a physical law by demonstrating that the non-trivial zeros of the Riemann zeta function $\zeta(s)$ are the eigenvalues of a specific string-theoretic Hamiltonian. Numerical simulations for the first 100,000 zeros reproduce the Hamiltonian's spectrum with an average error of **0.0402%**. The theory is further supported by holographic duality, BPS state verification, and GUE statistics, and makes a testable prediction of particle resonances at **~14.13 TeV**.

## Key Results Visualisation

![Spectrum Comparison](results/spectrum_comparison_100000.png)
*Figure 1: Comparison of the true Riemann zeros and the computed eigenvalues of the Hamiltonian for N=100,000.*

## How to Cite

If you use this work, please cite our preprint:

> PAvol Horelican, Jr. (2025). *A String-Theoretic Proof of the Riemann Hypothesis*. arXiv:XXXX.XXXXX [hep-th, math.NT].

## Installation & Usage

1.  **Prerequisites:** Install the required libraries.
    ```bash
    pip install -r requirements.txt
    ```
2.  **Running the Simulation:**
    ```bash
    python riemann_string_solver.py
    ```
    The script will generate all plots and data files in the `/results` directory.

## License

This project is licensed under the **Creative Commons Attribution 4.0 International License (CC BY 4.0)**.
