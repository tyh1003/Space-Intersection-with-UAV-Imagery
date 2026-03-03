# Space Intersection with UAV Imagery (Photogrammetry, MATLAB)

## Overview
This project implements a complete Space Intersection (SI) pipeline for estimating 3D ground point coordinates from multi-view UAV imagery. It includes forward image measurement simulation using collinearity equations, Gaussian noise injection, iterative least-squares adjustment, covariance-based uncertainty estimation, and gross-error testing.

## What This Project Does
Given:
- UAV camera positions and orientations per image
- Image measurements (or simulated image points)

It estimates:
- 3D ground point coordinates (X, Y, Z)
- Posterior unit weight standard deviation (σ₀)
- Parameter covariance matrix Σxx
- 2D/3D uncertainty visualization (error ellipse / error ellipsoid)
- Robustness under gross errors (3σ, 5σ, 10σ, 100σ scenarios)

## Methodology

### 1) Image Point Generation (Collinearity Equations)
- Compute ideal image coordinates (xp, yp) from collinearity equations
- Used as reference “truth” for simulation

### 2) Measurement Noise Simulation
- Generate N(0, σ) Gaussian noise (σ = 20 μm)
- Validate distribution using summary statistics and histogram

### 3) Space Intersection via Iterative Least Squares
- Linearize observation equations and solve normal equations:
  V = AX − L  
  N = AᵀPA,  U = AᵀPL,  NX = U
- Solvers used:
  - Cholesky factorization
  - Gauss-Jordan elimination
- Iterate until convergence:
  max |Δ| ≤ 0.002 m

### 4) Uncertainty Quantification
- Compute σ₀ (posterior unit weight std.)
- Derive Σxx and perform eigen-decomposition
- Visualize:
  - 3D error ellipsoid (λ₁, λ₂, λ₃)
  - 2D error ellipse (λ₁, λ₂)

### 5) Gross Error Injection and Detection
- Inject outliers into selected image measurements
- Compare convergence behavior and residual patterns
- Use standardized residuals to identify potential gross errors

## Repository Structure
- `HW3.m` — main MATLAB implementation
- `pho.xlsx`, `pho2.xlsx`, `UAV_Coordinates_and_Angles.xlsx` — UAV trajectory, orientation, and measurement inputs
- `*.tif` — error ellipsoid/ellipse plots and histogram outputs
- `HW3.pdf` — report with derivations and analysis

## Notes on Reproducibility
- Input spreadsheets define the test configuration (UAV pose + observations)
- Output figures summarize distribution checks and uncertainty geometry

## References
- 尤瑞哲，《基礎測量平差法》，2018，第5.6節
- Course slides: Space Intersection
