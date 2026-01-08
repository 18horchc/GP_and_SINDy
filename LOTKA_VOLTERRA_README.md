# Lotka-Volterra GPSINDy Comparison

This script compares two approaches for sparse noisy data system identification, implementing the GPSINDy methodology from Hsin et al. (2025):

## Data Generation

- Generates true Lotka-Volterra solution with parameters: α=1.5, β=1.0, γ=3.0, δ=1.0
- Samples 6 sparse time points evenly spaced from t=0 to t=10
- Creates 3 independent noisy measurements at each time point (18 total data points)
- Adds 10% Gaussian noise to measurements

## Path 1: Direct ESINDy (Discrete-Time)

- **Input**: 6 averaged data points (averages of 3 replicates at each of 6 unique time points)
- **Formulation**: Discrete-time: X(t+1) = f(X(t))
- **Library**: [1, x, y, x², y², xy] (6 functions)
- **Optimization**: ADMM for LASSO (matching Hsin et al. 2025)
- **Ensemble**: 1000 ensemble members, each using 5 out of 6 library functions
- **Note**: No derivatives computed - uses next-state prediction directly

## Path 2: GP-Enhanced ESINDy (Continuous-Time)

- **Input**: 18 sparse noisy points (3 replicates at 6 time points)
- **GP Fitting**: 
  - Fits squared exponential GP to all 18 points
  - Generates 50 dense time points for augmentation
  - Samples from GP posterior distribution
- **Data Augmentation**:
  - Combines GP mean and GP samples (interleaved)
  - Replaces values at original time points with mean of replicates
  - Results in 50 combined data points
- **Derivatives**: 
  - Uses **GP analytical derivatives** for ALL points (matching Hsin et al. 2025)
  - Computed via `gp_derivative_analytical()` function
  - No finite differences used
- **Library**: Expanded polynomial library [1, x, y, x², y², xy, x³, y³, x²y, xy²] (10 functions)
- **Optimization**: ADMM for LASSO (matching Hsin et al. 2025)
- **Ensemble**: 1000 ensemble members, each using 9 out of 10 library functions
- **Formulation**: Continuous-time: dX/dt = f(X)

## Key Differences from Original README

1. **Optimization Method**: Now uses **ADMM** instead of STLS (Sequential Thresholded Least Squares)
2. **GP Derivatives**: Uses **analytical GP derivatives** for all points, not finite differences
3. **Data Structure**: Path 1 uses 6 averaged points; Path 2 uses 18 raw points for GP fitting, then 50 augmented points
4. **Methodology**: Aligned with Hsin et al. (2025) GPSINDy approach

## Usage

Simply run:
```matlab
lotka_volterra_comparison
```

## Output

The script generates:
1. **Figure 1**: True Lotka-Volterra solution vs sparse noisy data (18 points)
2. **Figure 2**: GP fit and sampling visualization
3. **Figure 3**: Model predictions comparison (Path 1 vs Path 2 vs True)
4. **Console output**: 
   - Model coefficients for both paths
   - RMSE metrics comparing predictions to ground truth
   - True Lotka-Volterra coefficients for reference

## Parameters to Adjust

- `alpha, beta, gamma, delta`: Lotka-Volterra parameters (lines ~15-18)
- `noise_level`: Gaussian noise level (line ~44, default: 0.1 = 10%)
- `num_replicates`: Number of measurements per time point (line ~43, default: 3)
- `N`: Number of ESINDy ensemble members (line ~116, default: 1000)
- `t_dense`: Number of GP-sampled points (line ~237, default: 50)
- `lambda`: LASSO sparsity parameter (lines ~139, ~437, default: 0.01)
- `admm_options`: ADMM parameters (rho, tolerances, max iterations)

## ADMM Configuration

Default ADMM settings (matching Hsin et al. 2025):
- `rho = 1.0`: ADMM penalty parameter
- `max_iterations = 1000`: Maximum ADMM iterations
- `abs_tol = 1e-4`: Absolute convergence tolerance
- `rel_tol = 1e-2`: Relative convergence tolerance
- `verbose = false`: Suppress convergence messages

## Notes

- **Path 1** uses discrete-time ESINDy (next state prediction) with ADMM optimization
- **Path 2** uses continuous-time approach with GP analytical derivatives and ADMM optimization
- **GP derivatives** are computed analytically from the GP kernel, providing denoised, uncertainty-quantified derivatives
- The expanded library in Path 2 allows for more complex models when more data is available
- Both paths use ensemble methodology (1000 members) to handle uncertainty
- Implementation follows Hsin et al. (2025) GPSINDy methodology:
  - GP analytical derivatives for all augmented data points
  - ADMM for LASSO optimization
  - Consistent approach across both paths

## True Lotka-Volterra Coefficients

For reference, the true system is:
- dx/dt = 1.5x - 1.0xy
- dy/dt = 1.0xy - 3.0y

In library form [1, x, y, x², y², xy]:
- dx/dt: [0, 1.5, 0, 0, 0, -1.0]
- dy/dt: [0, 0, -3.0, 0, 0, 1.0]
