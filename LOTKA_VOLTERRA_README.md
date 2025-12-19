# Lotka-Volterra ESINDy Comparison

This script compares two approaches for sparse noisy data system identification:

## Path 1: Direct ESINDy
- Takes 6 sparse noisy points from Lotka-Volterra equations
- Runs ESINDy ensemble directly on these 6 points
- Uses discrete-time formulation: X(t+1) = f(X(t))
- Library: [1, x, y, x², y², xy] (6 functions)

## Path 2: GP-enhanced ESINDy
- Fits Gaussian Process to the 6 sparse noisy points
- Samples additional data points from the GP (50 total points)
- Computes GP derivative using finite differences
- Combines original 6 points with GP-sampled data
- Runs ESINDy ensemble with GP derivatives (continuous-time)
- Expanded library: [1, x, y, x², y², xy, x³, y³, x²y, xy²] (10 functions)

## Usage

Simply run:
```matlab
lotka_volterra_comparison
```

## Output

The script generates:
1. **Figure 1**: True Lotka-Volterra solution vs sparse noisy data
2. **Figure 2**: GP fit and sampling visualization
3. **Figure 3**: Model predictions comparison (Path 1 vs Path 2 vs True)
4. **Console output**: Model coefficients, RMSE metrics, and comparison

## Parameters to Adjust

- `alpha, beta, gamma, delta`: Lotka-Volterra parameters (line ~15)
- `noise_level`: Gaussian noise level (line ~30)
- `N`: Number of ESINDy ensemble members (line ~92)
- `t_dense`: Number of GP-sampled points (line ~162)
- `lambda`: Sparse regression threshold (line ~127)

## Notes

- Path 1 uses discrete-time ESINDy (next state prediction)
- Path 2 uses continuous-time approach with GP derivatives
- The expanded library in Path 2 allows for more complex models
- Both paths use ensemble methodology to handle uncertainty

