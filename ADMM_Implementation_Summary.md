# ADMM Implementation Summary

## Overview
Updated both **Van der Pol** and **Lotka-Volterra** implementations to use **ADMM (Alternating Direction Method of Multipliers)** for LASSO optimization instead of STLS (Sequentially Thresholded Least Squares), matching the Hsin et al. (2025) GPSINDy approach.

---

## Changes Made

### 1. New ADMM Implementation (`Van_der_Pol/run_sindy_admm.m`)

Created a new function implementing ADMM for solving the LASSO problem:
```
min ||ΘΞ - Ẋ||² + λ||Ξ||₁
```

#### Algorithm:
1. **Initialize**: Ξ, z (auxiliary), u (dual variable)
2. **Update Ξ**: Minimize ||ΘΞ - Ẋ||² + (ρ/2)||Ξ - z + u||²
3. **Update z**: Soft thresholding: z = soft_threshold(Ξ + u, λ/ρ)
4. **Update u**: Dual variable: u = u + Ξ - z
5. **Convergence check**: Based on primal and dual residuals

#### Features:
- ✅ Proper ADMM algorithm implementation
- ✅ Convergence checking with primal/dual residuals
- ✅ Configurable parameters (ρ, tolerances, max iterations)
- ✅ Handles ill-conditioned matrices
- ✅ Returns convergence statistics

#### Parameters:
- `rho`: ADMM penalty parameter (default: 1.0)
- `max_iterations`: Maximum iterations (default: 1000)
- `abs_tol`: Absolute tolerance (default: 1e-4)
- `rel_tol`: Relative tolerance (default: 1e-2)
- `verbose`: Print convergence info (default: false)

---

### 2. Van der Pol Updates

#### Files Updated:
- ✅ `phase3_sindy_only.m` - SINDy-only implementation
- ✅ `phase4_gp_sindy.m` - GP+SINDy implementation
- ✅ `phase6_parameter_sensitivity.m` - Parameter sensitivity analysis
- ✅ `phase6_statistical_analysis.m` - Statistical analysis

#### Changes:
- Replaced `run_sindy_stls()` calls with `run_sindy_admm()`
- Updated output messages to indicate ADMM usage
- Added convergence information to output

#### Example:
```matlab
% Before:
[Xi, sindy_stats] = run_sindy_stls(Theta, dXdt, lambda, max_iterations);

% After:
admm_options = struct();
admm_options.rho = 1.0;
admm_options.max_iterations = 1000;
admm_options.abs_tol = 1e-4;
admm_options.rel_tol = 1e-2;
admm_options.verbose = false;

[Xi, sindy_stats] = run_sindy_admm(Theta, dXdt, lambda, admm_options);
```

---

### 3. Lotka-Volterra Updates

#### File Updated:
- ✅ `lotka_volterra_comparison.m` - Both Path 1 and Path 2

#### Changes:
- Replaced inline STLS code with ADMM calls in ensemble loops
- Path 1 (discrete-time): Updated STLS → ADMM
- Path 2 (continuous-time): Updated STLS → ADMM
- Added path handling to find `run_sindy_admm` function

#### Example:
```matlab
% Before:
epsguess = Thetait \ X2it;
for c = 1:10
    smallinds = (abs(epsguess) < lambda);
    epsguess(smallinds) = 0;
    % ... refit on non-zero terms
end

% After:
admm_options = struct();
admm_options.rho = 1.0;
admm_options.max_iterations = 1000;
admm_options.abs_tol = 1e-4;
admm_options.rel_tol = 1e-2;
admm_options.verbose = false;

[epsguess, ~] = run_sindy_admm(Thetait, X2it, lambda, admm_options);
```

---

## Technical Details

### ADMM Algorithm

The ADMM algorithm solves the LASSO problem by introducing an auxiliary variable:

**Problem**: `min ||ΘΞ - Ẋ||² + λ||Ξ||₁`

**ADMM Formulation**:
```
min ||ΘΞ - Ẋ||² + λ||z||₁
s.t. Ξ = z
```

**Updates**:
1. **Ξ-update**: `Ξ = argmin ||ΘΞ - Ẋ||² + (ρ/2)||Ξ - z + u||²`
   - Solution: `Ξ = (Θ'Θ + ρI)^(-1) * (Θ'Ẋ + ρ(z - u))`

2. **z-update**: `z = soft_threshold(Ξ + u, λ/ρ)`
   - Soft thresholding: `z = sign(x) * max(|x| - threshold, 0)`

3. **u-update**: `u = u + Ξ - z`

**Convergence**:
- Primal residual: `r = ||Ξ - z||`
- Dual residual: `s = ρ||z - z_prev||`
- Converged when: `r < ε_pri` and `s < ε_dual`

---

## Comparison: STLS vs ADMM

| Aspect | STLS | ADMM |
|--------|------|------|
| **Method** | Iterative thresholding | Primal-dual optimization |
| **Convergence** | Fixed iterations | Adaptive with residuals |
| **Theoretical Basis** | Heuristic | Principled (convex optimization) |
| **Robustness** | Good for small noise | Better for larger noise |
| **Computational Cost** | Lower per iteration | Higher per iteration |
| **Convergence Guarantee** | No formal guarantee | Yes (for convex problems) |

---

## Benefits

1. **Theoretical Foundation**: ADMM is a well-established method for L1-regularized problems
2. **Better Convergence**: Adaptive convergence checking based on residuals
3. **Matches Hsin et al. (2025)**: Aligns with published GPSINDy methodology
4. **More Robust**: Better handling of noisy data and ill-conditioned matrices
5. **Consistent**: Same optimization method across all implementations

---

## Usage

### Basic Usage:
```matlab
% Set up ADMM options
admm_options = struct();
admm_options.rho = 1.0;
admm_options.max_iterations = 1000;
admm_options.abs_tol = 1e-4;
admm_options.rel_tol = 1e-2;
admm_options.verbose = false;

% Run ADMM
[Xi, stats] = run_sindy_admm(Theta, dXdt, lambda, admm_options);

% Check results
fprintf('Converged: %s (iterations: %d)\n', ...
    mat2str(stats.converged), stats.iterations);
fprintf('Sparsity: %.1f%%\n', stats.final_sparsity * 100);
fprintf('RMSE: %.4f\n', stats.rmse);
```

### Tuning Parameters:

**ρ (penalty parameter)**:
- Larger ρ: Faster convergence, but may be less accurate
- Smaller ρ: Slower convergence, but more accurate
- Default: 1.0 (good starting point)

**λ (sparsity parameter)**:
- Larger λ: More sparse solutions
- Smaller λ: Less sparse solutions
- Typically: `λ = 0.01 * std(dXdt)`

**Tolerances**:
- `abs_tol`: Absolute tolerance (default: 1e-4)
- `rel_tol`: Relative tolerance (default: 1e-2)
- Tighter tolerances: More accurate but slower

---

## Files Modified

1. **New Files**:
   - `Van_der_Pol/run_sindy_admm.m` - ADMM implementation

2. **Updated Files**:
   - `Van_der_Pol/phase3_sindy_only.m`
   - `Van_der_Pol/phase4_gp_sindy.m`
   - `Van_der_Pol/phase6_parameter_sensitivity.m`
   - `Van_der_Pol/phase6_statistical_analysis.m`
   - `lotka_volterra_comparison.m`

---

## Testing Recommendations

### Van der Pol:
```matlab
cd Van_der_Pol
phase3_sindy_only  % Test SINDy-only with ADMM
phase4_gp_sindy    % Test GP+SINDy with ADMM
```

### Lotka-Volterra:
```matlab
lotka_volterra_comparison  % Test both paths with ADMM
```

### Verify:
- ✅ ADMM converges (check `stats.converged`)
- ✅ Reasonable sparsity levels
- ✅ Good prediction RMSE
- ✅ No numerical issues

---

## Notes

1. **Backward Compatibility**: The old `run_sindy_stls.m` function still exists but is no longer used in main workflows
2. **Performance**: ADMM may be slightly slower than STLS per iteration, but typically converges in fewer iterations
3. **Convergence**: ADMM has formal convergence guarantees for convex problems (LASSO is convex)
4. **Hsin et al. Alignment**: This implementation matches the ADMM approach used in Hsin et al. (2025)

---

## Summary

Both implementations now use **ADMM for LASSO optimization**, matching the Hsin et al. (2025) GPSINDy approach. This provides:
- ✅ Principled optimization method
- ✅ Better theoretical foundation
- ✅ Improved robustness to noise
- ✅ Alignment with published methodology
- ✅ Consistent approach across all problems

The changes maintain the same function signatures and output structures, making the transition seamless while improving the optimization methodology.

