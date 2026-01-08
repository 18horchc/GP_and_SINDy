# GPSINDy Comparison: Hsin et al. (2025) vs Your Implementation

## Overview

This document compares the GPSINDy approach from **Hsin et al. (2025) - "Symbolic Regression on Sparse and Noisy Data with Gaussian Processes"** with your implementations for **Van der Pol** and **Lotka-Volterra** problems.

---

## Key Similarities

### 1. **Core Philosophy**
Both approaches share the fundamental idea:
- Use GP regression to smooth and denoise sparse, noisy data
- Compute derivatives analytically from the GP model
- Use these derivatives with SINDy to discover governing equations

### 2. **GP Kernel Choice**
- **Hsin et al.**: Uses squared exponential kernel (implied from their derivative formulation)
- **Your implementation**: Uses squared exponential kernel (`'squaredexponential'`) for both Van der Pol and Lotka-Volterra
- **Similarity**: ✅ **Identical**

### 3. **Analytical Derivative Computation**
Both use the analytical derivative of the GP mean:

**Hsin et al. approach** (from paper):
- Computes derivative of GP mean: `dμ/dt = K'(X_new, X) * K(X,X)^(-1) * y`
- Uses kernel derivative: `dk/dt = -(t_star - t_train) / L² * k(t_star, t_train)`

**Your implementation** (Van der Pol, Lotka-Volterra):
```matlab
% From fit_gp_vdp.m and run_GP_Interpolation.m
dk_dt = -(dist / (L^2)) .* k_star;
dXdt_gp(j, i) = dk_dt' * alpha;
```
- **Similarity**: ✅ **Identical mathematical approach**

### 4. **Data Augmentation Strategy**
- **Hsin et al.**: Samples from GP posterior to augment sparse data
- **Your Van der Pol**: Uses `augment_data_with_gp_samples()` which appends GP samples to original sparse data
- **Your Lotka-Volterra**: Uses GP mean and samples to create dense time grid
- **Similarity**: ✅ **Conceptually similar** (both augment sparse data with GP predictions)

### 5. **SINDy Library Construction**
- **Hsin et al.**: Uses polynomial library (order not explicitly stated, but likely 2-3)
- **Your Van der Pol**: Uses polynomial library up to order 3 (`polyOrder = 3`)
- **Your Lotka-Volterra**: Uses expanded library with cubic terms
- **Similarity**: ✅ **Similar polynomial libraries**

---

## Key Differences

### 1. **Sparsity Optimization Method**

**Hsin et al. (2025)**:
- Uses **ADMM (Alternating Direction Method of Multipliers)** to solve LASSO problem
- Solves: `min ||ΘΞ - Ẋ||² + λ||Ξ||₁`
- More sophisticated optimization with cross-validation for λ selection
- Handles each state variable's λ separately

**Your Implementation**:
- Uses **STLS (Sequentially Thresholded Least Squares)**
- Iteratively thresholds small coefficients and refits
- Simpler, more direct approach
- Single λ for all state variables (or data-adaptive: `λ = λ_base * std(dXdt)`)

**Difference**: ⚠️ **Different optimization algorithms**
- ADMM is more principled for L1 regularization
- STLS is simpler but may be less robust to noise

### 2. **Derivative Combination Strategy**

**Hsin et al. (2025)**:
- Uses **only GP analytical derivatives** for all augmented data points
- No mention of combining finite differences with GP derivatives

**Your Van der Pol Implementation**:
- **Hybrid approach**: Uses finite differences for sparse original points, GP analytical derivatives for dense GP-augmented points
- From `augment_data_with_gp_samples.m`:
  ```matlab
  dXdt_sparse_fd = compute_finite_differences(t_sparse, X_sparse, 'central');
  dXdt_aug = [dXdt_sparse_fd; dXdt_gp_combined];
  ```
- Then sorts by time and removes duplicates

**Your Lotka-Volterra Implementation**:
- Uses **only GP analytical derivatives** (similar to Hsin et al.)
- No finite differences for original sparse points

**Difference**: ⚠️ **Van der Pol uses hybrid, Lotka-Volterra matches Hsin et al.**

### 3. **Data Augmentation Details**

**Hsin et al. (2025)**:
- Samples from GP posterior distribution
- Uses smoothed state measurements to construct kernel
- Focuses on denoising through GP smoothing

**Your Van der Pol**:
- Can sample multiple GP curves (`num_gp_samples`)
- Combines original sparse data with GP samples via `'append'` or `'replace'` methods
- More explicit handling of duplicate time points

**Your Lotka-Volterra**:
- Uses mix of GP mean and GP samples (interleaved)
- Replaces original noisy data at original time points with mean of replicates
- More complex combination strategy

**Difference**: ⚠️ **Different augmentation strategies** (yours is more flexible)

### 4. **Hyperparameter Selection**

**Hsin et al. (2025)**:
- Uses cross-validation for λ selection in LASSO
- Tests multiple λ values and selects best fit
- More systematic hyperparameter tuning

**Your Implementation**:
- Uses data-adaptive threshold: `λ = λ_base * std(dXdt)` where `λ_base = 0.01`
- Fixed `max_iterations = 10` for STLS
- Simpler, heuristic-based approach

**Difference**: ⚠️ **Hsin et al. uses cross-validation, yours uses heuristics**

### 5. **Library Expansion**

**Hsin et al. (2025)**:
- Library size not explicitly stated, but appears to use standard polynomial library
- Focus on robustness rather than library expansion

**Your Lotka-Volterra**:
- **Explicitly expands library** when using GP-augmented data:
  - Basic: `[1, x, y, x², y², xy]` (6 terms)
  - Expanded: `[1, x, y, x², y², xy, x³, y³, x²y, xy²]` (10 terms)
- Rationale: More data allows more complex models

**Your Van der Pol**:
- Uses fixed polynomial order 3 library
- No explicit expansion based on data availability

**Difference**: ⚠️ **Lotka-Volterra explicitly expands library, Hsin et al. doesn't mention this**

### 6. **Ensemble Methods**

**Hsin et al. (2025)**:
- Does not mention ensemble methods
- Single model identification

**Your Lotka-Volterra**:
- Uses **ESINDy (Ensemble SINDy)** approach
- Runs 1000 ensemble members with different library subsets
- More robust to noise through ensemble averaging

**Your Van der Pol**:
- Uses single SINDy run (no ensemble)
- Simpler, faster approach

**Difference**: ⚠️ **Lotka-Volterra uses ensemble, Hsin et al. and Van der Pol don't**

### 7. **Time Formulation**

**Hsin et al. (2025)**:
- Continuous-time formulation: `ẋ = f(x, u)`
- Uses derivatives directly

**Your Van der Pol**:
- Continuous-time formulation: `dX/dt = f(X)`
- Uses derivatives directly

**Your Lotka-Volterra (Path 1)**:
- **Discrete-time formulation**: `X(t+1) = f(X(t))`
- No derivatives needed
- Different from Hsin et al. approach

**Difference**: ⚠️ **Lotka-Volterra Path 1 uses discrete-time, which is fundamentally different**

---

## Summary Table

| Aspect | Hsin et al. (2025) | Your Van der Pol | Your Lotka-Volterra |
|--------|-------------------|------------------|-------------------|
| **GP Kernel** | Squared Exponential | Squared Exponential | Squared Exponential |
| **Derivative Method** | GP Analytical | Hybrid (FD + GP) | GP Analytical |
| **Sparsity Optimization** | ADMM (LASSO) | STLS | STLS |
| **Data Augmentation** | GP Posterior Sampling | Append/Replace | Mean + Samples Mix |
| **Library** | Polynomial (order?) | Polynomial (order 3) | Expanded (3rd order) |
| **Ensemble** | No | No | Yes (ESINDy) |
| **Time Formulation** | Continuous | Continuous | Discrete (Path 1) / Continuous (Path 2) |
| **λ Selection** | Cross-validation | Data-adaptive | Fixed/Heuristic |

---

## Recommendations

### 1. **Consider ADMM for Sparsity Optimization**
Your STLS approach works, but ADMM (as in Hsin et al.) is more principled for L1 regularization and may provide better results with noisy data.

### 2. **Standardize Derivative Combination**
Your Van der Pol hybrid approach (FD + GP) is interesting but different from Hsin et al. Consider:
- **Option A**: Use only GP analytical derivatives (like Hsin et al. and your Lotka-Volterra Path 2)
- **Option B**: Keep hybrid but document rationale (e.g., preserving original data fidelity)

### 3. **Add Cross-Validation for λ**
Hsin et al.'s cross-validation approach for λ selection is more systematic than your data-adaptive heuristic.

### 4. **Consider Ensemble for Van der Pol**
Your Lotka-Volterra ESINDy ensemble approach shows robustness. Consider applying similar ensemble methods to Van der Pol for better noise handling.

### 5. **Document Library Expansion Rationale**
Your Lotka-Volterra library expansion is a good idea when you have more data. Consider making this explicit in Van der Pol as well.

---

## Code-Specific Observations

### Van der Pol Implementation
- ✅ Clean separation of phases (Phase 1-4)
- ✅ Good documentation and visualization
- ⚠️ Hybrid derivative approach differs from paper
- ⚠️ Could benefit from ensemble methods

### Lotka-Volterra Implementation
- ✅ Explicit comparison of discrete vs continuous-time
- ✅ Ensemble methods (ESINDy)
- ✅ Library expansion strategy
- ⚠️ More complex code structure
- ⚠️ Path 1 (discrete-time) is fundamentally different from GPSINDy approach

---

## Conclusion

Your implementations are **conceptually aligned** with Hsin et al. (2025) in the core GPSINDy approach:
1. GP smoothing of sparse noisy data ✅
2. Analytical derivative computation ✅
3. SINDy for model discovery ✅

**Key differences** are in:
1. **Optimization method** (STLS vs ADMM)
2. **Derivative combination** (hybrid vs pure GP)
3. **Ensemble methods** (your Lotka-Volterra uses ESINDy)
4. **Library expansion** (your explicit expansion strategy)

These differences are not necessarily weaknesses—your hybrid derivative approach and ensemble methods may actually provide benefits in certain scenarios. However, aligning more closely with Hsin et al.'s ADMM approach and cross-validation could improve robustness.

