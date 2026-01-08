# GPSINDy Derivative Update Summary

## Overview
Updated both **Van der Pol** and **Lotka-Volterra** implementations to match the Hsin et al. (2025) approach: **using GP analytical derivatives for ALL augmented data points** instead of mixing finite differences with GP derivatives.

---

## Changes Made

### 1. Van der Pol Implementation (`Van_der_Pol/augment_data_with_gp_samples.m`)

#### Before:
- Used **finite differences** for original sparse time points
- Used **GP analytical derivatives** for dense GP-augmented time points
- Mixed approach: `dXdt_aug = [dXdt_sparse_fd; dXdt_gp_combined]`

#### After:
- Uses **GP analytical derivatives for ALL points** (both sparse and dense)
- Computes derivatives at all time points upfront: `t_all = [t_sparse; t_dense]` (for 'append' method)
- Consistent approach: `dXdt_aug = [dXdt_sparse_gp; dXdt_dense_gp]` where both are from GP

#### Key Changes:
1. **Lines 74-116**: Added computation of GP analytical derivatives for ALL time points (sparse + dense)
2. **Lines 159-163**: Updated to use GP derivatives for sparse points instead of finite differences
3. **Lines 174-175**: Fixed 'replace' case to use all GP derivatives correctly
4. **Updated documentation**: Added note about following Hsin et al. (2025) approach

#### Benefits:
- ✅ Consistent denoised derivatives throughout
- ✅ Better handling of noisy sparse data (GP smooths noise)
- ✅ Matches Hsin et al. (2025) methodology
- ✅ No mixing of different derivative computation methods

---

### 2. Lotka-Volterra Implementation (`lotka_volterra_comparison.m`)

#### Before:
- Computed GP derivatives only for dense grid (`X_dense`)
- Comments didn't explicitly state that derivatives are for all points

#### After:
- **Already computes GP derivatives for all points** (since `t_combined = t_dense`)
- Updated comments to explicitly state this matches Hsin et al. (2025) approach
- Clarified that GP analytical derivatives are used for all augmented data points

#### Key Changes:
1. **Lines 299-318**: Updated comments to explicitly state GP derivatives are computed for ALL points
2. **Lines 383-385**: Clarified that derivatives match the combined dataset structure

#### Note:
The Lotka-Volterra implementation was already mostly correct since `t_combined = t_dense`, meaning GP derivatives computed for `X_dense` automatically apply to all points in the combined dataset. The main update was clarifying documentation and ensuring consistency with Hsin et al. approach.

---

## Technical Details

### GP Analytical Derivative Computation

Both implementations use the same mathematical approach:

For squared exponential kernel: `k(t, t') = σ_f² * exp(-0.5 * (t-t')² / L²)`

**Derivative of kernel:**
```
dk/dt = -(t_star - t_train) / L² * k(t_star, t_train)
```

**GP mean derivative:**
```
dμ/dt = Σ(α_i * dk/dt)
```

Where:
- `α` = GP model's Alpha vector (from `gp_model.Alpha`)
- `L` = lengthscale (from kernel parameters)
- `σ_f` = signal standard deviation (from kernel parameters)

### Implementation Consistency

Both implementations now:
1. ✅ Compute GP analytical derivatives at **all** time points (sparse + dense)
2. ✅ Use GP derivatives consistently throughout (no finite differences)
3. ✅ Match Hsin et al. (2025) methodology
4. ✅ Provide denoised, uncertainty-quantified derivatives

---

## Testing Recommendations

To verify the changes work correctly:

### Van der Pol:
```matlab
% Run Phase 4
cd Van_der_Pol
phase4_gp_sindy
% Check that dXdt_aug uses GP derivatives for all points
% Verify no finite differences are used
```

### Lotka-Volterra:
```matlab
% Run comparison script
lotka_volterra_comparison
% Verify GP derivatives are computed for all points in t_combined
% Check that Path 2 uses GP analytical derivatives consistently
```

---

## Comparison with Hsin et al. (2025)

| Aspect | Before | After | Hsin et al. |
|--------|--------|-------|-------------|
| **Sparse points derivatives** | Finite differences | GP analytical | GP analytical |
| **Dense points derivatives** | GP analytical | GP analytical | GP analytical |
| **Consistency** | Mixed methods | Consistent GP | Consistent GP |
| **Noise handling** | Partial (sparse noisy) | Full GP smoothing | Full GP smoothing |

---

## Files Modified

1. `Van_der_Pol/augment_data_with_gp_samples.m` - Major update
2. `lotka_volterra_comparison.m` - Documentation/clarification update

---

## Next Steps (Optional)

Consider these additional improvements to further align with Hsin et al.:

1. **ADMM for Sparsity Optimization**: Replace STLS with ADMM (as in Hsin et al.)
2. **Cross-Validation for λ**: Use cross-validation for sparsity threshold selection
3. **Ensemble Methods**: Consider applying ensemble methods to Van der Pol (like Lotka-Volterra)

---

## Summary

Both implementations now use **GP analytical derivatives for all augmented data points**, matching the Hsin et al. (2025) GPSINDy approach. This provides:
- Consistent, denoised derivatives throughout
- Better noise handling for sparse data
- Alignment with published methodology

The changes maintain backward compatibility (same function signatures) while improving the derivative computation approach.

