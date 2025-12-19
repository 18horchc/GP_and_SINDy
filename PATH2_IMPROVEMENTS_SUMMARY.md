# Path 2 Improvements Summary

## Changes Implemented

### 1. ✅ GP Sampling: One Curve Per Run
- **Status**: Implemented
- **Details**: `num_gp_samples = 1` - samples one GP curve per run
- **Location**: Line ~238 in `lotka_volterra_comparison.m`

### 2. ✅ Composite Kernel: Rational Quadratic
- **Status**: Implemented (using Rational Quadratic as built-in flexible kernel)
- **Details**: 
  - Switched from Matern32 to Rational Quadratic kernel
  - Rational Quadratic is more flexible and can capture multiple scales
  - **STABILITY ISSUE FLAGGED**: Rational Quadratic derivative computation may have issues
  - X derivative has std=0.0002 (essentially zero), causing all-zero coefficients
  - Added finite differences fallback when derivative variance is too low
- **Location**: Lines ~199-211 in `lotka_volterra_comparison.m`
- **Derivative Support**: Added to `gp_derivative_analytical.m` (lines ~72-90)

### 3. ✅ Data Usage: GP Mean + Samples
- **Status**: Implemented
- **Details**: 
  - Interleaves GP mean (even indices) and GP samples (odd indices)
  - This gives a mix of smooth mean estimates and variability
  - Original data points (mean of replicates) are preserved at their time points
- **Location**: Lines ~271-300 in `lotka_volterra_comparison.m`

## Stability Issues Flagged

### ⚠️ CRITICAL: Rational Quadratic Derivative Issue
- **Problem**: X derivative has std=0.0002 (essentially zero)
- **Impact**: All Path 2 coefficients for X are zero
- **Cause**: Likely incorrect derivative formula for Rational Quadratic kernel
- **Mitigation**: Added finite differences fallback (should trigger automatically)
- **Status**: Needs investigation - derivative formula may need correction

### ⚠️ Potential Issues:
1. **Composite Kernel Derivatives**: Custom composite kernels would require more complex derivative formulas
2. **Numerical Stability**: Rational Quadratic derivatives may be numerically unstable
3. **GP Fit Quality**: R² is good (0.98), but derivatives may not reflect this

## Recommendations

1. **For Rational Quadratic**: Verify derivative formula or use finite differences
2. **For True Composite Kernels**: Would need custom kernel function + custom derivative implementation
3. **Alternative**: Try Matern52 or other built-in kernels that might work better

## Next Steps

1. Fix Rational Quadratic derivative computation
2. Test if finite differences fallback resolves the zero coefficients issue
3. Consider trying other kernel options if RQ continues to have issues

