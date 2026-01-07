# Van der Pol Oscillator: GP+SINDy Comparison

This folder contains the implementation for testing whether combining Gaussian Processes (GPs) and SINDy can improve data-driven model accuracy on sparse, noisy datasets using the Van der Pol oscillator as a toy problem.

## Problem Setup

**Van der Pol Oscillator:**
- Equations: `dx/dt = y`, `dy/dt = μ(1-x²)y - x`
- Parameters: μ = 1.0 (moderate nonlinearity)
- Initial conditions: x(0) = 2, y(0) = 0
- Time span: t ∈ [0, 10]

**Data Characteristics:**
- Sparse data: 8-10 points (default: 9)
- Noise level: 10% Gaussian noise
- Sampling: Uniform time spacing

## Implementation Phases

### Phase 1: Data Generation ✅
- `generate_vdp_ground_truth.m` - Generate true Van der Pol solution
- `generate_sparse_noisy_data.m` - Create sparse noisy observations
- `phase1_data_generation.m` - Main script for Phase 1

**Usage:**
```matlab
cd Van_der_Pol
phase1_data_generation
```

### Phase 2: GP Implementation ✅
- `fit_gp_vdp.m` - Fit GP models and compute analytical derivatives
- `phase2_gp_implementation.m` - Main script for Phase 2
- GP fitting on sparse noisy data
- GP prediction on dense time grid
- Analytical derivative computation
- GP quality validation (R² metrics)

**Usage:**
```matlab
cd Van_der_Pol
% Make sure Phase 1 has been run first
phase2_gp_implementation
```

### Phase 3: SINDy-Only Implementation ✅
- `compute_finite_differences.m` - Compute derivatives using finite differences
- `build_library_vdp.m` - Build polynomial library matrix (up to order 3)
- `run_sindy_stls.m` - Run Sequential Thresholded Least Squares algorithm
- `forward_integrate_sindy.m` - Forward integrate discovered model
- `phase3_sindy_only.m` - Main script for Phase 3
- Finite difference derivatives from sparse noisy data
- STLS algorithm for sparse model discovery
- Forward integration and trajectory comparison

**Usage:**
```matlab
cd Van_der_Pol
% Make sure Phase 1 has been run first
phase3_sindy_only
```

### Phase 4: GP+SINDy Implementation ✅
- `augment_data_with_gp_samples.m` - Sample from GP posterior and augment sparse data
- `phase4_gp_sindy.m` - Main script for Phase 4
- GP posterior sampling to generate additional data points
- Data augmentation (combining sparse + GP samples)
- SINDy on enhanced dataset with improved derivatives
- Forward integration and comparison with SINDy-only

**Usage:**
```matlab
cd Van_der_Pol
% Make sure Phases 1, 2, and 3 have been run first
phase4_gp_sindy
```

### Phase 5: Comparison & Analysis ✅
- `compute_true_vdp_coefficients.m` - Compute true Van der Pol coefficients for comparison
- `compute_comprehensive_metrics.m` - Compute comprehensive metrics for all methods
- `phase5_comparison_analysis.m` - Main script for Phase 5
- Comprehensive metrics computation (trajectory, derivative, coefficient metrics)
- Statistical comparison between methods
- Improvement analysis (GP+SINDy vs SINDy-only)
- Summary visualizations and tables

**Usage:**
```matlab
cd Van_der_Pol
% Make sure Phases 1-4 have been run first
phase5_comparison_analysis
```

### Phase 6: Refinement & Analysis ✅
- `phase6_parameter_sensitivity.m` - Parameter sensitivity analysis
- `phase6_statistical_analysis.m` - Statistical robustness analysis (50 trials)
- Parameter sensitivity: Vary sparse points (5, 7, 9, 11, 15) and noise levels (5%, 10%, 15%, 20%)
- Statistical analysis: 50 noise realizations with fixed parameters
- Comprehensive visualizations: heatmaps, distributions, box plots
- Summary statistics and improvement analysis

**Usage:**
```matlab
cd Van_der_Pol
% Make sure all phase functions are available
% Parameter sensitivity analysis
phase6_parameter_sensitivity

% Statistical robustness analysis
phase6_statistical_analysis
```

**Note:** Phase 6 does NOT modify previous phase code. All refinements use separate files.

## Files

- `vdp_data.mat` - Saved data from all phases (generated after running phase scripts)
  - Phase 1: Ground truth, sparse noisy data
  - Phase 2: GP models, predictions, and derivatives
  - Phase 3: SINDy coefficients, discovered model, and trajectory
  - Phase 4: Augmented data, GP+SINDy coefficients, and trajectory
  - Phase 5: Comprehensive metrics and comparison results

## Results

Results and visualizations will be generated as each phase is completed.

