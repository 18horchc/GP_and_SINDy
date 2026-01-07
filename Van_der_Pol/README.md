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

### Phase 4: GP+SINDy Implementation (Planned)
- GP posterior sampling
- Data augmentation
- SINDy on enhanced dataset

### Phase 5: Comparison & Analysis (Planned)
- Metrics computation
- Statistical comparison
- Visualizations

### Phase 6: Refinement (Planned)
- Parameter sensitivity
- Ablation studies

## Files

- `vdp_data.mat` - Saved data from all phases (generated after running phase scripts)
  - Phase 1: Ground truth, sparse noisy data
  - Phase 2: GP models, predictions, and derivatives
  - Phase 3: SINDy coefficients, discovered model, and trajectory

## Results

Results and visualizations will be generated as each phase is completed.

