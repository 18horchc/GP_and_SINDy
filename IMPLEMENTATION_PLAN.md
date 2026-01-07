# GP+SINDy Implementation Plan

## Overview
This document outlines the implementation plan for testing whether combining Gaussian Processes (GPs) and SINDy can improve data-driven model accuracy on sparse, noisy datasets.

---

## Stage 1: Toy Problem Selection

### Criteria for Good Toy Problem
1. **Known analytical solution** - We need ground truth to evaluate accuracy
2. **Sparse data challenge** - Problem should be difficult with few points
3. **Non-trivial dynamics** - Should test GP's ability to capture structure
4. **Derivative-dependent** - SINDy needs good derivative estimates
5. **Computational tractability** - Should run quickly for iteration

### Recommended Toy Problems (Ranked)

#### Option 1: **Van der Pol Oscillator** ⭐ (Recommended)
- **Equation**: `ẍ - μ(1-x²)ẋ + x = 0` (2D system: `ẋ = y`, `ẏ = μ(1-x²)y - x`)
- **Why it's good**:
  - Well-known, well-studied
  - Oscillatory behavior tests GP's ability to capture periodic structure
  - Nonlinear coupling between states
  - Can be made challenging with sparse sampling
- **Parameters**: μ = 1.0 (moderate nonlinearity)
- **Time span**: t ∈ [0, 10] with ~5-10 sparse points
- **Initial conditions**: x(0) = 2, y(0) = 0

#### Option 2: **Duffing Oscillator**
- **Equation**: `ẍ + δẋ + αx + βx³ = γcos(ωt)`
- **Why it's good**:
  - Forced nonlinear oscillator
  - Can have chaotic regimes
  - Tests GP on driven systems
- **Complexity**: More complex than Van der Pol

#### Option 3: **Lorenz System** (3D)
- **Equations**: Standard Lorenz attractor
- **Why it's good**:
  - Classic chaotic system
  - Tests multi-dimensional GP
- **Complexity**: 3D makes it more challenging

#### Option 4: **Lotka-Volterra** (You already have this!)
- **Current status**: You have implementations
- **Consideration**: Good baseline, but you may want a different problem to test generalizability

### Recommendation
Start with **Van der Pol** - it's a good balance of complexity and tractability, and tests oscillatory dynamics that GPs handle well.

---

## Stage 2: Data Generation

### 2.1 Ground Truth Generation
- Generate high-resolution solution using `ode45` or `ode15s`
- Time span: Sufficient to capture key dynamics (e.g., 2-3 periods for oscillators)
- Resolution: Fine enough for smooth derivatives (dt = 0.01-0.05)

### 2.2 Sparse Noisy Data Creation
- **Sparsity level**: 5-15 points (start with 8-10)
- **Sampling strategy**: 
  - Option A: Uniform random sampling
  - Option B: Strategic sampling (more points where dynamics change rapidly)
  - Option C: Uniform time spacing (simplest, good baseline)
- **Noise model**: 
  - Additive Gaussian noise: `y_noisy = y_true + ε`, where `ε ~ N(0, σ²)`
  - Noise level: `σ = noise_fraction * std(y_true)` where `noise_fraction ∈ [0.05, 0.20]`
  - Start with 10% noise (moderate challenge)

### 2.3 Data Validation
- Visualize: True solution + sparse noisy points
- Check: Noise level is visible but signal is recoverable
- Document: Number of points, noise level, time span

---

## Stage 3: GP-Only Implementation

### 3.1 GP Fitting Best Practices

#### Kernel Selection
- **Primary**: Squared Exponential (RBF) - good default
- **Alternatives to test**:
  - Matern 3/2 or 5/2 (less smooth, more flexible)
  - Periodic kernel (if oscillatory behavior is known)
  - Composite kernels (e.g., RBF + Periodic)

#### Hyperparameter Optimization
- **Method**: Maximum Likelihood Estimation (MLE) via `fitrgp`
- **Initialization**: 
  - Length scale: ~10% of time span
  - Signal variance: ~variance of data
  - Noise variance: Start with estimated noise level
- **Constraints**: Set reasonable bounds to avoid overfitting
- **Validation**: Check log-likelihood and R²

#### Basis Function
- **Options**: Constant, Linear, Quadratic
- **Recommendation**: Start with constant, try linear if needed
- **Note**: Your code uses constant - this is fine for most cases

#### Standardization
- **Decision**: `Standardize = true` for numerical stability
- **Exception**: If you need exact derivative formulas, may need `false`

### 3.2 GP Prediction
- **Time grid**: Dense grid (e.g., 100-200 points) for smooth interpolation
- **Output**: Mean prediction `μ(t)` and uncertainty `σ(t)`
- **Quality checks**:
  - R² between GP mean and true solution (on dense grid)
  - Prediction intervals should contain most true data
  - Visual inspection: GP should capture general trend

### 3.3 GP Derivative Computation
- **Method**: Analytical derivative (you already have this!)
- **Formula**: `dμ/dt = sum(αᵢ * dK/dt)` where K is kernel
- **Validation**: 
  - Compare GP derivative to true derivative (on dense grid)
  - Compute RMSE of derivative estimates
  - Visualize: True derivative vs GP derivative

### 3.4 Metrics for GP-Only
- **Fit quality**:
  - R² (coefficient of determination)
  - RMSE (root mean squared error)
  - Log-likelihood
- **Derivative quality**:
  - RMSE of derivatives
  - Correlation coefficient between true and GP derivatives

---

## Stage 4: SINDy-Only Implementation

### 4.1 SINDy Setup

#### Library Construction
- **Polynomial order**: Start with order 2 (cubic can overfit)
- **Library functions**: `[1, x, y, x², y², xy, ...]`
- **For Van der Pol**: Should include `[1, x, y, x², y², xy, x³, x²y, xy², y³]` to capture `μ(1-x²)y` term

#### Derivative Estimation (SINDy-only path)
- **Method**: Finite differences on sparse data
  - Central difference for interior points
  - Forward/backward difference for boundaries
- **Challenge**: Noisy, sparse data → poor derivative estimates
- **This is the key difficulty SINDy faces alone!**

### 4.2 SINDy Algorithm
- **Method**: STLS (Sequential Thresholded Least Squares) - you have this
- **Threshold parameter (λ)**:
  - Start with `λ = 0.01 * std(X_dot)` (data-adaptive)
  - Cross-validation: Try multiple λ values
  - Goal: Balance sparsity vs accuracy
- **Iterations**: 10-20 iterations of thresholding (you use 10-100)

### 4.3 Model Validation
- **Forward integration**: Use discovered model to predict trajectory
- **Initial conditions**: Use first data point (or true IC if known)
- **Time span**: Same as training data
- **Solver**: `ode45` or `ode15s` (stiff if needed)

### 4.4 Metrics for SINDy-Only
- **Coefficient accuracy**:
  - Compare discovered coefficients to true coefficients
  - Compute coefficient error (L2 norm of difference)
  - Identify which terms were correctly/incorrectly identified
- **Trajectory accuracy**:
  - RMSE between predicted and true trajectories
  - R² of predicted vs true
  - Long-term prediction error (if applicable)
- **Sparsity**: Number of non-zero coefficients (should match true model)

---

## Stage 5: GP+SINDy Combined Implementation

### 5.1 GP Posterior Sampling Strategy

#### Sampling Approach
- **Method**: Sample full curves from GP posterior (you have `sample_gp_posterior.m`)
- **Number of samples**: 
  - Option A: Single sample per run (simpler)
  - Option B: Multiple samples (e.g., 5-10), then aggregate SINDy results
  - **Recommendation**: Start with 1, then try ensemble of 5-10
- **Sampling points**: Dense grid (e.g., 50-100 points)

#### Data Combination Strategy
- **Option 1**: Replace sparse data entirely with GP samples
  - Risk: Loses information from original data
- **Option 2**: Combine original + GP samples
  - Weight original data more heavily (e.g., duplicate original points)
  - **Recommendation**: Use this approach
- **Option 3**: Use GP mean + add noise based on GP uncertainty
  - More conservative, respects GP uncertainty

### 5.2 Enhanced Data Set
- **Original data**: Keep all N sparse points
- **GP-sampled data**: Add M additional points (e.g., M = 50-100)
- **Total**: N + M points for SINDy training
- **Time grid**: Uniform or adaptive (more points where uncertainty is high)

### 5.3 SINDy on Enhanced Data
- **Derivatives**: Use GP analytical derivatives (not finite differences!)
- **Library**: Same as SINDy-only (for fair comparison)
- **Threshold**: Same λ as SINDy-only (for fair comparison)
- **Note**: More data should help, but need to avoid overfitting

### 5.4 Metrics for GP+SINDy
- **Same as SINDy-only**:
  - Coefficient accuracy
  - Trajectory accuracy
  - Sparsity
- **Additional**:
  - Improvement over SINDy-only (percentage reduction in error)
  - Comparison of derivative quality (GP vs finite difference)

---

## Stage 6: Comparison Framework

### 6.1 Metrics Summary Table

| Metric | GP-Only | SINDy-Only | GP+SINDy | True Model |
|--------|---------|------------|----------|------------|
| **Fit Quality** |
| R² (trajectory) | ✓ | ✓ | ✓ | 1.0 |
| RMSE (trajectory) | ✓ | ✓ | ✓ | 0.0 |
| **Derivative Quality** |
| RMSE (derivatives) | ✓ | N/A | ✓ | 0.0 |
| **Model Discovery** |
| Coefficient error | N/A | ✓ | ✓ | 0.0 |
| Terms identified | N/A | ✓ | ✓ | ✓ |
| **Prediction** |
| Long-term RMSE | ✓ | ✓ | ✓ | 0.0 |

### 6.2 Statistical Comparison
- **Multiple runs**: Run each method 10-20 times with different noise realizations
- **Report**: Mean ± std for each metric
- **Statistical tests**: 
  - Paired t-test: GP+SINDy vs SINDy-only
  - Effect size: How much improvement?

### 6.3 Visualization Plan

#### Figure 1: Data Overview
- True solution (smooth curve)
- Sparse noisy data points (scatter)
- Noise level annotation

#### Figure 2: GP Fit Quality
- True solution
- GP mean prediction
- GP uncertainty bands (2σ)
- Sparse data points
- Subplot: GP derivative vs true derivative

#### Figure 3: SINDy-Only Results
- True solution
- SINDy-discovered model prediction
- Sparse data points
- Subplot: Discovered coefficients vs true coefficients (bar chart)

#### Figure 4: GP+SINDy Results
- True solution
- GP+SINDy-discovered model prediction
- GP-sampled data (different marker/style)
- Original sparse data
- Subplot: Discovered coefficients vs true coefficients

#### Figure 5: Comparison Summary
- Side-by-side trajectories (True, SINDy-only, GP+SINDy)
- Metrics table/bar chart
- Coefficient comparison (all three)

#### Figure 6: Ablation Study (Optional)
- Vary number of GP samples (10, 50, 100, 200)
- Show how SINDy accuracy improves with more samples
- Identify "sweet spot" for sample size

---

## Stage 7: Best Practices & Overfitting Prevention

### 7.1 GP Overfitting Prevention
- **Hyperparameter constraints**: Set reasonable bounds
- **Cross-validation**: 
  - Leave-one-out or k-fold on sparse data
  - Check if GP generalizes to held-out points
- **Regularization**: Noise variance prevents overfitting to noise
- **Visual inspection**: GP should be smooth, not wiggly

### 7.2 SINDy Overfitting Prevention
- **Library size**: Don't use unnecessarily large libraries
- **Threshold selection**: 
  - Cross-validation on held-out data
  - Information criteria (AIC, BIC)
  - Stability: Model should be stable across noise realizations
- **Ensemble methods**: Your ESINDy helps here!
- **Validation**: Check if discovered model generalizes beyond training time

### 7.3 GP+SINDy Overfitting Prevention
- **GP quality check**: If GP fit is poor (R² < 0.5), don't trust GP samples
- **Sample quality**: GP samples should respect original data (pass through or near original points)
- **Data balance**: Don't overwhelm original data with too many GP samples
- **Independent validation**: Test discovered models on new time intervals

### 7.4 Parameter Selection Workflow
1. **GP hyperparameters**: Optimize via MLE, validate with R²
2. **SINDy threshold (λ)**: 
   - Try grid: `[0.001, 0.01, 0.1, 1.0] * std(X_dot)`
   - Choose λ that gives best coefficient accuracy (if known) or best trajectory fit
3. **Number of GP samples**: Start with 50, increase if needed
4. **Library order**: Start with 2, increase only if necessary

---

## Stage 8: Implementation Checklist

### Phase 1: Setup & Data Generation
- [ ] Choose toy problem (recommend Van der Pol)
- [ ] Implement ground truth ODE solver
- [ ] Generate sparse noisy data function
- [ ] Create visualization of data

### Phase 2: GP Implementation
- [ ] Set up GP fitting with proper hyperparameters
- [ ] Implement GP prediction on dense grid
- [ ] Compute GP analytical derivatives
- [ ] Validate GP quality (R², visual inspection)
- [ ] Create GP visualization

### Phase 3: SINDy-Only Implementation
- [ ] Implement finite difference derivatives
- [ ] Set up SINDy library construction
- [ ] Implement STLS algorithm (you have this)
- [ ] Forward integrate discovered model
- [ ] Compute SINDy-only metrics
- [ ] Create SINDy visualization

### Phase 4: GP+SINDy Implementation
- [ ] Implement GP posterior sampling (you have this)
- [ ] Combine original + GP-sampled data
- [ ] Run SINDy on enhanced dataset
- [ ] Compute GP+SINDy metrics
- [ ] Create combined visualization

### Phase 5: Comparison & Analysis
- [ ] Implement metrics computation for all methods
- [ ] Create comparison visualizations
- [ ] Run multiple trials (different noise realizations)
- [ ] Statistical analysis (mean, std, significance tests)
- [ ] Generate summary report

### Phase 6: Refinement
- [ ] Parameter sensitivity analysis
- [ ] Ablation studies (vary sample size, noise level, sparsity)
- [ ] Documentation and code comments
- [ ] Final visualizations and summary

---

## Stage 9: Key Decisions to Make

### 9.1 Toy Problem
- **Decision**: Van der Pol vs Lotka-Volterra vs Other
- **Rationale**: Van der Pol is good test, but you already have LV working

### 9.2 Sparsity Level
- **Decision**: How many sparse points? (5, 8, 10, 15?)
- **Rationale**: Too few = impossible, too many = not challenging
- **Recommendation**: Start with 8-10 points

### 9.3 Noise Level
- **Decision**: What noise fraction? (5%, 10%, 15%, 20%?)
- **Rationale**: Need challenge but recoverable signal
- **Recommendation**: Start with 10%

### 9.4 GP Sampling Strategy
- **Decision**: Single sample vs ensemble of samples?
- **Rationale**: Ensemble may be more robust but more complex
- **Recommendation**: Start with single, then try ensemble

### 9.5 SINDy Threshold Selection
- **Decision**: Fixed λ vs cross-validated λ?
- **Rationale**: CV is better but more work
- **Recommendation**: Start with fixed (data-adaptive), then try CV

### 9.6 Evaluation Focus
- **Decision**: Primary metric? (Coefficient accuracy vs trajectory accuracy)
- **Rationale**: Both matter, but coefficient accuracy is more direct
- **Recommendation**: Report both, but emphasize coefficient accuracy

---

## Stage 10: Potential Challenges & Solutions

### Challenge 1: GP Fits Poorly on Sparse Data
- **Symptom**: Low R², straight line predictions
- **Solutions**: 
  - Try different kernels (Matern, Periodic)
  - Adjust hyperparameter initialization
  - Use linear/quadratic basis functions
  - Consider sparse GP methods

### Challenge 2: GP Samples Don't Respect Original Data
- **Symptom**: GP samples far from original points
- **Solutions**:
  - Ensure GP passes through or near original data
  - Use conditional sampling that fixes original points
  - Weight original data more heavily

### Challenge 3: SINDy Discovers Wrong Model
- **Symptom**: Wrong coefficients, missing terms
- **Solutions**:
  - Adjust threshold λ
  - Try different library (add/remove terms)
  - Use ensemble SINDy (you have ESINDy!)
  - Check derivative quality (may be the issue)

### Challenge 4: GP+SINDy Not Better Than SINDy-Only
- **Symptom**: No improvement or worse performance
- **Solutions**:
  - Check GP quality (may be poor)
  - Check if GP samples are adding noise
  - Try different sampling strategies
  - May need more GP samples

### Challenge 5: Overfitting
- **Symptom**: Good on training, poor on validation
- **Solutions**:
  - Reduce library size
  - Increase SINDy threshold
  - Use cross-validation
  - Check GP hyperparameters (may be overfit)

---

## Next Steps

1. **Review this plan** - Does it align with your goals?
2. **Choose toy problem** - Van der Pol or stick with Lotka-Volterra?
3. **Clarify decisions** - Make choices on key parameters
4. **Start implementation** - Begin with Phase 1 (data generation)
5. **Iterate** - Test each phase before moving to next

---

## Questions for Discussion

1. Do you want to use Van der Pol or stick with Lotka-Volterra?
2. How sparse should the data be? (I recommend 8-10 points)
3. What noise level? (I recommend 10%)
4. Single GP sample or ensemble? (Start with single, then try ensemble)
5. Primary evaluation metric? (Coefficient accuracy vs trajectory)
6. Do you want to test multiple toy problems or focus on one?

---

## References & Resources

- **GP**: Rasmussen & Williams, "Gaussian Processes for Machine Learning"
- **SINDy**: Brunton et al., "Discovering governing equations from data"
- **Your existing code**: Good foundation to build on!

