# MATLAB Implementation Plan: GP Experimental Design for Logistic Growth

**Source:** Experimental Design for GPs.pdf  
**Scope:** Logistic growth curve only; full factorial design with 400 trial types × K realizations.

---

## 1. Design Summary (from PDF)

| Factor | Levels | Values |
|--------|--------|--------|
| **A: Kernel** | 5 | Squared Exponential, Matern 1/2, 3/2, 5/2, Rational Quadratic |
| **B: Optimization** | 2 | Default vs. OptimizeHyperparameters: auto |
| **C: Noise (σ)** | 5 | 0%, 1%, 5%, 10%, 20% |
| **D: Sampling Mode** | 2 | Regular (fixed intervals) vs. Irregular (random times) |
| **E: Sparsity (N)** | 4 | 5, 10, 25, 50 points |

- **Total trial types:** 5 × 2 × 5 × 2 × 4 = **400**
- **Realizations per type:** K = 20 (different random seeds)
- **Total GP fits:** 400 × 20 = **8,000**

---

## 2. Implementation Steps

### 2.1 Ground truth (Step 1)

- Use existing [logistic_growth.m](logistic_growth.m) with `plot_results = false`.
- Generate **500-point** high-resolution curve over **t = [0, 50]**:
  - `params.tspan = [0 50]`; call `ode45` with output time vector `linspace(0, 50, 500)'` (or get 500 points via interpolation from ode45 output).
- Store as vectors `t_gt` (500×1) and `y_gt` (500×1). This is the “true” f(t) for all metrics.

### 2.2 Sampling logic (Step 2)

- **Regular (N):** `t_obs = linspace(0, 50, N)'` (N evenly spaced points including 0 and 50).
- **Irregular (N):** Draw N unique times from U(0, 50), e.g. `sort(rand(N,1)) * 50` or `randperm`-based sampling so endpoints can be included if desired (PDF says “uniform distribution U(0,50)”).
- **Noise:** For each selected t, evaluate ground truth (interpolate from `t_gt`, `y_gt`), then:
  - `y_sampled = y_true + eps`, where `eps ~ N(0, (noise_pct * Amplitude)^2)`.
  - **Amplitude:** PDF does not define it (see Section 4). Recommend: **Amplitude = std(y_gt)** so “10% noise” = 0.1×std(y_gt). Alternative: max(y_gt) or range(y_gt); must be fixed once and documented.
- Use `rng(seed)` before drawing irregular times and before adding noise so each realization is reproducible.

### 2.3 Design matrix and evaluation loop (Step 3)

- Build a **design matrix** (400 rows): each row = (Kernel, Optimization, Noise_pct, Sampling, N).
- For each row (trial type) and each realization k = 1,…,K:
  1. Set `rng(seed)` with a unique seed, e.g. `seed = (trial_index - 1)*K + k`.
  2. Generate sample: N times (regular or irregular) + noisy y.
  3. Fit GP: `fitrgp(t_obs, y_obs, ...)` with the row’s Kernel and Optimization.
  4. Predict on ground-truth grid: `[ymu, ystd, yint] = predict(gprMdl, t_gt)`.
  5. Compute all metrics (Section 2.4) and store in a table row.

- **Kernel mapping (MATLAB names):**
  - Squared Exponential → `'squaredexponential'`
  - Matern 1/2 → `'exponential'`
  - Matern 3/2 → `'matern32'`
  - Matern 5/2 → `'matern52'`
  - Rational Quadratic → `'rationalquadratic'`

- **Optimization levels:**
  - **Default:** Do not set `'OptimizeHyperparameters'` (use fitrgp’s default fitting of sigma/kernel scale).
  - **Auto:** `'OptimizeHyperparameters','auto'` with `HyperparameterOptimizationOptions` (e.g. `ShowPlots', false`, `Verbose`, 0) for reproducibility.

### 2.4 Metrics (Step 3)

For each trial, compare GP predictions on the **500-point ground-truth grid** (not just the N training points).

| Metric | Definition / MATLAB |
|--------|----------------------|
| **RMSE** | `sqrt(mean((ymu - y_gt).^2))` |
| **Max Error** | `max(abs(ymu - y_gt))` |
| **LML** | Log marginal likelihood of the GP model. **Issue:** fitrgp does not expose LML directly when not using `OptimizeHyperparameters`; when using `'auto'`, check `gprMdl.HyperparameterOptimizationResults` or compute LML from kernel + sigma (see Section 4). |
| **PICP** | Proportion of the 500 ground-truth points that lie inside the GP 95% prediction interval. Use `yint` from `predict(..., 'Alpha', 0.05)` or equivalent. `PICP = mean(y_gt >= yint(:,1) & y_gt <= yint(:,2))`. |
| **MPIW** | Mean width of the 95% interval: `mean(yint(:,2) - yint(:,1))`. |
| **CRPS** | Continuous Ranked Probability Score. For Gaussian predictive at each point, CRPS has a closed form; average over the 500 points. Formula: e.g. `sigma * (z*(2*normcdf(z)-1) + 2*normpdf(z) - 1/sqrt(pi))` with `z = (y_gt - ymu)./sigma`, then `mean(crps_i)`. |

### 2.5 Storage and analysis (Section 4 of PDF)

- **Table columns (per run):** TrialTypeId, Kernel, Optimization, NoisePct, Sampling, N, Realization, RMSE, MaxError, LML, PICP, MPIW, CRPS, (optional: runtime_sec).
- Save to **MAT** and/or **CSV** (e.g. `results_gp_logistic_factorial.mat` and `.csv`) so you can use `groupsummary` and `fitlm` in MATLAB.
- **Main effects:** For each factor, average metrics across all other factors; plot factor level vs mean RMSE (and optionally other metrics).
- **Interaction (Noise × Sampling):** Plot RMSE vs Noise with two lines (Regular vs Irregular); crossing or different slopes indicate interaction.
- **Kernel comparison:** Boxplot of RMSE (or other metric) by Kernel across all trials; or mean RMSE per kernel averaged over all conditions.

---

## 3. Suggested File Structure

| File | Purpose |
|------|--------|
| `gp_logistic_design/design_matrix.m` | Returns the 400×5 design matrix (or table) of factor levels. |
| `gp_logistic_design/ground_truth_logistic.m` | Returns (t_gt, y_gt) for t in [0,50], 500 points, using logistic_growth. |
| `gp_logistic_design/sample_observations.m` | Given (t_gt, y_gt), N, sampling type, noise_pct, Amplitude, seed → (t_obs, y_obs). |
| `gp_logistic_design/fit_and_metrics.m` | Given (t_obs, y_obs), (t_gt, y_gt), kernel, optimization_option → struct of RMSE, MaxError, LML, PICP, MPIW, CRPS. |
| `gp_logistic_design/run_one_trial.m` | Given one design row + realization index: set rng, sample, fit, compute metrics, return row. |
| `gp_logistic_design/run_full_factorial.m` | Loop over 400 types × K realizations; accumulate table; save MAT/CSV. |
| `gp_logistic_design/analyze_main_effects.m` | Load results, compute main effects, plot (e.g. RMSE vs Noise, vs N, vs Kernel). |
| `gp_logistic_design/analyze_interactions.m` | Noise × Sampling plot; optional other 2-way plots. |
| `gp_logistic_design/compute_CRPS.m` | CRPS for Gaussian predictive (vectorized over 500 points). |
| `gp_logistic_design/compute_LML.m` | Compute log marginal likelihood from fitted GP (if not in model object). |

---

## 4. Flagged Issues and Oversights

1. **“Amplitude” for noise (PDF Step 2)**  
   The PDF defines noise as ε ~ N(0, (Noise% × Amplitude)²) but does not define Amplitude.  
   - **Recommendation:** Define **Amplitude = std(y_gt)** and document it. Then 10% noise = 0.1×std(y_gt). Alternatively use max(y_gt) or range(y_gt); keep one definition consistent.

2. **Factor B: “Default” vs “OptimizeHyperparameters: auto”**  
   In MATLAB, fitrgp already fits hyperparameters (sigma, kernel scale) by default. “Default” may mean “no extra hyperparameter optimization,” while “auto” turns on Bayesian optimization over hyperparameters.  
   - **Recommendation:** Define explicitly: e.g. **Default** = no `OptimizeHyperparameters` argument (standard MLE fit); **Auto** = `'OptimizeHyperparameters','auto'` with fixed options. Document in code and in the results table.

3. **LML not directly available**  
   fitrgp does not expose the log marginal likelihood as a property when using the default fit. When using `'OptimizeHyperparameters','auto'`, the optimization result may contain a loss; that loss may be negative log likelihood.  
   - **Recommendation:** Implement a small helper that computes LML from the fitted model (using alpha, L from the internal representation, or by reconstructing K and σ² from kernel parameters and Sigma). Alternatively, document that LML is “NA” for Default and only stored for Auto from HyperparameterOptimizationResults if available.

4. **Endpoint inclusion for Irregular sampling**  
   PDF says “draw N unique time-steps from U(0, 50).” So 0 and 50 are not guaranteed.  
   - **Recommendation:** Use exactly U(0,50) as stated; if you want to include endpoints for comparability with Regular, document a variant (e.g. “irregular_with_endpoints”) and run it as a separate check rather than changing the main design.

5. **Total runtime**  
   8,000 GP fits can be slow (minutes to hours depending on N and optimizer).  
   - **Recommendation:** Run in batches (e.g. by trial type or by N), save intermediate results, and optionally reduce K for a quick test (e.g. K=5) before full K=20.

6. **Random seed scheme**  
   PDF says “rng(trial_number)” but each trial type has K realizations.  
   - **Recommendation:** Use a deterministic seed per (trial_type, realization), e.g. `rng((trial_type_index - 1)*K + realization_index)` so results are reproducible and each of the 8,000 runs has a unique seed.

7. **CRPS and prediction interval format**  
   fitrgp’s `predict` can return prediction intervals; ensure the interval is 95% (Alpha 0.05) for PICP and MPIW. For CRPS, use the predictive standard deviation (second output of predict) and the Gaussian CRPS formula.

8. **Possible failure of fitrgp**  
   Some kernel/optimizer/noise/sparsity combinations may lead to numerical issues or failed fits.  
   - **Recommendation:** Wrap each fit in try/catch; store NaN for metrics on failure and a “failed” flag in the table so you can count and optionally exclude or report failures.

---

## 5. Order of Implementation

1. Implement ground truth and sampling (Section 2.1–2.2).  
2. Implement metrics: RMSE, MaxError, PICP, MPIW, CRPS (Section 2.4); then LML (helper or from optimization results).  
3. Implement design matrix and single-trial function (Section 2.3).  
4. Run full factorial (optionally first with small K or subset of design).  
5. Implement analysis scripts (main effects, interaction, kernel boxplot).  
6. Document Amplitude, Default vs Auto, and seed scheme in a short README in the same folder.
