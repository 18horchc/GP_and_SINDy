# Bio-Informed GP Design

GP encoding experiments for biologically informed priors (acute transient and future toy problems).

## Quick Start

```matlab
addpath('gp_bio_informed_design');
addpath('..');

% Run with defaults: acute_transient, 10 irregular points, 2 replicates, 5% noise
[T, path] = run_bio_experiment('acute_transient');

% Custom data config
data_cfg.n_points = 15;
data_cfg.n_replicates = 3;
[T, path] = run_bio_experiment('acute_transient', data_cfg);

% Specific encodings only
[T, path] = run_bio_experiment('acute_transient', [], {'vanilla','mean_only','deriv_only'});

% All 16 encoding combinations
[T, path] = run_bio_experiment('acute_transient', [], 'all');
```

## Encoding Registry (Acute Transient)

Mix-and-match encodings via `encoding_registry`:

| Flag | Description |
|------|-------------|
| `informative_mean` | Pseudo-obs at t=0 and t=20 with y=b (0.1) to encode baseline |
| `kernel_matern32` | Use Matérn 3/2 kernel instead of Squared Exponential |
| `late_pseudo_obs` | Pseudo-obs at t=16,18,20 near baseline |
| `derivative_constraints` | Shape-guide points (prior-based) to encourage rise-then-fall |

**Predefined config IDs:** `vanilla`, `mean_only`, `matern_only`, `pseudo_only`, `deriv_only`, `mean_matern`, `mean_pseudo`, `mean_deriv`, `matern_pseudo`, `matern_deriv`, `pseudo_deriv`, `mean_matern_pseudo`, `mean_matern_deriv`, `mean_pseudo_deriv`, `matern_pseudo_deriv`, `full`.

**Custom combination:**
```matlab
flags.informative_mean = true;
flags.kernel_matern32 = false;
flags.late_pseudo_obs = true;
flags.derivative_constraints = false;
config = encoding_registry('custom', flags);
```

## Baseline

**Vanilla:** SE kernel, constant mean (fitted), no pseudo-observations, no derivative encoding. All encodings are compared against this on the same data.

## Data Config (Default)

- 10 irregularly sampled time points in [0, 20]
- 2 replicates per time point (= 20 observations)
- 5% proportional noise: σ = 0.05 × |y_true|

## Output

Results saved to `gp_bio_informed_design/results/`:
- `results_bio_acute_transient_N10_2rep_5noise.mat`
- `results_bio_acute_transient_N10_2rep_5noise.csv`

Metrics: RMSE, MAE, R², NLPD, MSLL, CRPS, sMSE, Coverage, NLML (from `metric_helpers`).

## Files

| File | Purpose |
|------|---------|
| `encoding_registry.m` | Mix-and-match encoding configs |
| `fit_gp_bio.m` | Fit GP with encoding (fitrgp-based) |
| `ground_truth_bio.m` | Ground truth wrapper for bio_informed_gp_priors |
| `generate_bio_data.m` | Generate (t,y) observations with noise |
| `run_bio_experiment.m` | Main experiment runner |

## GPML Fixed Prior Mean

For a **true fixed prior mean** m(x)=b (no pseudo-observations), use GPML:

```matlab
[T, gpr_naive, hyp_gpml] = run_naive_vs_gpml_fixed_mean();
% With GPML path:
[T, ~, ~] = run_naive_vs_gpml_fixed_mean(struct('gpml_path', 'C:\path\to\gpml-matlab'));
```

**Requires:** [GPML toolbox](https://gitlab.com/hnickisch/gpml-matlab). Download and unpack, then pass `opts.gpml_path` to the run script or `fit_gp_gpml_fixed_mean`. Uses `meanConst(5)` with fixed hyperparameter; only covariance and likelihood are optimized.

## Notes

- **Derivative constraints:** Implemented as shape-guide virtual points using the prior formula f(t)=b+A·t·exp(-ct). For true derivative observations (df/dt at points), GPML or custom inference would be needed.
- **Informative mean:** fitrgp: pseudo-observations. GPML: fixed `meanConst` (proper prior).
- Uses MATLAB `fitrgp` (Statistics and Machine Learning Toolbox); GPML optional.
