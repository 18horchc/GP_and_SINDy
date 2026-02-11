# GP SINDy+MCMC Microglia Experiments

GP experiments on the **SINDy+MCMC microglial cell dynamics** model (Amato & Arnold, 2025). Two states: M1, M2. Same structure as `gp_logistic_design`, `gp_LV_design`, and `gp_kinetics_design`. Parameters match **sindy_mcmc_microglia.m** and **all_toy_problems.m** (default MCMC posterior means, M0=[5;5], tspan=[0 50]).

## Streamlined workflow (recommended)

Use the parameterized runner (same interface as other GP design folders):

```matlab
% Run all 20 experiments (aggregation is a separate call)
run_all_SINDy_experiments();

% When ready, build the summary table from all existing result files
run('gp_SINDy_design/aggregate_SINDy_metrics.m');
```

**Files:**
- **run_all_SINDy_experiments.m** — Config-driven driver: runs experiments. Supports `config_filter` and `run_subset`.
- **run_SINDy_experiment.m** — Core runner: takes config `{exp_id, noise_pct, is_regular, is_auto}` and runs one experiment.
- **master_run_scripts.m** — Convenience: calls `run_all_SINDy_experiments()`.
- **aggregate_SINDy_metrics.m** — Call separately; compiles all `results_SINDy_exp_*.mat` into `SINDy_metrics_all.csv` / `.mat`.

### Config filter examples (run subsets)

```matlab
run_all_SINDy_experiments(struct('config_filter', @(c) c.noise_pct == 0));  % zero noise only
run_all_SINDy_experiments(struct('config_filter', @(c) c.is_regular));      % regular sampling only
run_all_SINDy_experiments(struct('config_filter', @(c) ~c.is_auto));        % default optimization only
```

**SINDy-specific:** Two states (M1, M2 microglia). Five built-in kernels only.

---

## Legacy: Individual experiment scripts

The original `SINDy_exp_NN_<scope>.m` files (01–20) are in the `old/` folder for reference.

## Shared code (in this folder)

- **sindy_mcmc_microglia.m** — SINDy+MCMC ODE model (same as used by all_toy_problems).
- **ground_truth_SINDy.m** — High-res M1, M2 curves (500 pts on [0, 50]).
- **metric_helpers.m** — Lives in the **project root**. Runners use `addpath('..')` to access it.

## Experiments (all 20)

| Exp | Sampling | Noise | Optimization |
|-----|----------|-------|--------------|
| 01–05 | Regular | 0%, 1%, 5%, 10%, 20% | Default |
| 06–10 | Irregular | 0%, 1%, 5%, 10%, 20% | Default |
| 11–15 | Regular | 0%, 1%, 5%, 10%, 20% | Auto |
| 16–20 | Irregular | 0%, 1%, 5%, 10%, 20% | Auto |

Each varies N = 5, 10, 25, 50 and kernels: SqExp, Matern 1/2, 3/2, 5/2, Rational Quadratic. One GP per state (M1, M2).

## How to run

From the **project root** or this folder:

```matlab
cd('gp_SINDy_design')
run_all_SINDy_experiments
% When done, build the summary table:
run('aggregate_SINDy_metrics.m')
```

Or run one experiment only:

```matlab
cfg = struct('exp_id', 1, 'noise_pct', 0, 'is_regular', true, 'is_auto', false);
[T_point, T_prob, T_calib] = run_SINDy_experiment(cfg);
```
