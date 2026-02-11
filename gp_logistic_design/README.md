# GP Logistic Growth Experiments

Experiments for the **GP-on-logistic-growth** design (see [GP_Logistic_Experimental_Design_Plan.md](../GP_Logistic_Experimental_Design_Plan.md) in the project root).

## Streamlined workflow (recommended)

Instead of running 20 individual experiment files, use the parameterized runner:

```matlab
% Run all 60 experiments (aggregation is a separate call)
run_all_logistic_experiments();

% When ready, build the summary table from all existing result files
run('gp_logistic_design/aggregate_logistic_metrics.m');
```

**Files:**
- **run_all_logistic_experiments.m** — Config-driven driver: runs experiments. Supports `config_filter` and `run_subset`.
- **run_logistic_experiment.m** — Core runner: takes config `{exp_id, noise_pct, is_regular, is_auto}` and runs one experiment.
- **master_run_scripts.m** — Convenience: calls `run_all_logistic_experiments()`.
- **aggregate_logistic_metrics.m** — Call separately; compiles all `results_exp_*.mat` into `logistic_metrics_all.csv` / `.mat`.

**Benefits:** One place to fix bugs, add kernels, or change N_list. Add a new condition by extending `build_logistic_configs()`.

### Config filter examples (run subsets)

```matlab
% Run only base experiments (1 obs per time)
run_all_logistic_experiments(struct('config_filter', @(c) c.n_replicates == 1));

% Run only 3-replicate experiments
run_all_logistic_experiments(struct('config_filter', @(c) c.n_replicates == 3));
```

Use `config_filter` to run only experiments matching criteria. Aggregation is **separate** — it never overwrites until you explicitly call it.

```matlab
% Zero noise only (exp 01, 06, 11, 16)
run_all_logistic_experiments(struct('config_filter', @(c) c.noise_pct == 0));

% Regular sampling only (exp 01-05, 11-15)
run_all_logistic_experiments(struct('config_filter', @(c) c.is_regular));

% Default optimization only (exp 01-10)
run_all_logistic_experiments(struct('config_filter', @(c) ~c.is_auto));

% Low noise (0% and 1%)
run_all_logistic_experiments(struct('config_filter', @(c) c.noise_pct <= 0.01));

% Struct match: noise_pct == 0
run_all_logistic_experiments(struct('config_filter', struct('noise_pct', 0)));

% After running subsets, aggregate when ready (merges all existing .mat files)
run('gp_logistic_design/aggregate_logistic_metrics.m');
```

Index-based subset (by position in filtered list):

```matlab
run_all_logistic_experiments(struct('config_filter', @(c) c.is_regular, 'run_subset', [1 3 5]));
```

---

## Legacy: Individual experiment scripts

The original `log_exp_NN_<scope>.m` files (01–20) still work and produce identical outputs. They are kept for reference; new work can use the streamlined runner.

## Naming convention

- **Runner scripts:** `log_exp_NN_<scope>.m` (legacy) or `run_logistic_experiment(cfg)`.
- **Results:** `results_exp_NN_<scope>.mat` and `results_exp_NN_<scope>.csv` (same stem; e.g. `results_exp_01_regular_0noise`).
- **Shared helpers:** `ground_truth_logistic.m`, `logistic_growth.m`, `metric_helpers.m` (project root).

## Experiments (60 total)

| Exp | Sampling | Noise | Optimization | Replicates |
|-----|----------|-------|--------------|------------|
| 01–20 | Regular/irregular | 0–20% | Default/auto | 1 per time |
| 21–40 | Same | Same | Same | 3 per time |
| 41–60 | Same | Same | Same | 8 per time |

N = 5, 10, 25, 50 time points. Kernels: SqExp, Matern 1/2, 3/2, 5/2, Rational Quadratic. Replicate noise is proportional to the ground truth at each time (realistic for growth data).

## How to run

From the **project root** or this folder:

```matlab
cd('gp_logistic_design')
run_all_logistic_experiments
% When done, build the summary table:
run('aggregate_logistic_metrics.m')
```

Or run one experiment only:

```matlab
cfg = struct('exp_id', 1, 'noise_pct', 0, 'is_regular', true, 'is_auto', false, 'n_replicates', 1);
[T_point, T_prob, T_calib] = run_logistic_experiment(cfg);
```

## Shared code (in this folder)

- **logistic_growth.m** — Logistic ODE model (same as used by all_toy_problems).
- **ground_truth_logistic.m** — Returns high-res logistic curve (default 500 pts on [0, 50]).
- **metric_helpers.m** — Lives in the **project root** (`GP_and_SINDy/metric_helpers.m`) and is shared with `gp_LV_design`.
