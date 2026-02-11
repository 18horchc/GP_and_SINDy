# GP Logistic Growth Experiments

Experiments for the **GP-on-logistic-growth** design (see [GP_Logistic_Experimental_Design_Plan.md](../GP_Logistic_Experimental_Design_Plan.md) in the project root).

## Streamlined workflow (recommended)

Instead of running 20 individual experiment files, use the parameterized runner:

```matlab
% Run all 20 experiments and aggregate into logistic_metrics_all.csv/.mat
run_all_logistic_experiments();

% Or with options (no plots, subset, skip aggregate)
run_all_logistic_experiments(struct('make_plots', false, 'run_subset', [1 5 10], 'do_aggregate', false));
```

**Files:**
- **run_all_logistic_experiments.m** — Config-driven driver: runs all experiments and aggregates.
- **run_logistic_experiment.m** — Core runner: takes config `{exp_id, noise_pct, is_regular, is_auto}` and runs one experiment.
- **master_run_scripts.m** — Convenience: calls `run_all_logistic_experiments()`.
- **aggregate_logistic_metrics.m** — Unchanged; compiles all `results_exp_*.mat` into `logistic_metrics_all.csv` / `.mat`.

**Benefits:** One place to fix bugs, add kernels, or change N_list. Add a new condition by extending `build_logistic_configs()`.

---

## Legacy: Individual experiment scripts

The original `log_exp_NN_<scope>.m` files (01–20) still work and produce identical outputs. They are kept for reference; new work can use the streamlined runner.

## Naming convention

- **Runner scripts:** `log_exp_NN_<scope>.m` (legacy) or `run_logistic_experiment(cfg)`.
- **Results:** `results_exp_NN_<scope>.mat` (for aggregate) and `results_exp_NN_*_metrics.csv`.
- **Shared helpers:** `ground_truth_logistic.m`, `logistic_growth.m`, `metric_helpers.m` (project root).

## Experiments (all 20)

| Exp | Sampling | Noise | Optimization |
|-----|----------|-------|--------------|
| 01–05 | Regular | 0%, 1%, 5%, 10%, 20% | Default |
| 06–10 | Irregular | 0%, 1%, 5%, 10%, 20% | Default |
| 11–15 | Regular | 0%, 1%, 5%, 10%, 20% | Auto |
| 16–20 | Irregular | 0%, 1%, 5%, 10%, 20% | Auto |

Each varies N = 5, 10, 25, 50 and kernels: SqExp, Matern 1/2, 3/2, 5/2, Rational Quadratic.

## How to run

From the **project root** or this folder:

```matlab
cd('gp_logistic_design')
run_all_logistic_experiments
```

Or run one experiment only:

```matlab
cfg = struct('exp_id', 1, 'noise_pct', 0, 'is_regular', true, 'is_auto', false);
[T_point, T_prob, T_calib] = run_logistic_experiment(cfg);
```

## Shared code (in this folder)

- **logistic_growth.m** — Logistic ODE model (same as used by all_toy_problems).
- **ground_truth_logistic.m** — Returns high-res logistic curve (default 500 pts on [0, 50]).
- **metric_helpers.m** — Lives in the **project root** (`GP_and_SINDy/metric_helpers.m`) and is shared with `gp_LV_design`.
