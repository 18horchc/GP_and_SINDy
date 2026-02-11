# GP Reaction Kinetics Experiments

GP experiments on the **Robertson stiff reaction kinetics** system (3 states: X, Y, Z), mirroring the structure of `gp_logistic_design` and `gp_LV_design`. Parameters match **reaction_kinetics.m** and **all_toy_problems.m**: default `x0=1`, `y0=0`, `z0=0`, `tspan=[0 50]`.

## Streamlined workflow (recommended)

Use the parameterized runner (same interface as gp_logistic_design and gp_LV_design):

```matlab
% Run all 20 experiments and aggregate into kinetics_metrics_all.csv/.mat
run_all_kinetics_experiments();

% Or with options (no plots, subset, skip aggregate)
run_all_kinetics_experiments(struct('make_plots', false, 'run_subset', 1:10, 'do_aggregate', false));
```

**Files:**
- **run_all_kinetics_experiments.m** — Config-driven driver: runs all experiments and aggregates.
- **run_kinetics_experiment.m** — Core runner: takes config `{exp_id, noise_pct, is_regular, is_auto}` and runs one experiment.
- **master_run_scripts.m** — Convenience: calls `run_all_kinetics_experiments()`.
- **aggregate_kinetics_metrics.m** — Unchanged; compiles all `results_kin_exp_*.mat` into `kinetics_metrics_all.csv` / `.mat`.

**Kinetics-specific:** Three states (X, Y, Z). Five built-in kernels only (no periodic).

---

## Legacy: Individual experiment scripts

The original `kin_exp_NN_<scope>.m` files (01–20) are in the `old/` folder for reference.

## Shared code (in this folder)

- **reaction_kinetics.m** — Robertson stiff ODE model (same as used by all_toy_problems).
- **ground_truth_kinetics.m** — High-res kinetics curves (500 pts on [0, 50]) for X, Y, Z.
- **metric_helpers.m** — Lives in the **project root**. Runners use `addpath('..')` to access it.

## Experiments (all 20)

| Exp | Sampling | Noise | Optimization |
|-----|----------|-------|--------------|
| 01–05 | Regular | 0%, 1%, 5%, 10%, 20% | Default |
| 06–10 | Irregular | 0%, 1%, 5%, 10%, 20% | Default |
| 11–15 | Regular | 0%, 1%, 5%, 10%, 20% | Auto |
| 16–20 | Irregular | 0%, 1%, 5%, 10%, 20% | Auto |

Each varies N = 5, 10, 25, 50 and kernels: SqExp, Matern 1/2, 3/2, 5/2, Rational Quadratic. One GP per state (X, Y, Z).

## How to run

From the **project root** or this folder:

```matlab
cd('gp_kinetics_design')
run_all_kinetics_experiments
```

Or run one experiment only:

```matlab
cfg = struct('exp_id', 1, 'noise_pct', 0, 'is_regular', true, 'is_auto', false);
[T_point, T_prob, T_calib] = run_kinetics_experiment(cfg);
```
