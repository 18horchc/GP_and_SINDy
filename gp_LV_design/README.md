# GP Lotka-Volterra Experiments

GP experiments on the **Lotka-Volterra** predator–prey model, mirroring the structure of `gp_logistic_design`. LV parameters match **lotka_volterra_model.m** and **all_toy_problems.m**: default `alpha=1`, `beta=0.2`, `delta=0.5`, `gamma=0.2`, `x0=1`, `y0=2`, `tspan=[0 50]`.

## Streamlined workflow (recommended)

Use the parameterized runner (same interface as gp_logistic_design):

```matlab
% Run all 20 experiments and aggregate into LV_metrics_all.csv/.mat
run_all_LV_experiments();

% Or with options (no plots, subset, skip aggregate)
run_all_LV_experiments(struct('make_plots', false, 'run_subset', 1:10, 'do_aggregate', false));
```

**Files:**
- **run_all_LV_experiments.m** — Config-driven driver: runs all experiments and aggregates.
- **run_LV_experiment.m** — Core runner: takes config `{exp_id, noise_pct, is_regular, is_auto}` and runs one experiment.
- **master_run_scripts.m** — Convenience: calls `run_all_LV_experiments()`.
- **aggregate_LV_metrics.m** — Unchanged; compiles all `results_exp_*.mat` into `LV_metrics_all.csv` / `.mat`.

**LV-specific:** Two states (Prey, Predator); six kernels (5 built-in + custom `periodicKernel`). Periodic kernel uses `KernelParameters` (no auto optimization). Tables include a **State** column.

---

## Legacy: Individual experiment scripts

The original `LV_exp_NN_<scope>.m` files (01–20) are in the `old/` folder for reference.

## Shared code (in this folder)

- **lotka_volterra_model.m** — LV ODE model (same as used by all_toy_problems).
- **ground_truth_LV.m** — High-res LV curves (500 pts on [0, 50]) for prey and predator.
- **periodicKernel.m** — Custom periodic kernel for fitrgp (period = t_range/3, length scale = period/4, sigmaF = std(Y)).
- **metric_helpers.m** — Lives in the **project root**. Runners use `addpath('..')` to access it.

## Experiments (all 20)

| Exp | Sampling | Noise | Optimization |
|-----|----------|-------|--------------|
| 01–05 | Regular | 0%, 1%, 5%, 10%, 20% | Default |
| 06–10 | Irregular | 0%, 1%, 5%, 10%, 20% | Default |
| 11–15 | Regular | 0%, 1%, 5%, 10%, 20% | Auto |
| 16–20 | Irregular | 0%, 1%, 5%, 10%, 20% | Auto |

Each varies N = 5, 10, 25, 50 and kernels: SqExp, Matern 1/2, 3/2, 5/2, Rational Quadratic, **Periodic**. One GP per state (Prey, Predator).

## How to run

From the **project root** or this folder:

```matlab
cd('gp_LV_design')
run_all_LV_experiments
```

Or run one experiment only:

```matlab
cfg = struct('exp_id', 1, 'noise_pct', 0, 'is_regular', true, 'is_auto', false);
[T_point, T_prob, T_calib] = run_LV_experiment(cfg);
```
