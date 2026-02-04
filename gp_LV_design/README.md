# GP Lotka-Volterra Experiments

GP experiments on the **Lotka-Volterra** predator–prey model, mirroring the structure of `gp_logistic_design` (same metrics, same kernels, same N and noise ideas). LV parameters match **lotka_volterra_model.m** and **all_toy_problems.m**: default `alpha=1`, `beta=0.2`, `delta=0.5`, `gamma=0.2`, `x0=1`, `y0=2`, `tspan=[0 50]`.

## Shared code

- **ground_truth_LV.m** — High-res LV curves (default 500 pts on [0, 50]) for prey and predator. Calls `lotka_volterra_model(params, false)` and interpolates.
- **metric_helpers.m** — Lives in the **project root** (`GP_and_SINDy/metric_helpers.m`). Same as for logistic; runners use `addpath('..')` so both `gp_logistic_design` and `gp_LV_design` use it.

## Experiments

| Script | Scope |
|--------|--------|
| **run_exp_01_regular_0noise.m** | Ground truth LV; regular sampling; 0% noise; N = 5, 10, 25, 50; all 5 kernels; one GP per state (Prey, Predator). |

Tables include a **State** column (Prey / Predator). Figures: one RMSE vs N (Prey and Predator side by side); one figure per kernel with 2×4 subplots (row 1: Prey for N=5,10,25,50; row 2: Predator).

## How to run

From the **project root** `GP_and_SINDy`:

```matlab
run('gp_LV_design/run_exp_01_regular_0noise.m')
```

Or:

```matlab
cd('gp_LV_design')
run_exp_01_regular_0noise
```

Ensure the project root is on the path (or use the script’s `addpath('..')`) so `metric_helpers` and `lotka_volterra_model` are found.
