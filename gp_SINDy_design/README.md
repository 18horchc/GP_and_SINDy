# GP SINDy+MCMC Microglia Experiments

GP experiments on the **SINDy+MCMC microglial cell dynamics** model (Amato & Arnold, 2025). Two states: M1, M2. Same structure as `gp_logistic_design` and `gp_LV_design`: same kernels, N list, and metric tables. Parameters match **sindy_mcmc_microglia.m** and **all_toy_problems.m** (default MCMC posterior means, M0=[5;5], tspan=[0 50]).

## Shared code (in this folder)

- **sindy_mcmc_microglia.m** — SINDy+MCMC ODE model (same as used by all_toy_problems). Moved here from project root.
- **ground_truth_SINDy.m** — High-res M1, M2 curves (default 500 pts on [0, 50]). Calls `sindy_mcmc_microglia(params, false)` and interpolates.
- **metric_helpers.m** — Lives in the **project root** (`GP_and_SINDy/metric_helpers.m`). Runners use `addpath('..')` to access it.

## Naming convention

- **Runner scripts:** `SINDy_exp_NN_<scope>.m` (e.g. `SINDy_exp_01_regular_0noise.m`).

## Experiments

| Script | Scope |
|--------|--------|
| **SINDy_exp_01_regular_0noise.m** | Ground truth SINDy; regular sampling; 0% noise; N = 5, 10, 25, 50; 5 kernels (SqExp, Matern 1/2, 3/2, 5/2, RatQuad); one GP per state (M1, M2). |

Tables include a **State** column (M1 / M2). Figures: RMSE vs N (2 subplots for M1, M2); one figure per kernel with 2×4 subplots (row 1: M1 for N=5,10,25,50; row 2: M2).

## How to run

From the **project root** `GP_and_SINDy`:

```matlab
run('gp_SINDy_design/SINDy_exp_01_regular_0noise.m')
```

Or:

```matlab
cd('gp_SINDy_design')
SINDy_exp_01_regular_0noise
```

Ensure the project root is on the path so `metric_helpers` is found.
