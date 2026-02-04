# GP Reaction Kinetics Experiments

GP experiments on the **Robertson stiff reaction kinetics** system (3 states: X, Y, Z), mirroring the structure of `gp_logistic_design` and `gp_LV_design`. Parameters match **reaction_kinetics.m** and **all_toy_problems.m**: default `x0=1`, `y0=0`, `z0=0`, `tspan=[0 50]`.

## Shared code (in this folder)

- **reaction_kinetics.m** — Robertson stiff ODE model (same as used by all_toy_problems). Moved here from project root.
- **ground_truth_kinetics.m** — High-res kinetics curves (default 500 pts on [0, 50]) for X, Y, Z. Calls `reaction_kinetics(params, false)` and interpolates.
- **metric_helpers.m** — Lives in the **project root** (`GP_and_SINDy/metric_helpers.m`). Runners use `addpath('..')` to access it.

## Naming convention

- **Runner scripts:** `kin_exp_NN_<scope>.m` (e.g. `kin_exp_01_regular_0noise.m`).

## Experiments

| Script | Scope |
|--------|--------|
| **kin_exp_01_regular_0noise.m** | Ground truth kinetics; regular sampling; 0% noise; N = 5, 10, 25, 50; 5 kernels; one GP per state (X, Y, Z). |
| **kin_exp_02_regular_1noise.m** | Same as 01 with 1% Gaussian noise (sigma = 0.01×std per state); rng(42) for reproducibility. |
| **kin_exp_03_regular_5noise.m** | Same as 02 with 5% Gaussian noise (sigma = 0.05×std per state). |
| **kin_exp_04_regular_10noise.m** | Same as 03 with 10% Gaussian noise (sigma = 0.10×std per state). |
| **kin_exp_05_regular_20noise.m** | Same as 04 with 20% Gaussian noise (sigma = 0.20×std per state). |

Tables include a **State** column (X / Y / Z). Figures: RMSE vs N (3 subplots for X, Y, Z); one figure per kernel with 3×4 subplots (rows: X, Y, Z; cols: N = 5, 10, 25, 50).

## How to run

From the **project root** `GP_and_SINDy`:

```matlab
run('gp_kinetics_design/kin_exp_01_regular_0noise.m')
```

Or:

```matlab
cd('gp_kinetics_design')
kin_exp_01_regular_0noise
```

Ensure the project root is on the path (or use the script’s `addpath('..')`) so `metric_helpers` is found.
