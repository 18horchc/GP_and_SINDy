# GP Lotka-Volterra Experiments

GP experiments on the **Lotka-Volterra** predator–prey model, mirroring the structure of `gp_logistic_design` (same metrics, same kernels, same N and noise ideas). LV parameters match **lotka_volterra_model.m** and **all_toy_problems.m**: default `alpha=1`, `beta=0.2`, `delta=0.5`, `gamma=0.2`, `x0=1`, `y0=2`, `tspan=[0 50]`.

## Shared code (in this folder)

- **lotka_volterra_model.m** — LV ODE model (same as used by all_toy_problems). Moved here from project root.
- **ground_truth_LV.m** — High-res LV curves (default 500 pts on [0, 50]) for prey and predator. Calls `lotka_volterra_model(params, false)` and interpolates.
- **periodicKernel.m** — Custom periodic kernel for fitrgp; used by LV_exp_01 (initial period = t_range/3, length scale = period/4, sigmaF = std(Y)).
- **metric_helpers.m** — Lives in the **project root** (`GP_and_SINDy/metric_helpers.m`). Runners use `addpath('..')` to access it.

## Naming convention

- **Runner scripts:** `LV_exp_NN_<scope>.m` (e.g. `LV_exp_01_regular_0noise.m`). Logistic design uses `log_exp_NN_...` in its folder.

## Experiments

| Script | Scope |
|--------|--------|
| **LV_exp_01_regular_0noise.m** | Ground truth LV; regular sampling; 0% noise; N = 5, 10, 25, 50; all 6 kernels; one GP per state (Prey, Predator). |
| **LV_exp_02_regular_1noise.m** | Same as 01 with 1% Gaussian noise on sampled prey and predator (sigma = 0.01*std per state); rng(42) for reproducibility. |
| **LV_exp_03_regular_5noise.m** | Same as 02 with 5% Gaussian noise (sigma = 0.05*std per state). |
| **LV_exp_04_regular_10noise.m** | Same as 03 with 10% Gaussian noise (sigma = 0.10*std per state). |
| **LV_exp_05_regular_20noise.m** | Same as 04 with 20% Gaussian noise (sigma = 0.20*std per state). |

Tables include a **State** column (Prey / Predator). Figures: one RMSE vs N (Prey and Predator side by side); one figure per kernel with 2×4 subplots (row 1: Prey for N=5,10,25,50; row 2: Predator).

## How to run

From the **project root** `GP_and_SINDy`:

```matlab
run('gp_LV_design/LV_exp_01_regular_0noise.m')
```

Or:

```matlab
cd('gp_LV_design')
LV_exp_01_regular_0noise
```

Ensure the project root is on the path (or use the script’s `addpath('..')`) so `metric_helpers` and `lotka_volterra_model` are found.
