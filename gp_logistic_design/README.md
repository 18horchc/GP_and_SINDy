# GP Logistic Growth Experiments

Experiments for the **GP-on-logistic-growth** design (see [GP_Logistic_Experimental_Design_Plan.md](../GP_Logistic_Experimental_Design_Plan.md) in the project root).

## Naming convention

- **Runner scripts:** `run_exp_NN_<scope>.m`  
  - `NN` = two-digit experiment number (01, 02, 03, …).  
  - `<scope>` = short description of what varies in this run (e.g. `regular_0noise`, `irregular_0noise`, `regular_noise`).

- **Results:** `results_exp_NN_<scope>.mat` / `.csv` (when saved).

- **Shared helpers:** Descriptive names, e.g. `ground_truth_logistic.m`, `sample_observations.m` (when added).

This keeps experiments ordered, easy to extend, and clear in the big picture.

## Experiments

| Script | Scope |
|--------|--------|
| **run_exp_01_regular_0noise.m** | Ground truth; regular sampling; 0% noise; N = 5, 10, 25, 50; all 5 kernels. |

*(Add new rows as you add steps, e.g. 02 = irregular 0% noise, 03 = regular + noise, etc.)*

## How to run

From the **project root** `GP_and_SINDy`:

```matlab
run('gp_logistic_design/run_exp_01_regular_0noise.m')
```

Or from this folder (parent is added to path automatically):

```matlab
cd('gp_logistic_design')
run_exp_01_regular_0noise
```

## Shared code

- **ground_truth_logistic.m** — Returns high-res logistic curve (default 500 pts on [0, 50]). Used by all experiments.
- **metric_helpers.m** — Computes the three metric groups: point-estimate (RMSE, MAE, R²), probabilistic (NLPD, MSLL, CRPS), and calibration (sMSE, Coverage, NLML). Used by experiment runners.
