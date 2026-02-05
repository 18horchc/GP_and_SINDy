# Exploring aggregated GP metrics

**explore_gp_metrics.m** (in the project root) loads the combined metrics from each design folder, applies filters, and plots one figure. Use it after you have run the aggregate scripts so that `*_metrics_all.mat` (or `.csv`) exist in `gp_logistic_design/`, `gp_LV_design/`, `gp_kinetics_design/`, and `gp_SINDy_design/`.

## How to use

1. Open **explore_gp_metrics.m**.
2. Set the **CONFIG** block at the top:
   - **MODEL** — which system: `'logistic'`, `'LV'`, `'kinetics'`, `'SINDy'`
   - **METRIC** — what to plot on y-axis: `'RMSE'`, `'MAE'`, `'R2'`, `'NLPD'`, `'MSLL'`, `'CRPS'`, `'sMSE'`, `'Coverage'`, `'NLML'`
   - **X_VAR** — what on x-axis: usually `'Noise'` or `'N'`
   - **LINE_VAR** — what defines each line: usually `'Kernel'` (or `'State'` for multi-state)
   - **FILTER_*** — restrict to certain N, Regular, Optimization, N_per_Time, State, Noise (use `[]` for “all”)
   - **SUBPLOT_BY_STATE** — for LV/kinetics/SINDy: `true` = one subplot per state, `false` = single plot
3. Run the script. One figure is produced.

## Example configurations

Copy the CONFIG values below into the CONFIG block in **explore_gp_metrics.m**, then run.

### Example 1: RMSE vs Noise, logistic, N=50, one line per kernel

- **Goal:** See how RMSE changes with noise level for regularly sampled data with 50 points; compare kernels.
- **CONFIG:**  
  `MODEL = 'logistic';`  
  `METRIC = 'RMSE';`  
  `X_VAR = 'Noise';`  
  `LINE_VAR = 'Kernel';`  
  `FILTER_N = 50;`  
  `FILTER_REGULAR = 'Yes';`  
  `FILTER_OPT = 'Default';`  
  `FILTER_NPT = 1;`  
  `FILTER_STATE = [];`  
  `FILTER_NOISE = [];`  
  `SUBPLOT_BY_STATE = true;`
- **Result:** One figure: y = RMSE, x = Noise, five lines (one per kernel).

---

### Example 2: RMSE vs N, logistic, 0% noise only

- **Goal:** Compare kernels as the number of points (N) changes, with no noise.
- **CONFIG:**  
  `MODEL = 'logistic';`  
  `METRIC = 'RMSE';`  
  `X_VAR = 'N';`  
  `LINE_VAR = 'Kernel';`  
  `FILTER_N = [];`  
  `FILTER_NOISE = 0;`  
  (other filters as in Example 1.)
- **Result:** One figure: RMSE vs N at 0% noise, one line per kernel.

---

### Example 3: Coverage vs Noise, SINDy, N=50, one subplot per state (M1, M2)

- **Goal:** Check 95% coverage for M1 and M2 as noise increases.
- **CONFIG:**  
  `MODEL = 'SINDy';`  
  `METRIC = 'Coverage';`  
  `X_VAR = 'Noise';`  
  `LINE_VAR = 'Kernel';`  
  `FILTER_N = 50;`  
  `FILTER_NOISE = [];`  
  `SUBPLOT_BY_STATE = true;`
- **Result:** One figure with two subplots (M1, M2); each shows Coverage vs Noise, one line per kernel.

---

### Example 4: NLPD vs N, kinetics, 5% noise, one subplot per state (X, Y, Z)

- **Goal:** Compare kernels on NLPD as N changes, at 5% noise, for each chemical.
- **CONFIG:**  
  `MODEL = 'kinetics';`  
  `METRIC = 'NLPD';`  
  `X_VAR = 'N';`  
  `LINE_VAR = 'Kernel';`  
  `FILTER_N = [];`  
  `FILTER_NOISE = 0.05;`  
  `SUBPLOT_BY_STATE = true;`
- **Result:** One figure with three subplots (X, Y, Z); each has NLPD vs N at 5% noise, one line per kernel.

---

### Example 5: R² vs Noise, LV, N=25, Prey only

- **Goal:** Single plot for Prey only: how R² varies with noise.
- **CONFIG:**  
  `MODEL = 'LV';`  
  `METRIC = 'R2';`  
  `X_VAR = 'Noise';`  
  `LINE_VAR = 'Kernel';`  
  `FILTER_N = 25;`  
  `FILTER_STATE = {'Prey'};`  
  `FILTER_NOISE = [];`  
  `SUBPLOT_BY_STATE = false;`
- **Result:** One plot: R² vs Noise for Prey, one line per kernel.

---

### Example 6: MAE vs N, logistic, 10% noise

- **Goal:** Compare kernels by MAE as N changes, at 10% noise.
- **CONFIG:**  
  `MODEL = 'logistic';`  
  `METRIC = 'MAE';`  
  `X_VAR = 'N';`  
  `LINE_VAR = 'Kernel';`  
  `FILTER_N = [];`  
  `FILTER_NOISE = 0.10;`
- **Result:** One figure: MAE vs N at 10% noise, one line per kernel.

---

More example blocks are in the comments at the bottom of **explore_gp_metrics.m**.
