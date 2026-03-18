# Bio-Informed GP Encoding: Design Plan

Design for encoding biological information into Gaussian processes, with a flexible framework for comparing encoding variants and a clear baseline.

---

## 1. Baseline: “Vanilla” GP

**The vanilla GP is your comparison anchor.** Every encoding variant should be compared against it on the same data.

| Component | Vanilla (baseline) | Purpose |
|-----------|--------------------|---------|
| **Mean function** | Constant (`'constant'` in fitrgp) | No prior shape; let data speak |
| **Kernel** | Squared Exponential (or Matérn 3/2) | Generic smoothness only |
| **Pseudo-observations** | None | No anchor points |
| **Constraints** | None | No derivative / monotonicity / bounds |
| **Likelihood** | Homoscedastic Gaussian | Standard observation model |

**Implementation:** `fitrgp(...)` with SE kernel, constant mean. Acute transient implemented in `gp_bio_informed_design/`.

**Metrics (from `metric_helpers`):** RMSE, MAE, R², NLPD, MSLL, CRPS, sMSE, Coverage, NLML.

---

## 2. Suggested File / Folder Structure

Use a **dedicated design folder** to stay consistent with your other designs:

```
GP_and_SINDy/
├── bio_informed_gp_priors.m          % Toy problems + ground truth (already exists)
├── metric_helpers.m                 % Shared metrics
├── gp_bio_informed_design/          % NEW folder
│   ├── README.md                    % Overview, how to run, encoding registry
│   ├── run_bio_experiment.m         % Main runner: toy_problem × encoding_config × data_config
│   ├── encoding_registry.m          % Registry of encoding configs (human-readable IDs)
│   ├── encodings/                   % Modular encoding components (optional subfolder)
│   │   ├── mean_functions.m         % Mean function builders (baseline, rise-and-fall, etc.)
│   │   ├── kernels.m               % Kernel builders (periodic, quasi-periodic, changepoint)
│   │   ├── pseudo_observations.m    % Anchor / pseudo-obs builders
│   │   └── constraints.m           % Placeholder for monotonicity, derivative signs (if using GPML/custom)
│   ├── ground_truth_bio.m           % Wrapper for bio_informed_gp_priors (consistent interface)
│   └── results/                    % Output .mat, .csv (or use out_dir from opts)
```

**Alternative (flatter):** If you prefer fewer files, fold `encodings/` into a single `gp_encoding_config.m` that returns a struct describing how to build the GP (mean, kernel, pseudo-obs, etc.).

---

## 3. Encoding Config: Modular and Extensible

Each **encoding variant** is a struct (or function handle) that describes *how to encode*:

```matlab
% Example: encoding config struct
enc = struct();
enc.id = 'circadian_periodic';           % Short, unique ID for tracking
enc.label = 'Circadian: periodic kernel'; % Human-readable label
enc.mean_fun = 'constant';               % or handle to custom mean
enc.kernel = 'periodic';                 % or handle / name for custom kernel
enc.pseudo_obs = [];                     % [] = none; or struct with t_pseudo, y_pseudo, sigma_pseudo
enc.constraints = [];                    % [] = none; or struct for monotonicity, etc.
```

**Registry idea:** `encoding_registry.m` returns a cell array or struct array of such configs. When you add a new variant, you add one entry. You never edit run logic—only the registry.

```matlab
% encoding_registry.m
function configs = encoding_registry()
    configs = {
        vanilla_config(),
        circadian_periodic_config(),
        circadian_quasi_periodic_config(),
        acute_transient_mean_only_config(),
        acute_transient_mean_plus_anchors_config(),
        delayed_changepoint_config(),
        ...
    };
end
```

---

## 4. Experiment Layout: What You Track

For each **toy problem** (acute_transient, bounded_activation, delayed_activation, circadian) you run:

- **Data config:** N points, regular/irregular, noise %, replicates
- **Encoding variants:** vanilla + your encoding configs

**Comparison matrix (conceptually):**

| Toy problem      | Data config | Vanilla | Encoding A | Encoding B | ... |
|------------------|-------------|---------|------------|------------|-----|
| acute_transient  | N=10, 5%    | RMSE=X  | RMSE=Y     | RMSE=Z     | ... |
| acute_transient  | N=25, 5%    | ...     | ...        | ...        | ... |
| circadian        | ...         | ...     | ...        | ...        | ... |

**Output naming** so you can trace exactly what was run:

```
results_bio_{toy}_{data_id}_{encoding_id}.mat
results_bio_acute_transient_N10_5noise_vanilla.mat
results_bio_acute_transient_N10_5noise_risefall_mean_only.mat
results_bio_acute_transient_N10_5noise_risefall_mean_plus_anchors.mat
results_bio_circadian_N25_1noise_vanilla.mat
results_bio_circadian_N25_1noise_periodic.mat
results_bio_circadian_N25_1noise_quasi_periodic.mat
```

---

## 5. Flexibility for Iterations

### 5.1 Adding a new encoding variant

1. Implement the encoding logic (custom mean, kernel, or pseudo-obs) in the appropriate module.
2. Add a config function (e.g. `acute_monotonicity_plus_anchors_config()`).
3. Register it in `encoding_registry.m`.
4. Re-run experiments—vanilla is always included automatically.

### 5.2 Comparing variants on the same problem

Your runner should:

1. Load or generate ground truth for the chosen toy problem.
2. Generate training data (t_obs, y_obs) from the data config.
3. For each encoding in the registry: fit GP, predict on t_gt, compute metrics.
4. Save results with the naming above.
5. Optionally produce a comparison table (e.g. one row per encoding, columns = metrics).

### 5.3 A/B-style comparisons

To compare e.g. “monotonicity only” vs “monotonicity + anchors”:

- **Encoding A:** `monotonicity_only` (constraints only, no pseudo-obs).
- **Encoding B:** `monotonicity_plus_anchors` (same constraints + pseudo-obs at t=0, t=20).

Both live in the registry. The runner treats them identically; only the config differs.

---

## 6. Implementation Choices and Tradeoffs

### 6.1 fitrgp vs custom GP (GPML, etc.)

- **fitrgp:** Easy, built-in, but limited: custom mean functions and kernels require `KernelFunction` handle; no built-in monotonicity or constraints.
- **GPML / custom:** Full flexibility (constraints, heteroscedastic, etc.) but more setup.

**Recommendation:** Start with fitrgp for:
- Vanilla baseline
- Custom mean functions (via `BasisFunction` handle or `Beta` + custom basis)
- Custom kernels (via `KernelFunction` handle—you already have `periodicKernel.m` in gp_LV_design)
- Pseudo-observations (by augmenting t_obs, y_obs with virtual points and small sigma)

Add GPML or custom code later when you need monotonicity or other hard constraints.

### 6.2 Pseudo-observations

Treat as extra (t, y) with a chosen observation noise. Example for “late time ≈ baseline”:

```matlab
t_extra = [0; 20];   % anchor at start and end
y_extra = [b; b];   % baseline
sigma_extra = 0.01; % strong prior
t_aug = [t_obs; t_extra];
y_aug = [y_obs; y_extra];
% Use fitrgp with observation weights or Sigma per-point (if supported)
% Or: add small sigma to these points in a custom implementation
```

fitrgp uses a single `Sigma`; for per-point noise you’d need a custom fit or to approximate (e.g. duplicate points with scaled weights).

### 6.3 Experiment registry as single source of truth

Keep `encoding_registry.m` as the **single place** that defines what encodings exist. The runner:

- Reads the registry
- Iterates over configs
- Uses the same metrics for all

New users (or future-you) can see all variants and their IDs in one file.

---

## 7. Suggested Implementation Order

1. **Create `gp_bio_informed_design/`** with `run_bio_experiment.m` that:
   - Takes (toy_problem, data_config, encoding_id or list)
   - Uses `bio_informed_gp_priors` for ground truth and data generation
   - Runs **vanilla only** first, computes metrics, saves
2. **Add `encoding_registry.m`** with vanilla + 1–2 encodings (e.g. periodic for circadian).
3. **Implement those encodings** (custom kernel handle, or augmented data for pseudo-obs).
4. **Run full comparison** for one toy problem × a few data configs.
5. **Add more encoding variants** incrementally; keep registry updated.

---

## 8. Summary

| Question | Answer |
|----------|--------|
| New file? | Yes—new folder `gp_bio_informed_design/` with modular structure |
| Baseline? | Vanilla GP: constant mean, SqExp kernel, no pseudo-obs, no constraints |
| Compare encodings? | Registry of configs; runner evaluates each on same data, saves with `{toy}_{data}_{encoding_id}` |
| Track what you tried? | Encoding registry + consistent output filenames |
| New user friendly? | README + registry as central index of all variants |
