# Implementing Proper Derivative Observations in GPs

This guide walks through implementing a GP with derivative observations (observations of df/dt) rather than the virtual "shape guide" approximation.

---

## 1. Theory

A GP prior on f(x) implies that f and its derivatives are jointly Gaussian. For kernel k(x, x'):

| Block | Covariance |
|-------|------------|
| f and f | K_ff(i,j) = k(x_i, x_j) |
| f and df/dx' | K_fdf(i,j) = ∂k(x_i, x_j')/∂x_j' |
| df/dx and f | K_dff(i,j) = ∂k(x_i', x_j)/∂x_i' |
| df/dx and df/dx' | K_dfdf(i,j) = ∂²k(x_i', x_j')/∂x_i'∂x_j' |

**Observation model:**
- Function obs: y_i = f(x_i) + ε_i  (Gaussian noise σ²)
- Derivative obs: d_j = df(x_j)/dx + ε_j

**Joint vector:** z = [y_1, ..., y_n, d_1, ..., d_m]'

**Prediction:** Standard GP equations with the block-structured covariance.

---

## 2. SE Kernel Derivatives

For k(x, x') = σ_f² exp(-(x-x')²/(2ℓ²)) with r = x - x':

- **k(x, x')** = σ_f² exp(-r²/(2ℓ²))
- **∂k/∂x'** = k · (x - x')/ℓ²
- **∂k/∂x** = k · (x' - x)/ℓ²  = -∂k/∂x'
- **∂²k/∂x∂x'** = k · [1/ℓ² - (x-x')²/ℓ⁴]

---

## 3. Implementation Steps

### Step 1: Build the full covariance matrix K

```
Inputs:
  X_f   = [x_1; ...; x_n]     (function observation locations)
  X_df  = [x_{n+1}; ...; x_{n+m}]  (derivative observation locations)
  hyp   = [log(ℓ); log(σ_f); log(σ_n)]

Output: K_total ( (n+m) × (n+m) )
  K_total = [K_ff,   K_fdf;
             K_dff,  K_dfdf]
  + σ_n² * I
```

### Step 2: Build the mean vector

- Function obs: mean = m(x_i) for each x_i (e.g. constant baseline)
- Derivative obs: mean = dm/dx(x_j). For constant mean m(x)=c, dm/dx = 0.

### Step 3: Assemble observation vector

```
y_obs = [y_1; ...; y_n]     (observed function values)
d_obs = [d_1; ...; d_m]     (observed derivative values, e.g. +0.5 at t=2, -0.5 at t=10)

z = [y_obs; d_obs]
```

### Step 4: Solve for α

```
α = K_total \ (z - mean_vec)
```

### Step 5: Predictive mean at test point x*

```
k_star = [k(x*, X_f); ∂k(x*, X_df)/∂x*]   (for mean of f)
         or for derivatives: [∂k/∂x*; ∂²k/∂x*∂X_df]

μ* = mean(x*) + k_star' * α
```

### Step 6: Predictive variance

```
v* = k(x*, x*) - k_star' * (K_total \ k_star)
```

---

## 4. What You Need to Specify

1. **Derivative observation locations** (e.g. t = 2 and t = 10)
2. **Derivative values** (e.g. "df/dt > 0 at t=2" → use d = +0.5 as soft constraint, or treat as inequality)
3. **Encoding "slope decreases"**: 
   - Option A: Observe df/dt at t=2 and t=10 with values d_2 > d_10 (e.g. d_2 = 0.3, d_10 = -0.5)
   - Option B: Use inequality constraints (harder; requires constrained optimization or EP)

For a **soft** encoding, use Option A with chosen values. The GP will pull the posterior toward having slope ≈ d_2 at t=2 and ≈ d_10 at t=10.

---

## 5. Implementation Options

### A. Standalone MATLAB (no GPML)

Write `fit_gp_deriv.m` that:
- Takes (x_f, y_f, x_df, d_obs, x_test)
- Builds K using the formulas above
- Solves and predicts
- Uses `fminunc` to optimize hyp (log ℓ, log σ_f, log σ_n) via NLML

### B. GPML extension

GPML does not natively support derivative observations. You would need to:
- Create a new covariance that accepts "typed" inputs (function vs derivative)
- Or wrap the standalone implementation and use it instead of gp()

### C. Python (GPyTorch, GPflow)

Some Python GP libraries have derivative support. For MATLAB-only workflow, A is most practical.

---

## 6. Recommended Next Step

Implement the **standalone** `fit_gp_deriv.m` in this project. It will:
1. Build K_ff, K_fdf, K_dff, K_dfdf for SE kernel
2. Optimize hyp via NLML
3. Predict mean and variance at test points
4. Integrate with `run_naive_vs_deriv_encoding.m` as a third model (true derivative obs)

Would you like me to implement this standalone function?
