# Monotonicity via Expectation Propagation (Riihimäki & Vehtari 2010)

Implementation of Gaussian processes with monotonicity information for the acute transient problem.

## Algorithm Overview

1. **Virtual derivative locations** `Tm` in the monotone interval (e.g. [2,10])
2. **Prior**: f(t) ~ GP(0, k_SE) with k_SE(t,t'; η, ℓ)
3. **Covariance blocks**: K_ff, K_fd, K_dd, K_joint
4. **Likelihood**: p(Y|f) = N(Y|f, σ²I) for observations
5. **Monotonicity**: p(m_j|d_j) = Φ(-d_j/ν) at each virtual point (negative derivative)
6. **EP**: Approximate posterior with Gaussian sites; iterate until convergence
7. **Hyperparameters**: Optimize θ = {η, ℓ, σ} via EP approx. to marginal likelihood
8. **Prediction**: Standard GP prediction with virtual derivative sites
9. **Diagnostic**: Check P(f'(t) > 0) on grid; add points if too high

## Files

- `fit_gp_monotonic_ep.m` – Core EP implementation
- `predict_deriv_monotonic_ep.m` – Derivative posterior at arbitrary points
- `run_monotonic_ep_acute_transient.m` – Full pipeline for acute transient

## Usage

```matlab
% Basic usage
opts = struct('Tm', [2; 3.5; 5; 6.5; 8; 9.5], 'nu', 1e-6);
[ymu, ys2, hyp, post_d_mean, post_d_var, mu_site, sigma2_site] = ...
    fit_gp_monotonic_ep(t_obs, y_obs, [2, 10], t_gt, opts);

% Derivative diagnostic
t_check = linspace(2, 10, 100)';
[d_mean, d_var] = predict_deriv_monotonic_ep(t_obs, y_obs, Tm, hyp, mu_site, sigma2_site, t_check);
prob_positive = 1 - normcdf(0, d_mean, sqrt(d_var));
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Tm | [2, 3.5, 5, 6.5, 8, 9.5] | Virtual derivative locations |
| nu | 1e-6 | Monotonicity softness (→0 = stricter) |
| hyp_init | auto | [log(η); log(ℓ); log(σ)] |
| ep_max | 100 | Max EP iterations |
| ep_tol | 1e-4 | EP convergence tolerance |

## Reference

Riihimäki, J. & Vehtari, A. (2010). Gaussian processes with monotonicity information. AISTATS.
