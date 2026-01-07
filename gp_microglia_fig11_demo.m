function gp_microglia_fig11_demo()
% GP prior/posterior demo + microglia posterior plots
% Recreates Figure 1.1(a) style prior samples, then a Figure 1.1(b)-style
% posterior conditioned on microglia M1/M2 data (with noise -> "close", not exact).

clear; close all; clc;
rng(7); % reproducibility

%% --------------------- MICROGLIA DATA ---------------------
tpoints_M1 = [0,0,1,1,2,2,3,3,3,3,3,5,7,7,7,7,14,14,14]';
datapointsM1 = [0,10,5,50,120,125,375,62,102,100,60,325,600,55,900,225,750,400,1300]';

tpoints_M2 = [0,0, 1, 1,1, 2, 2, 3,3,3,3,3, 5, 7, 7, 7, 7, 14,14,14]';
datapointsM2 = [0,10,170,15,50, 90,269,300,15,57,100,160, 800, 600,6,400,270,200,110,100]';

%% --------------------- KERNEL (SE/RBF) ---------------------
% k(x,x') = sf2 * exp( -0.5*(x-x')^2 / ell^2 )
kernel = @(xa, xb, ell, sf2) se_kernel(xa, xb, ell, sf2);

%% =====================================================================
%  FIGURE 1.1(a): PRIOR SAMPLES (generic x in [0,1])
% =====================================================================
x = linspace(0,1,250)';  % "input, x" like the textbook figure
m0 = zeros(size(x));     % zero mean prior

ell_prior = 0.20;  % lengthscale: smaller = wigglier, larger = smoother
sf2_prior = 1.0;   % prior variance (constant across x for SE kernel)

Kxx = kernel(x, x, ell_prior, sf2_prior);
Kxx = jitter(Kxx);

% Draw 4 samples from the GP prior
nsamp = 4;
Fprior = gp_sample(m0, Kxx, nsamp);

% Prior pointwise std (diagonal)
prior_std = sqrt(diag(Kxx));

figure('Color','w','Name','Figure 1.1(a) Prior');
hold on;

% Shaded ±2 std region
fill_band_lohi(x, m0 - 2*prior_std, m0 + 2*prior_std);

% Plot samples (solid)
plot(x, Fprior, 'LineWidth', 1.2);

xlabel('input, x');
ylabel('f(x)');
title('(a) prior');
box on; grid off;

%% =====================================================================
%  FIGURE 1.1(b) STYLE: POSTERIOR FOR MICROGLIA (M1 and M2)
% =====================================================================
tgrid = linspace(0,14,400)';   % time grid for smooth curves

% --- Hyperparameters ---
sf2_M1 = var(datapointsM1);
sf2_M2 = var(datapointsM2);

% --- Positivity via log transform ---
c1 = 1; z1 = log(datapointsM1 + c1);
c2 = 1; z2 = log(datapointsM2 + c2);

% Noise setup
sigma_n_M1 = 0.15 * std(z1);
sigma_n_M2 = 0.15 * std(z2);

% === NEW: AUTO-OPTIMIZE ELL ===
fprintf('Optimizing M1...\n');
ell_M1 = optimize_ell_1D(tpoints_M1, z1, var(z1), sigma_n_M1, kernel);

fprintf('Optimizing M2...\n');
ell_M2 = optimize_ell_1D(tpoints_M2, z2, var(z2), sigma_n_M2, kernel);
% ==============================

% Now use ell_M1 and ell_M2 in your posterior calls:
[mu1z, s21z, Fpost1z] = gp_posterior_and_samples( ...
    tpoints_M1, z1, tgrid, ell_M1, var(z1), sigma_n_M1, nsamp, kernel);

[mu2z, s22z, Fpost2z] = gp_posterior_and_samples( ...
    tpoints_M2, z2, tgrid, ell_M2, var(z2), sigma_n_M2, nsamp, kernel);

% transform mean + samples back to strictly-positive scale (y + c)
mu1 = exp(mu1z);
mu2 = exp(mu2z);

Fpost1 = exp(Fpost1z);
Fpost2 = exp(Fpost2z);

% 95% pointwise band on strictly-positive scale (y + c)
lo1 = exp(mu1z - 1.96*sqrt(s21z));
hi1 = exp(mu1z + 1.96*sqrt(s21z));

lo2 = exp(mu2z - 1.96*sqrt(s22z));
hi2 = exp(mu2z + 1.96*sqrt(s22z));


% Shift posterior back to original y-scale for visualization only
mu1_plot   = mu1   - c1;
lo1_plot   = lo1   - c1;
hi1_plot   = hi1   - c1;
Fpost1_plot = Fpost1 - c1;

mu2_plot   = mu2   - c2;
lo2_plot   = lo2   - c2;
hi2_plot   = hi2   - c2;
Fpost2_plot = Fpost2 - c2;



figure('Color','w','Name','Microglia GP Posterior (Figure 1.1(b) style)');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% --------- M1 subplot ----------
nexttile; hold on;
fill_band_lohi(tgrid, lo1_plot, hi1_plot);
plot(tgrid, mu1_plot, 'k-', 'LineWidth', 1.8);
plot(tgrid, Fpost1_plot, 'k--', 'LineWidth', 1.0);

% Plot ORIGINAL data (unchanged)
plot(tpoints_M1, datapointsM1, 'ko', 'MarkerFaceColor','w', 'LineWidth',1.2);

xlabel('time (days)');
ylabel('M1 count');
title('(b) posterior: M1');
box on;


% --------- M2 subplot ----------
nexttile; hold on;
fill_band_lohi(tgrid, lo2_plot, hi2_plot);
plot(tgrid, mu2_plot, 'k-', 'LineWidth', 1.8);
plot(tgrid, Fpost2_plot, 'k--', 'LineWidth', 1.0);

plot(tpoints_M2, datapointsM2, 'ko', 'MarkerFaceColor','w', 'LineWidth',1.2);

xlabel('time (days)');
ylabel('M2 count');
title('(b) posterior: M2');
box on;


%% Notes for you:
% - To force curves to pass *very* close to points: reduce sigma_n_M1/M2.
% - To make curves smoother: increase ell_data.
% - If extrapolation "blows up": usually increase noise a bit and/or increase ell.

end

%% ===================== Helper functions =====================

function K = se_kernel(xa, xb, ell, sf2)
% Squared exponential (RBF) kernel
xa = xa(:); xb = xb(:)';
d2 = (xa - xb).^2;
K = sf2 * exp(-0.5 * d2 / (ell^2));
end

function K = jitter(K)
% Small diagonal jitter for numerical stability
K = K + 1e-9 * eye(size(K));
end

function F = gp_sample(m, K, nsamp)
% Sample from N(m, K) using Cholesky
L = chol(K, 'lower');
Z = randn(size(K,1), nsamp);
F = m + L * Z;
end

function fill_band_lohi(x, lo, hi)
% Shade region between lo(x) and hi(x)
x  = x(:);
lo = lo(:);
hi = hi(:);

X = [x; flipud(x)];
Y = [hi; flipud(lo)];

h = fill(X, Y, [0.85 0.85 0.85], 'EdgeColor','none');
uistack(h,'bottom');
end


function [mu_s, var_s, Fsamp] = gp_posterior_and_samples( ...
    xtrain, ytrain, xtest, ell, sf2, sigma_n, nsamp, kernel)

xtrain = xtrain(:); ytrain = ytrain(:); xtest = xtest(:);

% Covariances
Ktt = kernel(xtrain, xtrain, ell, sf2);
Kts = kernel(xtrain, xtest,  ell, sf2);
Kss = kernel(xtest,  xtest,  ell, sf2);

% Add noise on training observations (preference to be close, not exact)
Kyy = Ktt + (sigma_n^2) * eye(length(xtrain));
Kyy = jitter(Kyy);

% Posterior mean: mu = K_{*t} (K_{tt}+sn^2 I)^{-1} y
L = chol(Kyy, 'lower');
alpha = L'\(L\ytrain);
mu_s = Kts' * alpha;

% Posterior covariance: Sigma = Kss - K_{*t} (Kyy)^{-1} K_{t*}
v = L \ Kts;
Sigma_s = Kss - v' * v;
Sigma_s = jitter((Sigma_s + Sigma_s')/2); % symmetrize + jitter

var_s = max(diag(Sigma_s), 0);

% Draw posterior samples
m = mu_s;
Fsamp = gp_sample(m, Sigma_s, nsamp);
end

function best_ell = optimize_ell_1D(xtrain, ytrain, current_sf2, current_sigma_n, kernel_func)
% OPTIMIZE_ELL_1D Finds the best lengthscale by minimizing Negative Log Likelihood
%
% Inputs:
%   xtrain, ytrain: The data (in the transformed log-space z)
%   current_sf2:    Signal variance (kept fixed)
%   current_sigma_n: Noise standard deviation (kept fixed)
%   kernel_func:    Function handle for the kernel (e.g. @se_kernel)

    % Define the objective function (NLL) as a function of ell
    % We search in the range [0.1, 10] days. Adjust bounds if needed.
    obj_fun = @(ell_guess) gp_neg_log_likelihood(ell_guess, xtrain, ytrain, ...
                                                 current_sf2, current_sigma_n, kernel_func);

    % Use fminbnd for simple 1D optimization
    options = optimset('Display','off');
    best_ell = fminbnd(obj_fun, 0.1, 10.0, options);

    fprintf('Optimized lengthscale: %.4f\n', best_ell);
end

function NLL = gp_neg_log_likelihood(ell, x, y, sf2, sigma_n, kernel)
    x = x(:); y = y(:);
    n = length(y);
    
    % 1. Construct Covariance Matrix K
    K = kernel(x, x, ell, sf2);
    
    % 2. Add Noise (and small jitter for stability)
    Kyy = K + (sigma_n^2 * eye(n)) + (1e-9 * eye(n));
    
    % 3. Cholesky Decomposition (L * L' = Kyy)
    [L, p] = chol(Kyy, 'lower');
    if p > 0
        NLL = inf; % Penalty if matrix is not positive definite
        return;
    end
    
    % 4. Compute NLL terms
    % Term 1: Data fit = 0.5 * y' * K^-1 * y
    alpha = L' \ (L \ y);
    data_fit = 0.5 * y' * alpha;
    
    % Term 2: Complexity penalty = 0.5 * log|K| = sum(log(diag(L)))
    complexity = sum(log(diag(L)));
    
    % Term 3: Normalization constant (0.5 * n * log(2pi))
    % (constant, so often ignored in optimization, but strictly part of NLL)
    normalization = 0.5 * n * log(2*pi);
    
    NLL = data_fit + complexity + normalization;
end
