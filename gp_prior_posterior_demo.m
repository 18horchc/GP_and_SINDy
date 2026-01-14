%% Gaussian Process Prior and Posterior Demonstration
% This script replicates the GP prior/posterior plots from the textbook
% using equation 2.16 for the covariance function and equation 2.19 for
% posterior sampling, with conditioning formulas from appendix A.2
%
% Reference: Equation 2.16: k(x_p, x_q) = exp(-1/2 |x_p - x_q|^2)

clear; close all; clc;

%% Parameters
x_min = -5;
x_max = 5;
n_grid = 200;  % Number of points on regular grid
n_samples = 4;  % Number of random function draws
n_train = 5;    % Number of training data points

%% Create regular grid of scalar inputs
x_grid = linspace(x_min, x_max, n_grid)';

%% Define covariance function (Equation 2.16)
% k(x_p, x_q) = exp(-1/2 |x_p - x_q|^2)
% This is the RBF/Squared Exponential kernel with length scale = 1
% Vectorized version for efficiency
compute_kernel = @(X1, X2) exp(-0.5 * (X1 - X2').^2);

%to change length scale it will be (X1-X2')/l where l is the length scale

%additionally, the overall vaiance (magnitude) could be changed by 
%multiplying by a positive pre-factor before the exp

%% ============================================
%% PART (a): PRIOR
%% ============================================
fprintf('Generating GP prior samples...\n');

% Compute prior covariance matrix K_prior on the grid
% K_prior[i,j] = k(x_i, x_j) for all pairs in x_grid
K_prior = compute_kernel(x_grid, x_grid);

% Ensure numerical stability (add small jitter)
K_prior = K_prior + 1e-9 * eye(n_grid);

% Prior mean is zero
mu_prior = zeros(n_grid, 1);

% Draw samples from prior: f ~ N(0, K_prior)
% Using Cholesky decomposition for numerical stability
L_prior = chol(K_prior, 'lower');
samples_prior = zeros(n_grid, n_samples);
for i = 1:n_samples
    z = randn(n_grid, 1);
    samples_prior(:, i) = mu_prior + L_prior * z;
end

% Compute uncertainty bands (mean ± 2*std)
std_prior = sqrt(diag(K_prior));
mu_upper_prior = mu_prior + 2 * std_prior;
mu_lower_prior = mu_prior - 2 * std_prior;

%% ============================================
%% PART (b): POSTERIOR
%% ============================================
fprintf('Generating GP posterior samples...\n');

% Invent some training data points
%these must be column vectors
% Based on the image, approximate locations:
% (-4, -2), (-2.5, 0.5), (-1, 1.5), (0, 1.8), (2.5, 0)
%x_train = [-4; -2.5; -1; 0; 2.5];
%y_train = [-2; 0.5; 1.5; 1.8; 0];

%new data: 
x_train = [-3.5;-2;0;1.53;2.1];
y_train = [-2; 0; 1.1; 1.8; -0.2];

% For noise-free observations, we assume zero observation noise
% (as mentioned in the image: "five noise free observations")
sigma_n = 1e-6;  % Very small noise for numerical stability

% Compute covariance matrices needed for posterior
% K_train: covariance between training points
K_train = compute_kernel(x_train, x_train);
K_train = K_train + sigma_n^2 * eye(n_train);  % Add observation noise

% K_cross: covariance between training points and grid points
K_cross = compute_kernel(x_train, x_grid);

% K_grid: covariance between grid points (same as K_prior)
K_grid = K_prior;

% Posterior mean (Equation 2.19 / A.2 conditioning formula)
% μ* = K_cross^T * K_train^(-1) * y_train
K_train_inv = K_train \ eye(n_train);
mu_posterior = K_cross' * (K_train_inv * y_train);

% Posterior covariance (Equation 2.19 / A.2 conditioning formula)
% Σ* = K_grid - K_cross^T * K_train^(-1) * K_cross
K_posterior = K_grid - K_cross' * (K_train_inv * K_cross);

% Ensure numerical stability
K_posterior = (K_posterior + K_posterior') / 2;
K_posterior = K_posterior + 1e-9 * eye(n_grid);

% Draw samples from posterior: f*|y ~ N(μ*, Σ*)
L_posterior = chol(K_posterior, 'lower');
samples_posterior = zeros(n_grid, n_samples);
for i = 1:n_samples
    z = randn(n_grid, 1);
    samples_posterior(:, i) = mu_posterior + L_posterior * z;
end

% Compute uncertainty bands (mean ± 2*std)
std_posterior = sqrt(diag(K_posterior));
mu_upper_posterior = mu_posterior + 2 * std_posterior;
mu_lower_posterior = mu_posterior - 2 * std_posterior;

%% ============================================
%% VISUALIZATION
%% ============================================
fprintf('Creating plots...\n');

figure('Color', 'w', 'Position', [100, 100, 1200, 500]);

% Plot (a): Prior
subplot(1, 2, 1);
hold on;

% Shaded uncertainty region (95% confidence: mean ± 2*std)
fill([x_grid; flipud(x_grid)], [mu_upper_prior; flipud(mu_lower_prior)], ...
    [0.85, 0.85, 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.7);

% Draw sample functions
colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], [0.5 0.2 0.8]};
for i = 1:n_samples
    plot(x_grid, samples_prior(:, i), '-', 'Color', colors{i}, ...
        'LineWidth', 1.5, 'DisplayName', sprintf('Sample %d', i));
end

xlabel('input, x', 'FontSize', 12);
ylabel('output, f(x)', 'FontSize', 12);
title('(a) prior', 'FontSize', 12, 'FontWeight', 'bold');
xlim([x_min, x_max]);
ylim([-3, 3]);
grid on;
box on;
legend('Location', 'best', 'FontSize', 10);

% Plot (b): Posterior
subplot(1, 2, 2);
hold on;

% Shaded uncertainty region (95% confidence: mean ± 2*std)
fill([x_grid; flipud(x_grid)], [mu_upper_posterior; flipud(mu_lower_posterior)], ...
    [0.85, 0.85, 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.7);

% Draw sample functions
for i = 1:n_samples
    plot(x_grid, samples_posterior(:, i), '-', 'Color', colors{i}, ...
        'LineWidth', 1.5, 'DisplayName', sprintf('Sample %d', i));
end

% Plot training data points
plot(x_train, y_train, 'k+', 'MarkerSize', 12, 'LineWidth', 2, ...
    'DisplayName', 'Observations');

xlabel('input, x', 'FontSize', 12);
ylabel('output, f(x)', 'FontSize', 12);
title('(b) posterior', 'FontSize', 12, 'FontWeight', 'bold');
xlim([x_min, x_max]);
ylim([-3, 3]);
grid on;
box on;
legend('Location', 'best', 'FontSize', 10);

sgtitle('Gaussian Process: Prior vs Posterior', 'FontSize', 14, 'FontWeight', 'bold');

%% Print summary
fprintf('\nSummary:\n');
fprintf('  Grid points: %d (from %.1f to %.1f)\n', n_grid, x_min, x_max);
fprintf('  Training points: %d\n', n_train);
fprintf('  Prior std range: [%.4f, %.4f]\n', min(std_prior), max(std_prior));
fprintf('  Posterior std range: [%.4f, %.4f]\n', min(std_posterior), max(std_posterior));
% Find grid points closest to training points for uncertainty comparison
[~, train_indices] = min(abs(x_grid - x_train'), [], 1);
fprintf('  Uncertainty reduction: %.1f%% (at training points)\n', ...
    100 * (1 - mean(std_posterior(train_indices)) / mean(std_prior(train_indices))));

fprintf('\nDone!\n');
