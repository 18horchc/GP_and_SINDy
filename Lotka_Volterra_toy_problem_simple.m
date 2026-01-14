%% Lotka-Volterra Gaussian Process Learning
% Based on: Ye, Dongwei, and Mengwu Guo. 2024. "Gaussian Process Learning 
% of Nonlinear Dynamics." Communications in Nonlinear Science and Numerical 
% Simulation 138 (November): 108184.
%
% This script:
% 1. Solves Lotka-Volterra system using implicit Euler method
% 2. Splits data into training [0, 0.4T] and testing (0.4T, T]
% 3. Randomly samples 25% of training data
% 4. Fits GP with squared exponential kernel, optimizing hyperparameters

clear; close all; clc;

%% System Parameters
alpha = 1.5;
beta = 1.0;
delta = 1.0;
gamma = 3.0;

fprintf('Lotka-Volterra System Parameters:\n');
fprintf('  alpha=%.2f, beta=%.2f, delta=%.2f, gamma=%.2f\n', alpha, beta, delta, gamma);

%% Numerical Integration Parameters
T = 20;           % Total time
dt = 0.001;       % Desired time step for output
t_all = 0:dt:T;   % Time vector for output

% Initial conditions
x0 = [1; 1];      % x1(0) = x2(0) = 1

fprintf('\nIntegration Parameters:\n');
fprintf('  Time span: [0, %.1f]\n', T);
fprintf('  Output time step: %.4f\n', dt);
fprintf('  Initial conditions: x1(0)=%.1f, x2(0)=%.1f\n', x0(1), x0(2));

%% Solve using ODE Solver (ode45)
fprintf('\nSolving ODE using ode45...\n');

% Lotka-Volterra ODE right-hand side
lv_ode = @(t, x) [alpha * x(1) - beta * x(1) * x(2);
                  delta * x(1) * x(2) - gamma * x(2)];

% Solve ODE using ode45 (adaptive step size, more accurate and faster)
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[t_sol, X_sol] = ode45(lv_ode, [0, T], x0, options);

% Interpolate to desired time points
X = interp1(t_sol, X_sol, t_all, 'pchip');
t_all = t_all(:);  % Ensure column vector

fprintf('  Integration complete.\n');
fprintf('  Number of output points: %d\n', length(t_all));
fprintf('  State range: x1 ∈ [%.3f, %.3f], x2 ∈ [%.3f, %.3f]\n', ...
    min(X(:,1)), max(X(:,1)), min(X(:,2)), max(X(:,2)));

%% Split Data into Training and Testing
T_split = 0.4 * T;  % Split point: 0.4T = 8
train_idx = t_all <= T_split;
test_idx = t_all > T_split;

t_train_full = t_all(train_idx);
X_train_full = X(train_idx, :);

t_test = t_all(test_idx);
X_test = X(test_idx, :);

fprintf('\nData Splitting:\n');
fprintf('  Training set: [0, %.1f], %d points\n', T_split, length(t_train_full));
fprintf('  Testing set: (%.1f, %.1f], %d points\n', T_split, T, length(t_test));

%% Randomly Sample 25% of Training Data
rng(42);  % Set seed for reproducibility
n_train_full = length(t_train_full);
n_samples = round(0.25 * n_train_full);  % 25% of training data

sample_idx = randperm(n_train_full, n_samples);
sample_idx = sort(sample_idx);  % Keep time order

t_train = t_train_full(sample_idx);
X_train = X_train_full(sample_idx, :);

fprintf('\nSampling:\n');
fprintf('  Sampled %d points (25%% of %d training points)\n', n_samples, n_train_full);

%% Add Noise (0% for now, can be adjusted later)
noise_level = 0.0;  % 0% noise
X_train_noisy = X_train + noise_level * randn(size(X_train));
X_train_noisy = max(0, X_train_noisy);  % Ensure non-negative

fprintf('  Noise level: %.1f%%\n', noise_level * 100);

%% Fit Gaussian Process Models
fprintf('\nFitting GP models with squared exponential kernel...\n');
fprintf('  Optimizing hyperparameters by maximizing negative log likelihood...\n');

% Variable names
var_names = {'x1 (Prey)', 'x2 (Predator)'};

% Initialize storage for GP models and predictions
gp_models = cell(2, 1);
X_gp_all = zeros(length(t_all), 2);
X_gp_lower = zeros(length(t_all), 2);
X_gp_upper = zeros(length(t_all), 2);

% Fit GP for each variable
for i = 1:2
    y_train = X_train_noisy(:, i);
    
    fprintf('\n  Fitting GP for %s...\n', var_names{i});
    
    % Fit GP model with squared exponential kernel
    % fitrgp optimizes hyperparameters by maximizing log likelihood by default
    gp_models{i} = fitrgp(t_train, y_train, ...
        'KernelFunction', 'squaredexponential', ...
        'BasisFunction', 'constant', ...
        'Standardize', false);
    
    % Predict on full time range with confidence intervals
    [y_pred, ~, y_int] = predict(gp_models{i}, t_all);
    X_gp_all(:, i) = y_pred;
    X_gp_lower(:, i) = y_int(:, 1);  % Lower bound of 95% prediction interval
    X_gp_upper(:, i) = y_int(:, 2);  % Upper bound of 95% prediction interval
    
    % Extract hyperparameters
    kernel_params = gp_models{i}.KernelInformation.KernelParameters;
    lengthscale = kernel_params(1);
    signal_std = kernel_params(2);
    noise_std = gp_models{i}.Sigma;
    
    % Compute R² on training data
    y_pred_train = predict(gp_models{i}, t_train);
    ss_res = sum((y_train - y_pred_train).^2);
    ss_tot = sum((y_train - mean(y_train)).^2);
    r2_train = 1 - ss_res / ss_tot;
    
    fprintf('    Hyperparameters:\n');
    fprintf('      Lengthscale: %.4f\n', lengthscale);
    fprintf('      Signal std: %.4f\n', signal_std);
    fprintf('      Noise std: %.4f\n', noise_std);
    fprintf('    Training R²: %.4f\n', r2_train);
end

%% Visualization
fprintf('\nGenerating visualizations...\n');

figure('Color', 'w', 'Position', [100, 100, 1200, 600]);
hold on;

% Ground truth: x1 and x2 - black solid lines
plot(t_all, X(:,1), '-', 'Color', [0, 0, 0], 'LineWidth', 2, ...
    'DisplayName', 'Ground Truth');
plot(t_all, X(:,2), '-', 'Color', [0, 0, 0], 'LineWidth', 2, ...
    'HandleVisibility', 'off');  % Same style, don't duplicate in legend

% Uncertainty band for x1 (Prey) - blue
fill([t_all; flipud(t_all)], ...
     [X_gp_upper(:,1); flipud(X_gp_lower(:,1))], ...
     [0, 0.4470, 0.7410], 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
     'DisplayName', 'GP 95% CI (x_1)');

% Uncertainty band for x2 (Predator) - orange/red
fill([t_all; flipud(t_all)], ...
     [X_gp_upper(:,2); flipud(X_gp_lower(:,2))], ...
     [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
     'DisplayName', 'GP 95% CI (x_2)');

% GP mean: x1 (Prey) - dotted blue line
plot(t_all, X_gp_all(:,1), ':', 'Color', [0, 0.4470, 0.7410], ...
    'LineWidth', 2, 'DisplayName', 'GP Mean (x_1)');

% GP mean: x2 (Predator) - dotted orange/red line
plot(t_all, X_gp_all(:,2), ':', 'Color', [0.8500, 0.3250, 0.0980], ...
    'LineWidth', 2, 'DisplayName', 'GP Mean (x_2)');

% Training data: x1 (Prey) - blue markers
scatter(t_train, X_train_noisy(:,1), 80, [0, 0.4470, 0.7410], ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1, ...
    'DisplayName', 'Training Data (x_1)');

% Training data: x2 (Predator) - orange/red markers
scatter(t_train, X_train_noisy(:,2), 80, [0.8500, 0.3250, 0.0980], ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1, ...
    'DisplayName', 'Training Data (x_2)');

% Train/test split line
xline(T_split, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Train/Test Split');

xlabel('t', 'FontSize', 12);
ylabel('x_i(t)', 'FontSize', 12);
title(sprintf('Lotka-Volterra GP Learning (Paper Implementation)\nTraining: [0, %.1f], Testing: (%.1f, %.1f]', ...
    T_split, T_split, T), 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
box on;

fprintf('\nVisualization complete!\n');
