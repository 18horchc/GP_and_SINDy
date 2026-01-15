%% Lotka-Volterra Gaussian Process Learning
% Based on: Ye, Dongwei, and Mengwu Guo. 2024. "Gaussian Process Learning 
% of Nonlinear Dynamics." Communications in Nonlinear Science and Numerical 
% Simulation 138 (November): 108184.
%
% This script implements Algorithm 1 from Ye & Guo (2024):
% 1. Solves Lotka-Volterra system using ode45
% 2. Splits data into training [0, 0.4T] and testing (0.4T, T]
% 3. Randomly samples 25% of training data
% 4. Estimates GP hyperparameters via maximum likelihood
% 5. Computes covariance matrices (K^uu, K^dd, K^du)
% 6. Smooths state data and estimates time derivatives
% 7. Estimates model parameters using weighted least squares with regularization

clear; close all; clc;

%% System Parameters
alpha = 1.5;
beta = 1.0;
delta = 1.0;
gamma = 3.0;

%% Experiment Configuration (match Figure 1 and Table 1)
% Noise levels are percentages of signal std (e.g., 0.1 = 10% noise)
noise_levels = [0.0, 0.10, 0.20];

% Sample counts for training data (absolute counts, not fractions)
sample_counts = [200, 100, 20];

% Posterior sampling for extrapolation uncertainty bands
posterior_samples = 200;        % Increase for smoother bands
posterior_max_state = 20;       % Discard unstable trajectories

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

%% Run Experiments (Table 1 + Figure 1-style visualization)
fprintf('\n=== Running LV GP experiments ===\n');

results = [];
num_cases = numel(noise_levels) * numel(sample_counts);
plot_results = cell(1, num_cases);
case_counter = 0;

for noise_level = noise_levels
    for sample_count = sample_counts
        case_counter = case_counter + 1;
        [case_result, plot_result] = run_lv_gp_case( ...
            t_all, X, t_train_full, X_train_full, t_test, X_test, T_split, ...
            alpha, beta, delta, gamma, ...
            noise_level, sample_count, posterior_samples, posterior_max_state, case_counter);
        results = [results; case_result];
        plot_results{case_counter} = plot_result;
    end
end

results_table = struct2table(results);
disp(results_table);

%% Visualization: subplots for all noise/sparsity cases
fprintf('\nGenerating subplot grid for all cases...\n');
num_rows = numel(noise_levels);
num_cols = numel(sample_counts);

figure('Color', 'w', 'Position', [100, 100, 1600, 900]);
tiledlayout(num_rows, num_cols, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(plot_results)
    pr = plot_results{i};
    nexttile;
    hold on;

    % Ground truth
    plot(t_all, X(:,1), '-', 'Color', [0, 0, 0], 'LineWidth', 1.5, ...
        'DisplayName', 'Ground Truth');
    plot(t_all, X(:,2), '-', 'Color', [0, 0, 0], 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');

    % GP mean/CI on training interval only
    fill([pr.t_train_all; flipud(pr.t_train_all)], ...
         [pr.gp_upper_train(:,1); flipud(pr.gp_lower_train(:,1))], ...
         [0, 0.4470, 0.7410], 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
         'DisplayName', 'GP 95% CI (x_1)');
    fill([pr.t_train_all; flipud(pr.t_train_all)], ...
         [pr.gp_upper_train(:,2); flipud(pr.gp_lower_train(:,2))], ...
         [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
         'DisplayName', 'GP 95% CI (x_2)');

    plot(pr.t_train_all, pr.gp_mean_train(:,1), '-', 'Color', [0, 0.4470, 0.7410], ...
        'LineWidth', 1.5, 'DisplayName', 'GP Mean (x_1, train)');
    plot(pr.t_train_all, pr.gp_mean_train(:,2), '-', 'Color', [0.8500, 0.3250, 0.0980], ...
        'LineWidth', 1.5, 'DisplayName', 'GP Mean (x_2, train)');

    % Forward simulation mean ± std (extrapolation band)
    fill([t_all; flipud(t_all)], ...
         [pr.pred_mean(:,1) + pr.pred_std(:,1); flipud(pr.pred_mean(:,1) - pr.pred_std(:,1))], ...
         [0, 0.4470, 0.7410], 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
         'DisplayName', 'Pred ±1σ (x_1)');
    fill([t_all; flipud(t_all)], ...
         [pr.pred_mean(:,2) + pr.pred_std(:,2); flipud(pr.pred_mean(:,2) - pr.pred_std(:,2))], ...
         [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.12, 'EdgeColor', 'none', ...
         'DisplayName', 'Pred ±1σ (x_2)');

    plot(t_all, pr.pred_mean(:,1), '--', 'Color', [0, 0.4470, 0.7410], ...
        'LineWidth', 1.2, 'DisplayName', 'Sim (x_1, est params)');
    plot(t_all, pr.pred_mean(:,2), '--', 'Color', [0.8500, 0.3250, 0.0980], ...
        'LineWidth', 1.2, 'DisplayName', 'Sim (x_2, est params)');

    % Training points and split
    scatter(pr.t_train, pr.X_train_noisy(:,1), 30, [0, 0.4470, 0.7410], ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, ...
        'DisplayName', 'Training Data (x_1)');
    scatter(pr.t_train, pr.X_train_noisy(:,2), 30, [0.8500, 0.3250, 0.0980], ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, ...
        'DisplayName', 'Training Data (x_2)');

    xline(T_split, 'r--', 'LineWidth', 1.0, 'DisplayName', 'Train/Test Split');

    title(sprintf('Noise %.0f%% | N=%d', pr.noise_level*100, pr.sample_count), ...
        'FontSize', 10, 'FontWeight', 'bold');
    set(gca, 'FontSize', 8);
    grid on;
    box on;

    if i == 1
        legend('Location', 'best', 'FontSize', 7);
    else
        legend('off');
    end
end

%% Local function: run one LV GP case
function [case_result, plot_result] = run_lv_gp_case( ...
    t_all, X, t_train_full, X_train_full, t_test, X_test, T_split, ...
    alpha, beta, delta, gamma, noise_level, sample_count, posterior_samples, posterior_max_state, case_counter)

% Variable names
var_names = {'x1 (Prey)', 'x2 (Predator)'};

% Reproducible sampling for each case
rng(42);
n_train_full = length(t_train_full);
n_samples = min(n_train_full, max(5, sample_count));
sample_idx = randperm(n_train_full, n_samples);
sample_idx = sort(sample_idx);
t_train = t_train_full(sample_idx);
X_train = X_train_full(sample_idx, :);

fprintf('\nCase %d: noise=%.1f%%, sample=%.2f%% (%d/%d)\n', ...
    case_counter, noise_level * 100, 100 * n_samples / n_train_full, n_samples, n_train_full);

% Add noise proportional to signal magnitude
X_train_noisy = X_train + noise_level * std(X_train, 0, 1) .* randn(size(X_train));
X_train_noisy = max(0, X_train_noisy);  % Ensure non-negative

% Fit GP models
fprintf('  Fitting GP models with squared exponential kernel...\n');

gp_models = cell(2, 1);
X_gp_all = zeros(length(t_all), 2);
X_gp_lower = zeros(length(t_all), 2);
X_gp_upper = zeros(length(t_all), 2);

min_noise_std = 0.002;
r2_train_all = zeros(2, 1);
r2_test_all = zeros(2, 1);

for i = 1:2
    y_train = X_train_noisy(:, i);

    fprintf('    GP for %s...\n', var_names{i});

    y_std = std(y_train);
    t_range = max(t_train) - min(t_train);
    initial_lengthscale = 0.3 * t_range;
    initial_signal_std = y_std;
    initial_noise_std = max(min_noise_std, noise_level * y_std);

    gp_models{i} = fitrgp(t_train, y_train, ...
        'KernelFunction', 'squaredexponential', ...
        'BasisFunction', 'constant', ...
        'Standardize', false, ...
        'KernelParameters', [initial_lengthscale, initial_signal_std], ...
        'Sigma', initial_noise_std, ...
        'SigmaLowerBound', max(1e-6, initial_noise_std - 0.002));

    kernel_params = gp_models{i}.KernelInformation.KernelParameters;
    lengthscale = kernel_params(1);
    signal_std = kernel_params(2);
    noise_std = gp_models{i}.Sigma;

    if noise_std < min_noise_std
        gp_models{i} = fitrgp(t_train, y_train, ...
            'KernelFunction', 'squaredexponential', ...
            'BasisFunction', 'constant', ...
            'Standardize', false, ...
            'KernelParameters', [lengthscale, signal_std], ...
            'Sigma', min_noise_std, ...
            'SigmaLowerBound', max(1e-6, min_noise_std - 0.002));
        noise_std = gp_models{i}.Sigma;
    end

    [y_pred, ~, y_int] = predict(gp_models{i}, t_all);
    X_gp_all(:, i) = y_pred;
    X_gp_lower(:, i) = y_int(:, 1);
    X_gp_upper(:, i) = y_int(:, 2);

    y_pred_train = predict(gp_models{i}, t_train);
    ss_res_train = sum((y_train - y_pred_train).^2);
    ss_tot_train = sum((y_train - mean(y_train)).^2);
    r2_train = 1 - ss_res_train / ss_tot_train;

    y_test = X_test(:, i);
    y_pred_test = predict(gp_models{i}, t_test);
    ss_res_test = sum((y_test - y_pred_test).^2);
    ss_tot_test = sum((y_test - mean(y_test)).^2);
    r2_test = 1 - ss_res_test / ss_tot_test;

    r2_train_all(i) = r2_train;
    r2_test_all(i) = r2_test;

    fprintf('      Lengthscale: %.4f | Signal std: %.4f | Noise std: %.4f\n', ...
        lengthscale, signal_std, noise_std);
    fprintf('      R2 train: %.4f | R2 test: %.4f\n', r2_train, r2_test);
end

% Algorithm 1: estimate derivatives and parameters
N = 2;
smoothed_states = zeros(length(t_all), N);
estimated_derivatives = zeros(length(t_all), N);
R_dd_inv = cell(N, 1);

for i = 1:N
    kernel_params = gp_models{i}.KernelInformation.KernelParameters;
    rho_i = kernel_params;
    chi_u = gp_models{i}.Sigma^2;

    T_vec = t_train;
    K_train = size(T_vec, 1);
    l = rho_i(1);
    sigma_f = rho_i(2);

    T_mat = repmat(T_vec, 1, K_train);
    dT = T_mat - T_mat';
    K_uu = sigma_f^2 * exp(-dT.^2 / (2*l^2));
    K_uu = K_uu + chi_u * eye(K_train);

    % Add small jitter to improve conditioning
    jitter_uu = max(1e-10, 1e-6 * trace(K_uu) / K_train);
    K_uu = K_uu + jitter_uu * eye(K_train);

    K_dd = (sigma_f^2 / l^2) * (1 - dT.^2 / l^2) .* exp(-dT.^2 / (2*l^2));
    K_dd = K_dd + chi_u * eye(K_train);
    jitter_dd = max(1e-10, 1e-6 * trace(K_dd) / K_train);
    K_dd = K_dd + jitter_dd * eye(K_train);

    K_du = -(sigma_f^2 / l^2) * dT .* exp(-dT.^2 / (2*l^2));

    % Avoid explicit matrix inverse
    K_uu_inv = K_uu \ eye(K_train);
    Schur = K_dd - K_du * K_uu_inv * K_du';
    Schur = (Schur + Schur') * 0.5;
    R_dd_inv{i} = Schur \ eye(K_train);

    T_all_mat = repmat(t_all, 1, K_train);
    T_train_mat = repmat(T_vec', length(t_all), 1);
    dT_pred = T_all_mat - T_train_mat;
    K_Tall_T = sigma_f^2 * exp(-dT_pred.^2 / (2*l^2));

    y_train_i = X_train_noisy(:, i);
    smoothed_states(:, i) = K_Tall_T * (K_uu_inv * y_train_i);
    K_du_Tall_T = -(sigma_f^2 / l^2) * dT_pred .* exp(-dT_pred.^2 / (2*l^2));
    estimated_derivatives(:, i) = K_du_Tall_T * (K_uu_inv * y_train_i);
end

X1_smooth = smoothed_states(:, 1);
X2_smooth = smoothed_states(:, 2);
dX1_est = estimated_derivatives(:, 1);
dX2_est = estimated_derivatives(:, 2);

lambda_reg = 1e-4;
X1_train_smooth = interp1(t_all, X1_smooth, t_train, 'linear');
X2_train_smooth = interp1(t_all, X2_smooth, t_train, 'linear');
dX1_train = interp1(t_all, dX1_est, t_train, 'linear');
dX2_train = interp1(t_all, dX2_est, t_train, 'linear');

Theta1_train = [X1_train_smooth, X1_train_smooth .* X2_train_smooth];
Theta2_train = [X1_train_smooth .* X2_train_smooth, X2_train_smooth];

R1_train = R_dd_inv{1};
R2_train = R_dd_inv{2};

Theta1_mat = Theta1_train' * R1_train * Theta1_train + lambda_reg * eye(2);
Theta2_mat = Theta2_train' * R2_train * Theta2_train + lambda_reg * eye(2);
Theta1_mat = (Theta1_mat + Theta1_mat') * 0.5;
Theta2_mat = (Theta2_mat + Theta2_mat') * 0.5;

theta1_cov = Theta1_mat \ eye(2);
theta2_cov = Theta2_mat \ eye(2);
theta1_est = theta1_cov * (Theta1_train' * R1_train * dX1_train);
theta2_est = theta2_cov * (Theta2_train' * R2_train * dX2_train);

alpha_est = theta1_est(1);
beta_est = -theta1_est(2);
delta_est = theta2_est(1);
gamma_est = -theta2_est(2);

rel_err_alpha = 100 * abs(alpha_est - alpha) / abs(alpha);
rel_err_beta = 100 * abs(beta_est - beta) / abs(beta);
rel_err_delta = 100 * abs(delta_est - delta) / abs(delta);
rel_err_gamma = 100 * abs(gamma_est - gamma) / abs(gamma);

fprintf('  Params: alpha=%.4f, beta=%.4f, delta=%.4f, gamma=%.4f\n', ...
    alpha_est, beta_est, delta_est, gamma_est);

case_result = struct( ...
    'noise_level', noise_level, ...
    'sample_count', n_samples, ...
    'n_samples', n_samples, ...
    'r2_train_x1', r2_train_all(1), ...
    'r2_train_x2', r2_train_all(2), ...
    'r2_test_x1', r2_test_all(1), ...
    'r2_test_x2', r2_test_all(2), ...
    'alpha_est', alpha_est, ...
    'beta_est', beta_est, ...
    'delta_est', delta_est, ...
    'gamma_est', gamma_est, ...
    'rel_err_alpha', rel_err_alpha, ...
    'rel_err_beta', rel_err_beta, ...
    'rel_err_delta', rel_err_delta, ...
    'rel_err_gamma', rel_err_gamma);

% Posterior sampling for extrapolation uncertainty
num_t = length(t_all);
traj_store = zeros(num_t, 2, posterior_samples);
valid_count = 0;

% Ensure covariance is positive definite for sampling
theta1_cov = (theta1_cov + theta1_cov') * 0.5;
theta2_cov = (theta2_cov + theta2_cov') * 0.5;

[L1, p1] = chol(theta1_cov, 'lower');
if p1 > 0
    jitter1 = max(1e-10, 1e-6 * trace(theta1_cov) / 2);
    [L1, p1] = chol(theta1_cov + jitter1 * eye(2), 'lower');
end

[L2, p2] = chol(theta2_cov, 'lower');
if p2 > 0
    jitter2 = max(1e-10, 1e-6 * trace(theta2_cov) / 2);
    [L2, p2] = chol(theta2_cov + jitter2 * eye(2), 'lower');
end

if p1 > 0 || p2 > 0
    warning('Posterior covariance not PD; using diagonal approximation for sampling.');
    L1 = diag(sqrt(max(diag(theta1_cov), 0)));
    L2 = diag(sqrt(max(diag(theta2_cov), 0)));
end

for s = 1:posterior_samples
    theta1_s = theta1_est + L1 * randn(2, 1);
    theta2_s = theta2_est + L2 * randn(2, 1);

    alpha_s = theta1_s(1);
    beta_s = -theta1_s(2);
    delta_s = theta2_s(1);
    gamma_s = -theta2_s(2);

    lv_ode_est = @(t, x) [alpha_s * x(1) - beta_s * x(1) * x(2); ...
                          delta_s * x(1) * x(2) - gamma_s * x(2)];
    [t_est, X_est] = ode45(lv_ode_est, [0, t_all(end)], X(1, :)', ...
        odeset('RelTol', 1e-8, 'AbsTol', 1e-10));
    X_est = interp1(t_est, X_est, t_all, 'pchip');

    if any(~isfinite(X_est(:))) || max(X_est(:)) > posterior_max_state
        continue;
    end

    valid_count = valid_count + 1;
    traj_store(:, :, valid_count) = X_est;
end

if valid_count == 0
    pred_mean = X_est;
    pred_std = zeros(num_t, 2);
else
    traj_store = traj_store(:, :, 1:valid_count);
    pred_mean = mean(traj_store, 3);
    pred_std = std(traj_store, 0, 3);
end

% GP train interval arrays for plotting
train_mask = t_all <= T_split;

plot_result = struct( ...
    'noise_level', noise_level, ...
    'sample_count', n_samples, ...
    't_train', t_train, ...
    'X_train_noisy', X_train_noisy, ...
    't_train_all', t_all(train_mask), ...
    'gp_mean_train', X_gp_all(train_mask, :), ...
    'gp_lower_train', X_gp_lower(train_mask, :), ...
    'gp_upper_train', X_gp_upper(train_mask, :), ...
    'pred_mean', pred_mean, ...
    'pred_std', pred_std);

end
