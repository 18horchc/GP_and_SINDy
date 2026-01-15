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
noise_levels = [0.10, 0.20];

% Sample counts for training data (absolute counts, not fractions)
sample_counts = [200, 100, 20];

% Select one case for the Figure 1-style plot
plot_case.noise_level = 0.10;
plot_case.sample_count = 20;

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
case_counter = 0;

for noise_level = noise_levels
    for sample_count = sample_counts
        case_counter = case_counter + 1;
        do_plot = (abs(noise_level - plot_case.noise_level) < 1e-12) && ...
                  (sample_count == plot_case.sample_count);
        case_result = run_lv_gp_case( ...
            t_all, X, t_train_full, X_train_full, t_test, X_test, T_split, ...
            alpha, beta, delta, gamma, ...
            noise_level, sample_count, do_plot, case_counter);
        results = [results; case_result];
    end
end

results_table = struct2table(results);
disp(results_table);

%% Local function: run one LV GP case
function case_result = run_lv_gp_case( ...
    t_all, X, t_train_full, X_train_full, t_test, X_test, T_split, ...
    alpha, beta, delta, gamma, noise_level, sample_count, do_plot, case_counter)

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

    K_dd = (sigma_f^2 / l^2) * (1 - dT.^2 / l^2) .* exp(-dT.^2 / (2*l^2));
    K_dd = K_dd + chi_u * eye(K_train);

    K_du = -(sigma_f^2 / l^2) * dT .* exp(-dT.^2 / (2*l^2));
    K_uu_inv = inv(K_uu);
    R_dd_inv{i} = inv(K_dd - K_du * K_uu_inv * K_du');

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

theta1_est = (Theta1_train' * R1_train * Theta1_train + lambda_reg * eye(2)) \ ...
             (Theta1_train' * R1_train * dX1_train);
theta2_est = (Theta2_train' * R2_train * Theta2_train + lambda_reg * eye(2)) \ ...
             (Theta2_train' * R2_train * dX2_train);

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

if do_plot
    fprintf('  Generating Figure 1-style plot...\n');
    figure('Color', 'w', 'Position', [100, 100, 1200, 600]);
    hold on;

    plot(t_all, X(:,1), '-', 'Color', [0, 0, 0], 'LineWidth', 2, ...
        'DisplayName', 'Ground Truth');
    plot(t_all, X(:,2), '-', 'Color', [0, 0, 0], 'LineWidth', 2, ...
        'HandleVisibility', 'off');

    train_mask = t_all <= T_split;
    t_train_all = t_all(train_mask);
    X_gp_upper_train = X_gp_upper(train_mask, :);
    X_gp_lower_train = X_gp_lower(train_mask, :);
    X_gp_mean_train = X_gp_all(train_mask, :);

    fill([t_train_all; flipud(t_train_all)], ...
         [X_gp_upper_train(:,1); flipud(X_gp_lower_train(:,1))], ...
         [0, 0.4470, 0.7410], 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
         'DisplayName', 'GP 95% CI (x_1)');
    fill([t_train_all; flipud(t_train_all)], ...
         [X_gp_upper_train(:,2); flipud(X_gp_lower_train(:,2))], ...
         [0.8500, 0.3250, 0.0980], 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
         'DisplayName', 'GP 95% CI (x_2)');

    plot(t_train_all, X_gp_mean_train(:,1), '-', 'Color', [0, 0.4470, 0.7410], ...
        'LineWidth', 2, 'DisplayName', 'GP Mean (x_1, train)');
    plot(t_train_all, X_gp_mean_train(:,2), '-', 'Color', [0.8500, 0.3250, 0.0980], ...
        'LineWidth', 2, 'DisplayName', 'GP Mean (x_2, train)');

    scatter(t_train, X_train_noisy(:,1), 80, [0, 0.4470, 0.7410], ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1, ...
        'DisplayName', 'Training Data (x_1)');
    scatter(t_train, X_train_noisy(:,2), 80, [0.8500, 0.3250, 0.0980], ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1, ...
        'DisplayName', 'Training Data (x_2)');

    xline(T_split, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Train/Test Split');

    % Forward simulate using estimated parameters for full horizon
    x0 = X(1, :)';
    lv_ode_est = @(t, x) [alpha_est * x(1) - beta_est * x(1) * x(2); ...
                          delta_est * x(1) * x(2) - gamma_est * x(2)];
    [t_est, X_est] = ode45(lv_ode_est, [0, t_all(end)], x0, ...
        odeset('RelTol', 1e-8, 'AbsTol', 1e-10));
    X_est = interp1(t_est, X_est, t_all, 'pchip');

    plot(t_all, X_est(:,1), '--', 'Color', [0, 0.4470, 0.7410], ...
        'LineWidth', 1.5, 'DisplayName', 'Sim (x_1, est params)');
    plot(t_all, X_est(:,2), '--', 'Color', [0.8500, 0.3250, 0.0980], ...
        'LineWidth', 1.5, 'DisplayName', 'Sim (x_2, est params)');

    xlabel('t', 'FontSize', 12);
    ylabel('x_i(t)', 'FontSize', 12);
    title(sprintf('Lotka-Volterra GP Learning (Noise %.1f%%, Sample %.2f%%)', ...
        noise_level * 100, 100 * n_samples / n_train_full), 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9);
    grid on;
    box on;

    param_str = sprintf(['Estimated Parameters:\n' ...
        'α=%.3f (true: %.2f)\nβ=%.3f (true: %.2f)\n' ...
        'δ=%.3f (true: %.2f)\nγ=%.3f (true: %.2f)'], ...
        alpha_est, alpha, beta_est, beta, delta_est, delta, gamma_est, gamma);
    text(0.02, 0.98, param_str, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'FontSize', 10, ...
        'BackgroundColor', [1, 1, 1, 0.8], 'EdgeColor', 'k');
end

end
