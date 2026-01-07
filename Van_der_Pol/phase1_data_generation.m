%% Phase 1: Data Generation for Van der Pol Oscillator
% This script generates ground truth solution and sparse noisy data
% for testing GP+SINDy approach

clear; close all; clc;

%% Parameters
% Van der Pol parameters
mu = 1.0;  % Moderate nonlinearity
x0 = [2; 0];  % Initial conditions

% Time span
t_start = 0;
t_end = 10;
dt_fine = 0.01;  % Fine resolution for ground truth
t_span = t_start:dt_fine:t_end;

% Sparse data parameters
n_sparse_points = 9;  % 8-10 points (using 9)
noise_fraction = 0.1;  % 10% noise
sampling_method = 'uniform';  % 'uniform', 'random', or 'strategic'

%% Generate Ground Truth
fprintf('Generating Van der Pol ground truth solution...\n');
fprintf('  Parameters: mu = %.2f, x0 = [%.2f; %.2f]\n', mu, x0(1), x0(2));
fprintf('  Time span: [%.2f, %.2f] with dt = %.3f\n', t_start, t_end, dt_fine);

[t_true, X_true, dXdt_true, params] = generate_vdp_ground_truth(t_span, x0, mu);

fprintf('  Generated %d time points\n', length(t_true));
fprintf('  State range: x ∈ [%.2f, %.2f], y ∈ [%.2f, %.2f]\n', ...
    min(X_true(:,1)), max(X_true(:,1)), min(X_true(:,2)), max(X_true(:,2)));

%% Generate Sparse Noisy Data
fprintf('\nGenerating sparse noisy data...\n');
fprintf('  Number of points: %d\n', n_sparse_points);
fprintf('  Noise level: %.1f%%\n', noise_fraction * 100);
fprintf('  Sampling method: %s\n', sampling_method);

[t_sparse, X_sparse, X_true_at_sparse, noise_level] = ...
    generate_sparse_noisy_data(t_true, X_true, n_sparse_points, noise_fraction, sampling_method);

fprintf('  Generated %d sparse points\n', length(t_sparse));
fprintf('  Noise std: [%.4f, %.4f]\n', noise_level(1), noise_level(2));

%% Visualization: Data Overview
fprintf('\nCreating visualization...\n');

figure('Color', 'w', 'Position', [100, 100, 1200, 400]);

% Plot 1: x(t)
subplot(1, 2, 1);
plot(t_true, X_true(:,1), '-', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2, ...
    'DisplayName', 'True x(t)'); hold on;
scatter(t_sparse, X_sparse(:,1), 80, [0.8500, 0.3250, 0.0980], 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', 'Noisy samples');
plot(t_sparse, X_true_at_sparse(:,1), 'o', 'Color', [0, 0.4470, 0.7410], ...
    'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'True at samples');
xlabel('Time'); ylabel('x(t)');
title(sprintf('Van der Pol: x(t) - %d sparse points, %.0f%% noise', ...
    n_sparse_points, noise_fraction*100));
legend('Location', 'best'); grid on;

% Plot 2: y(t)
subplot(1, 2, 2);
plot(t_true, X_true(:,2), '-', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2, ...
    'DisplayName', 'True y(t)'); hold on;
scatter(t_sparse, X_sparse(:,2), 80, [0, 0.4470, 0.7410], 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', 'Noisy samples');
plot(t_sparse, X_true_at_sparse(:,2), 'o', 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'True at samples');
xlabel('Time'); ylabel('y(t)');
title(sprintf('Van der Pol: y(t) - %d sparse points, %.0f%% noise', ...
    n_sparse_points, noise_fraction*100));
legend('Location', 'best'); grid on;

sgtitle('Phase 1: Ground Truth & Sparse Noisy Data', 'FontSize', 14, 'FontWeight', 'bold');

%% Save data for next phases
fprintf('\nSaving data...\n');
save('vdp_data.mat', 't_true', 'X_true', 'dXdt_true', 't_sparse', 'X_sparse', ...
    'X_true_at_sparse', 'params', 'noise_level', 'n_sparse_points', 'noise_fraction');

fprintf('Data saved to vdp_data.mat\n');
fprintf('\nPhase 1 complete! Ready for Phase 2 (GP Implementation).\n');

