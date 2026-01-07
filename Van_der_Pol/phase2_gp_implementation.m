%% Phase 2: GP Implementation for Van der Pol Oscillator
% This script fits Gaussian Process models to sparse noisy data,
% predicts on dense time grid, and computes analytical derivatives
%
% Requirements:
%   - Phase 1 must be run first (generates vdp_data.mat)
%   - Statistics and Machine Learning Toolbox (for fitrgp)

clear; close all; clc;

%% Check if Phase 1 data exists
if ~exist('vdp_data.mat', 'file')
    error('Phase 1 data not found. Please run phase1_data_generation.m first.');
end

%% Load Phase 1 data
fprintf('Loading Phase 1 data...\n');
load('vdp_data.mat', 't_true', 'X_true', 'dXdt_true', 't_sparse', ...
    'X_sparse', 'X_true_at_sparse', 'params', 'noise_level', ...
    'n_sparse_points', 'noise_fraction');

fprintf('  Loaded data: %d sparse points, %.1f%% noise\n', ...
    n_sparse_points, noise_fraction*100);
fprintf('  True solution: %d time points\n', length(t_true));

%% GP Fitting Parameters
fprintf('\n=== Phase 2: GP Implementation ===\n');

% Dense time grid for GP predictions (same as ground truth for comparison)
t_dense = t_true;

% GP options
gp_options = struct();
gp_options.kernel = 'squaredexponential';  % Squared exponential kernel
gp_options.standardize = false;  % Keep false for clean derivative math
gp_options.basis = 'constant';  % Constant basis function

fprintf('\nFitting GP models...\n');
fprintf('  Kernel: %s\n', gp_options.kernel);
fprintf('  Basis: %s\n', gp_options.basis);
fprintf('  Standardize: %s\n', mat2str(gp_options.standardize));

%% Fit GP Models
tic;
[gp_models, X_gp, dXdt_gp, gp_stats] = fit_gp_vdp(t_sparse, X_sparse, t_dense, gp_options);
gp_time = toc;

fprintf('  GP fitting completed in %.2f seconds\n', gp_time);

%% Display GP Quality Metrics
fprintf('\nGP Fit Quality:\n');
for i = 1:length(gp_stats.var_names)
    var_name = gp_stats.var_names{i};
    r2 = gp_stats.r2(i);
    L = gp_stats.hyperparams{i}.lengthscale;
    sigma_f = gp_stats.hyperparams{i}.signal_std;
    sigma_n = gp_stats.hyperparams{i}.noise_std;
    
    fprintf('  %s: R² = %.4f\n', var_name, r2);
    fprintf('    Lengthscale: %.4f, Signal std: %.4f, Noise std: %.4f\n', ...
        L, sigma_f, sigma_n);
    
    if r2 < 0.5
        warning('GP fit quality for %s is poor (R² = %.4f < 0.5)', var_name, r2);
    end
end

%% Compute Prediction Errors
fprintf('\nGP Prediction Errors (vs ground truth):\n');
for i = 1:length(gp_stats.var_names)
    var_name = gp_stats.var_names{i};
    
    % State prediction error
    state_error = X_gp(:, i) - X_true(:, i);
    rmse_state = sqrt(mean(state_error.^2));
    mae_state = mean(abs(state_error));
    
    % Derivative prediction error
    deriv_error = dXdt_gp(:, i) - dXdt_true(:, i);
    rmse_deriv = sqrt(mean(deriv_error.^2));
    mae_deriv = mean(abs(deriv_error));
    
    fprintf('  %s:\n', var_name);
    fprintf('    State: RMSE = %.4f, MAE = %.4f\n', rmse_state, mae_state);
    fprintf('    Derivative: RMSE = %.4f, MAE = %.4f\n', rmse_deriv, mae_deriv);
end

%% Visualization: GP Fit Overview
fprintf('\nCreating visualizations...\n');

figure('Color', 'w', 'Position', [100, 100, 1400, 900]);

% Plot 1: State x(t) - GP fit
subplot(2, 3, 1);
plot(t_true, X_true(:,1), '-', 'Color', [0, 0.4470, 0.7410], ...
    'LineWidth', 2, 'DisplayName', 'True x(t)'); hold on;
scatter(t_sparse, X_sparse(:,1), 100, [0.8500, 0.3250, 0.0980], ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Noisy samples');
plot(t_dense, X_gp(:,1), '--', 'Color', [0.9290, 0.6940, 0.1250], ...
    'LineWidth', 2, 'DisplayName', sprintf('GP fit (R²=%.3f)', gp_stats.r2(1)));
xlabel('Time'); ylabel('x(t)');
title('GP Fit: x(t)');
legend('Location', 'best'); grid on;

% Plot 2: State y(t) - GP fit
subplot(2, 3, 2);
plot(t_true, X_true(:,2), '-', 'Color', [0.8500, 0.3250, 0.0980], ...
    'LineWidth', 2, 'DisplayName', 'True y(t)'); hold on;
scatter(t_sparse, X_sparse(:,2), 100, [0, 0.4470, 0.7410], ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Noisy samples');
plot(t_dense, X_gp(:,2), '--', 'Color', [0.9290, 0.6940, 0.1250], ...
    'LineWidth', 2, 'DisplayName', sprintf('GP fit (R²=%.3f)', gp_stats.r2(2)));
xlabel('Time'); ylabel('y(t)');
title('GP Fit: y(t)');
legend('Location', 'best'); grid on;

% Plot 3: Phase space
subplot(2, 3, 3);
plot(X_true(:,1), X_true(:,2), '-', 'Color', [0, 0.4470, 0.7410], ...
    'LineWidth', 2, 'DisplayName', 'True trajectory'); hold on;
scatter(X_sparse(:,1), X_sparse(:,2), 100, [0.8500, 0.3250, 0.0980], ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Noisy samples');
plot(X_gp(:,1), X_gp(:,2), '--', 'Color', [0.9290, 0.6940, 0.1250], ...
    'LineWidth', 2, 'DisplayName', 'GP fit');
xlabel('x'); ylabel('y');
title('Phase Space');
legend('Location', 'best'); grid on; axis equal;

% Plot 4: Derivative dx/dt
subplot(2, 3, 4);
plot(t_true, dXdt_true(:,1), '-', 'Color', [0, 0.4470, 0.7410], ...
    'LineWidth', 2, 'DisplayName', 'True dx/dt'); hold on;
plot(t_dense, dXdt_gp(:,1), '--', 'Color', [0.9290, 0.6940, 0.1250], ...
    'LineWidth', 2, 'DisplayName', 'GP derivative');
xlabel('Time'); ylabel('dx/dt');
title('GP Derivative: dx/dt');
legend('Location', 'best'); grid on;

% Plot 5: Derivative dy/dt
subplot(2, 3, 5);
plot(t_true, dXdt_true(:,2), '-', 'Color', [0.8500, 0.3250, 0.0980], ...
    'LineWidth', 2, 'DisplayName', 'True dy/dt'); hold on;
plot(t_dense, dXdt_gp(:,2), '--', 'Color', [0.9290, 0.6940, 0.1250], ...
    'LineWidth', 2, 'DisplayName', 'GP derivative');
xlabel('Time'); ylabel('dy/dt');
title('GP Derivative: dy/dt');
legend('Location', 'best'); grid on;

% Plot 6: Prediction errors
subplot(2, 3, 6);
state_error_x = X_gp(:,1) - X_true(:,1);
state_error_y = X_gp(:,2) - X_true(:,2);
plot(t_dense, state_error_x, '-', 'Color', [0, 0.4470, 0.7410], ...
    'LineWidth', 1.5, 'DisplayName', 'x error'); hold on;
plot(t_dense, state_error_y, '-', 'Color', [0.8500, 0.3250, 0.0980], ...
    'LineWidth', 1.5, 'DisplayName', 'y error');
xlabel('Time'); ylabel('Prediction Error');
title('GP Prediction Errors');
legend('Location', 'best'); grid on;
yline(0, 'k--', 'LineWidth', 1);

sgtitle(sprintf('Phase 2: GP Implementation - %d sparse points, %.0f%% noise', ...
    n_sparse_points, noise_fraction*100), 'FontSize', 14, 'FontWeight', 'bold');

%% Save Phase 2 results
fprintf('\nSaving Phase 2 results...\n');
save('vdp_data.mat', 'gp_models', 'X_gp', 'dXdt_gp', 'gp_stats', ...
    'gp_time', 't_dense', '-append');

fprintf('Results saved to vdp_data.mat\n');
fprintf('\nPhase 2 complete! Ready for Phase 3 (SINDy-Only Implementation).\n');

