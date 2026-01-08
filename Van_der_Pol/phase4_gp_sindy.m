%% Phase 4: GP+SINDy Combined Implementation for Van der Pol Oscillator
% This script combines GP posterior sampling with SINDy to improve
% model discovery on sparse noisy data
%
% Requirements:
%   - Phase 1 must be run first (generates vdp_data.mat)
%   - Phase 2 must be run first (generates GP models)
%   - Statistics and Machine Learning Toolbox (for fitrgp)

clear; close all; clc;

%% Check if required data exists
if ~exist('vdp_data.mat', 'file')
    error('Phase 1 data not found. Please run phase1_data_generation.m first.');
end

%% Load data from previous phases
fprintf('Loading data from previous phases...\n');
load('vdp_data.mat', 't_true', 'X_true', 'dXdt_true', 't_sparse', ...
    'X_sparse', 'X_true_at_sparse', 'params', 'noise_level', ...
    'n_sparse_points', 'noise_fraction');

% Check for Phase 2 data
if ~exist('gp_models', 'var')
    load('vdp_data.mat', 'gp_models', 'X_gp', 'dXdt_gp', 't_dense', 'gp_stats');
    if ~exist('gp_models', 'var')
        error('Phase 2 data not found. Please run phase2_gp_implementation.m first.');
    end
end

% Load GP predictions if not already loaded
if ~exist('X_gp', 'var')
    load('vdp_data.mat', 'X_gp', 't_dense');
end

fprintf('  Loaded data: %d sparse points, %.1f%% noise\n', ...
    n_sparse_points, noise_fraction*100);
fprintf('  GP models available: %d state variables\n', length(gp_models));

%% Phase 4: GP+SINDy Implementation
fprintf('\n=== Phase 4: GP+SINDy Combined Implementation ===\n');

%% Step 1: Sample from GP Posterior and Augment Data
fprintf('\nStep 1: Sampling from GP posterior and augmenting data...\n');

num_gp_samples = 1;  % Number of GP sample curves (can increase for ensemble)
augment_options = struct();
augment_options.combine_method = 'append';  % Append GP samples to original data

fprintf('  Number of GP samples: %d\n', num_gp_samples);
fprintf('  Combine method: %s\n', augment_options.combine_method);

tic;
[t_aug, X_aug, dXdt_aug, aug_stats] = augment_data_with_gp_samples(...
    t_sparse, X_sparse, gp_models, t_dense, num_gp_samples, augment_options);
augmentation_time = toc;

fprintf('  Augmentation completed in %.4f seconds\n', augmentation_time);
fprintf('  Original points: %d\n', aug_stats.n_original);
fprintf('  GP sample points: %d\n', aug_stats.n_gp_samples);
fprintf('  Total augmented points: %d\n', aug_stats.n_augmented);
fprintf('  Augmentation factor: %.1fx\n', aug_stats.n_augmented / aug_stats.n_original);

%% Step 2: Build Library Matrix for Augmented Data
fprintf('\nStep 2: Building library matrix for augmented data...\n');
polyOrder = 3;  % Same as Phase 3
[Theta_aug, library_names] = build_library_vdp(X_aug, polyOrder);

fprintf('  Library order: %d\n', polyOrder);
fprintf('  Number of terms: %d\n', length(library_names));
fprintf('  Data points: %d\n', size(Theta_aug, 1));

%% Step 3: Run SINDy on Augmented Dataset (ADMM for LASSO)
fprintf('\nStep 3: Running SINDy on augmented dataset (ADMM for LASSO)...\n');

% Set sparsity threshold (data-adaptive)
lambda_base = 0.01;
lambda = lambda_base * std(dXdt_aug(:));

% ADMM options (matching Hsin et al. 2025 approach)
admm_options = struct();
admm_options.rho = 1.0;  % ADMM penalty parameter
admm_options.max_iterations = 1000;
admm_options.abs_tol = 1e-4;
admm_options.rel_tol = 1e-2;
admm_options.verbose = false;

fprintf('  Sparsity parameter: λ = %.6f (%.4f * std(dXdt))\n', lambda, lambda_base);
fprintf('  ADMM penalty: ρ = %.2f\n', admm_options.rho);
fprintf('  Max iterations: %d\n', admm_options.max_iterations);

tic;
[Xi_aug, sindy_stats_aug] = run_sindy_admm(Theta_aug, dXdt_aug, lambda, admm_options);
sindy_time_aug = toc;

fprintf('  SINDy (ADMM) completed in %.4f seconds\n', sindy_time_aug);
fprintf('  Converged: %s (iterations: %d)\n', ...
    mat2str(sindy_stats_aug.converged), sindy_stats_aug.iterations);
fprintf('  Final sparsity: %.1f%% (%.0f/%d nonzero terms)\n', ...
    sindy_stats_aug.final_sparsity * 100, sindy_stats_aug.num_nonzero, sindy_stats_aug.num_terms);
fprintf('  Prediction RMSE: %.4f\n', sindy_stats_aug.rmse);

% Display discovered coefficients
fprintf('\n  Discovered coefficients (GP+SINDy):\n');
fprintf('    Term        | dx/dt      | dy/dt\n');
fprintf('    ------------|------------|------------\n');
for i = 1:length(library_names)
    fprintf('    %-10s | %10.4f | %10.4f\n', library_names{i}, Xi_aug(i,1), Xi_aug(i,2));
end

%% Step 4: Forward Integration
fprintf('\nStep 4: Forward integrating GP+SINDy model...\n');

t_span = [min(t_true), max(t_true)];
x0 = params.x0;

fprintf('  Time span: [%.2f, %.2f]\n', t_span(1), t_span(2));
fprintf('  Initial conditions: x0 = [%.2f; %.2f]\n', x0(1), x0(2));

tic;
[t_sim_aug, X_sim_aug] = forward_integrate_sindy(Xi_aug, library_names, t_span, x0);
integration_time_aug = toc;

fprintf('  Integration completed in %.4f seconds\n', integration_time_aug);
fprintf('  Simulated %d time points\n', length(t_sim_aug));

% Interpolate to true time grid for comparison
X_sim_aug_interp = interp1(t_sim_aug, X_sim_aug, t_true, 'pchip');

% Compute trajectory errors
traj_error_aug = X_sim_aug_interp - X_true;
rmse_traj_aug = sqrt(mean(traj_error_aug(:).^2));
mae_traj_aug = mean(abs(traj_error_aug(:)));

fprintf('  Trajectory error (vs ground truth): RMSE = %.4f, MAE = %.4f\n', ...
    rmse_traj_aug, mae_traj_aug);

%% Step 5: Comparison with Previous Phases
fprintf('\nStep 5: Comparing with previous phases...\n');

% Load Phase 3 results if available
try
    load('vdp_data.mat', 'Xi', 'X_sim', 'X_sim_interp', 'rmse_traj', 'mae_traj', 't_sim');
    phase3_available = true;
catch
    phase3_available = false;
end

if phase3_available
    fprintf('\n  Trajectory Errors (RMSE):\n');
    fprintf('    SINDy-only:     %.4f\n', rmse_traj);
    fprintf('    GP+SINDy:       %.4f\n', rmse_traj_aug);
    fprintf('    Improvement:    %.1f%%\n', (1 - rmse_traj_aug/rmse_traj) * 100);
    
    fprintf('\n  Trajectory Errors (MAE):\n');
    fprintf('    SINDy-only:     %.4f\n', mae_traj);
    fprintf('    GP+SINDy:       %.4f\n', mae_traj_aug);
    fprintf('    Improvement:    %.1f%%\n', (1 - mae_traj_aug/mae_traj) * 100);
else
    fprintf('  Phase 3 results not found for comparison\n');
end

%% Visualization
fprintf('\nCreating visualizations...\n');

figure('Color', 'w', 'Position', [100, 100, 1600, 1000]);

% Plot 1: Augmented data overview
subplot(3, 3, 1);
plot(t_true, X_true(:,1), '-', 'Color', [0, 0.4470, 0.7410], ...
    'LineWidth', 2, 'DisplayName', 'True x(t)'); hold on;
scatter(t_sparse, X_sparse(:,1), 100, [0.8500, 0.3250, 0.0980], ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Sparse data');
plot(t_aug, X_aug(:,1), '-', 'Color', [0.9290, 0.6940, 0.1250], ...
    'LineWidth', 0.5, 'DisplayName', 'Augmented data');
xlabel('Time'); ylabel('x(t)');
title('Augmented Data: x(t)');
legend('Location', 'best'); grid on;

subplot(3, 3, 2);
plot(t_true, X_true(:,2), '-', 'Color', [0.8500, 0.3250, 0.0980], ...
    'LineWidth', 2, 'DisplayName', 'True y(t)'); hold on;
scatter(t_sparse, X_sparse(:,2), 100, [0, 0.4470, 0.7410], ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Sparse data');
plot(t_aug, X_aug(:,2), '-', 'Color', [0.9290, 0.6940, 0.1250], ...
    'LineWidth', 0.5, 'DisplayName', 'Augmented data');
xlabel('Time'); ylabel('y(t)');
title('Augmented Data: y(t)');
legend('Location', 'best'); grid on;

% Plot 2: Coefficients comparison
subplot(3, 3, 3);
if exist('Xi', 'var')
    % Compare SINDy-only vs GP+SINDy coefficients
    bar([Xi(:,1), Xi_aug(:,1)]);
    set(gca, 'XTickLabel', library_names);
    xtickangle(45);
    ylabel('Coefficient Value');
    title('Coefficients: dx/dt');
    legend('SINDy-only', 'GP+SINDy', 'Location', 'best');
    grid on;
else
    bar(Xi_aug(:,1));
    set(gca, 'XTickLabel', library_names);
    xtickangle(45);
    ylabel('Coefficient Value');
    title('GP+SINDy Coefficients: dx/dt');
    grid on;
end

% Plot 3: State trajectories
subplot(3, 3, 4);
plot(t_true, X_true(:,1), '-', 'Color', [0, 0.4470, 0.7410], ...
    'LineWidth', 2, 'DisplayName', 'True x(t)'); hold on;
scatter(t_sparse, X_sparse(:,1), 100, [0.8500, 0.3250, 0.0980], ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Sparse data');
plot(t_sim_aug, X_sim_aug(:,1), '--', 'Color', [0.9290, 0.6940, 0.1250], ...
    'LineWidth', 2, 'DisplayName', 'GP+SINDy model');
if phase3_available
    plot(t_sim, X_sim(:,1), ':', 'Color', [0.4940, 0.1840, 0.5560], ...
        'LineWidth', 2, 'DisplayName', 'SINDy-only');
end
xlabel('Time'); ylabel('x(t)');
title('GP+SINDy Model: x(t)');
legend('Location', 'best'); grid on;

subplot(3, 3, 5);
plot(t_true, X_true(:,2), '-', 'Color', [0.8500, 0.3250, 0.0980], ...
    'LineWidth', 2, 'DisplayName', 'True y(t)'); hold on;
scatter(t_sparse, X_sparse(:,2), 100, [0, 0.4470, 0.7410], ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Sparse data');
plot(t_sim_aug, X_sim_aug(:,2), '--', 'Color', [0.9290, 0.6940, 0.1250], ...
    'LineWidth', 2, 'DisplayName', 'GP+SINDy model');
if phase3_available
    plot(t_sim, X_sim(:,2), ':', 'Color', [0.4940, 0.1840, 0.5560], ...
        'LineWidth', 2, 'DisplayName', 'SINDy-only');
end
xlabel('Time'); ylabel('y(t)');
title('GP+SINDy Model: y(t)');
legend('Location', 'best'); grid on;

% Plot 4: Phase space
subplot(3, 3, 6);
plot(X_true(:,1), X_true(:,2), '-', 'Color', [0, 0.4470, 0.7410], ...
    'LineWidth', 2, 'DisplayName', 'True trajectory'); hold on;
scatter(X_sparse(:,1), X_sparse(:,2), 100, [0.8500, 0.3250, 0.0980], ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Sparse data');
plot(X_sim_aug(:,1), X_sim_aug(:,2), '--', 'Color', [0.9290, 0.6940, 0.1250], ...
    'LineWidth', 2, 'DisplayName', 'GP+SINDy model');
if phase3_available
    plot(X_sim(:,1), X_sim(:,2), ':', 'Color', [0.4940, 0.1840, 0.5560], ...
        'LineWidth', 2, 'DisplayName', 'SINDy-only');
end
xlabel('x'); ylabel('y');
title('Phase Space Comparison');
legend('Location', 'best'); grid on; axis equal;

% Plot 5: Prediction errors
subplot(3, 3, 7);
plot(t_true, traj_error_aug(:,1), '-', 'Color', [0, 0.4470, 0.7410], ...
    'LineWidth', 1.5, 'DisplayName', 'x error'); hold on;
plot(t_true, traj_error_aug(:,2), '-', 'Color', [0.8500, 0.3250, 0.0980], ...
    'LineWidth', 1.5, 'DisplayName', 'y error');
if phase3_available
    if ~exist('X_sim_interp', 'var')
        X_sim_interp = interp1(t_sim, X_sim, t_true, 'pchip');
    end
    traj_error_sindy = X_sim_interp - X_true;
    plot(t_true, traj_error_sindy(:,1), ':', 'Color', [0.4940, 0.1840, 0.5560], ...
        'LineWidth', 1.5, 'DisplayName', 'SINDy-only x error');
    plot(t_true, traj_error_sindy(:,2), ':', 'Color', [0.4940, 0.1840, 0.5560], ...
        'LineWidth', 1.5, 'DisplayName', 'SINDy-only y error');
end
xlabel('Time'); ylabel('Prediction Error');
title('Trajectory Errors');
legend('Location', 'best'); grid on;
yline(0, 'k--', 'LineWidth', 1);

% Plot 6: Data augmentation comparison
subplot(3, 3, 8);
bar([aug_stats.n_original, aug_stats.n_gp_samples, aug_stats.n_augmented]);
set(gca, 'XTickLabel', {'Original', 'GP Samples', 'Augmented'});
ylabel('Number of Points');
title('Data Augmentation');
grid on;

% Plot 7: Method comparison (if Phase 3 available)
subplot(3, 3, 9);
if phase3_available
    methods = {'SINDy-only', 'GP+SINDy'};
    rmse_values = [rmse_traj, rmse_traj_aug];
    bar(rmse_values);
    set(gca, 'XTickLabel', methods);
    ylabel('RMSE');
    title('Method Comparison: Trajectory RMSE');
    grid on;
else
    text(0.5, 0.5, 'Run Phase 3 for comparison', ...
        'HorizontalAlignment', 'center', 'FontSize', 12);
    axis off;
end

sgtitle(sprintf('Phase 4: GP+SINDy - %d sparse points, %.0f%% noise', ...
    n_sparse_points, noise_fraction*100), 'FontSize', 14, 'FontWeight', 'bold');

%% Create comprehensive comparison figure
fprintf('\nCreating comprehensive comparison figure...\n');

figure('Color', 'w', 'Position', [100, 100, 1600, 600]);

% Plot x(t) comparison
subplot(1, 2, 1);
plot(t_true, X_true(:,1), '-', 'Color', [0, 0.4470, 0.7410], ...
    'LineWidth', 2.5, 'DisplayName', 'True x(t)'); hold on;
scatter(t_sparse, X_sparse(:,1), 120, [0.8500, 0.3250, 0.0980], ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Sparse data');
plot(t_dense, X_gp(:,1), '--', 'Color', [0.4660, 0.6740, 0.1880], ...
    'LineWidth', 2, 'DisplayName', 'GP-only');
if phase3_available
    plot(t_sim, X_sim(:,1), ':', 'Color', [0.4940, 0.1840, 0.5560], ...
        'LineWidth', 2, 'DisplayName', 'SINDy-only');
end
plot(t_sim_aug, X_sim_aug(:,1), '-.', 'Color', [0.9290, 0.6940, 0.1250], ...
    'LineWidth', 2, 'DisplayName', 'GP+SINDy');
xlabel('Time'); ylabel('x(t)');
title('Comprehensive Comparison: x(t)');
legend('Location', 'best', 'FontSize', 10); grid on;

% Plot y(t) comparison
subplot(1, 2, 2);
plot(t_true, X_true(:,2), '-', 'Color', [0.8500, 0.3250, 0.0980], ...
    'LineWidth', 2.5, 'DisplayName', 'True y(t)'); hold on;
scatter(t_sparse, X_sparse(:,2), 120, [0, 0.4470, 0.7410], ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Sparse data');
plot(t_dense, X_gp(:,2), '--', 'Color', [0.4660, 0.6740, 0.1880], ...
    'LineWidth', 2, 'DisplayName', 'GP-only');
if phase3_available
    plot(t_sim, X_sim(:,2), ':', 'Color', [0.4940, 0.1840, 0.5560], ...
        'LineWidth', 2, 'DisplayName', 'SINDy-only');
end
plot(t_sim_aug, X_sim_aug(:,2), '-.', 'Color', [0.9290, 0.6940, 0.1250], ...
    'LineWidth', 2, 'DisplayName', 'GP+SINDy');
xlabel('Time'); ylabel('y(t)');
title('Comprehensive Comparison: y(t)');
legend('Location', 'best', 'FontSize', 10); grid on;

sgtitle(sprintf('All Methods Comparison - %d sparse points, %.0f%% noise', ...
    n_sparse_points, noise_fraction*100), 'FontSize', 14, 'FontWeight', 'bold');

%% Save Phase 4 results
fprintf('\nSaving Phase 4 results...\n');
save('vdp_data.mat', 't_aug', 'X_aug', 'dXdt_aug', 'aug_stats', ...
    'Theta_aug', 'Xi_aug', 'sindy_stats_aug', 't_sim_aug', 'X_sim_aug', ...
    'X_sim_aug_interp', 'rmse_traj_aug', 'mae_traj_aug', ...
    'augmentation_time', 'sindy_time_aug', 'integration_time_aug', '-append');

fprintf('Results saved to vdp_data.mat\n');
fprintf('\nPhase 4 complete! Ready for Phase 5 (Comparison & Analysis).\n');

