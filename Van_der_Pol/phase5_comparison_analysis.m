%% Phase 5: Comparison & Analysis for Van der Pol Oscillator
% This script performs comprehensive comparison and analysis of all methods:
% - GP-only
% - SINDy-only
% - GP+SINDy
%
% Requirements:
%   - Phases 1-4 must be run first

clear; close all; clc;

%% Load all data from previous phases
fprintf('=== Phase 5: Comparison & Analysis ===\n');
fprintf('\nLoading data from all previous phases...\n');

% Phase 1 data
if ~exist('vdp_data.mat', 'file')
    error('Phase 1 data not found. Please run phase1_data_generation.m first.');
end

load('vdp_data.mat', 't_true', 'X_true', 'dXdt_true', 't_sparse', ...
    'X_sparse', 'params', 'n_sparse_points', 'noise_fraction');

% Phase 2 data (GP-only)
try
    load('vdp_data.mat', 'gp_models', 'X_gp', 'dXdt_gp', 't_dense', 'gp_stats');
    phase2_available = true;
    fprintf('  Phase 2 (GP-only): ✓\n');
catch
    phase2_available = false;
    fprintf('  Phase 2 (GP-only): ✗\n');
end

% Phase 3 data (SINDy-only)
try
    load('vdp_data.mat', 'Xi', 'X_sim', 'X_sim_interp', 't_sim', ...
        'dXdt_fd', 'rmse_traj', 'mae_traj', 'rmse_fd');
    phase3_available = true;
    fprintf('  Phase 3 (SINDy-only): ✓\n');
catch
    phase3_available = false;
    fprintf('  Phase 3 (SINDy-only): ✗\n');
end

% Phase 4 data (GP+SINDy)
try
    load('vdp_data.mat', 'Xi_aug', 'X_sim_aug', 'X_sim_aug_interp', ...
        't_sim_aug', 'dXdt_aug', 'rmse_traj_aug', 'mae_traj_aug');
    phase4_available = true;
    fprintf('  Phase 4 (GP+SINDy): ✓\n');
catch
    phase4_available = false;
    fprintf('  Phase 4 (GP+SINDy): ✗\n');
end

if ~phase2_available && ~phase3_available && ~phase4_available
    error('No phase data available. Please run at least one of phases 2-4.');
end

%% Compute true coefficients for comparison
fprintf('\nComputing true Van der Pol coefficients...\n');
[Xi_true, library_names] = compute_true_vdp_coefficients(params.mu, 3);
fprintf('  True model: dx/dt = y, dy/dt = %.2f(1-x²)y - x\n', params.mu);

%% Compute metrics for each method
fprintf('\nComputing comprehensive metrics...\n');

all_metrics = struct();

% GP-only metrics
X_gp_interp = [];
if phase2_available
    fprintf('  Computing GP-only metrics...\n');
    % Interpolate GP predictions to true time grid
    X_gp_interp = interp1(t_dense, X_gp, t_true, 'pchip');
    dXdt_gp_interp = interp1(t_dense, dXdt_gp, t_true, 'pchip');
    
    all_metrics.gp_only = compute_comprehensive_metrics(...
        X_true, dXdt_true, X_gp_interp, dXdt_gp_interp, ...
        [], [], library_names, 'GP-only');
end

% SINDy-only metrics
if phase3_available
    fprintf('  Computing SINDy-only metrics...\n');
    % Interpolate SINDy simulation to true time grid
    if ~exist('X_sim_interp', 'var')
        X_sim_interp = interp1(t_sim, X_sim, t_true, 'pchip');
    end
    
    % Compute derivatives from SINDy model at sparse points
    % (We'll use the finite difference derivatives from Phase 3)
    dXdt_sindy_interp = [];
    if exist('dXdt_fd', 'var')
        % Interpolate finite difference derivatives
        dXdt_sindy_interp = interp1(t_sparse, dXdt_fd, t_true, 'pchip');
    end
    
    all_metrics.sindy_only = compute_comprehensive_metrics(...
        X_true, dXdt_true, X_sim_interp, dXdt_sindy_interp, ...
        Xi, Xi_true, library_names, 'SINDy-only');
end

% GP+SINDy metrics
if phase4_available
    fprintf('  Computing GP+SINDy metrics...\n');
    % Interpolate GP+SINDy simulation to true time grid
    if ~exist('X_sim_aug_interp', 'var')
        X_sim_aug_interp = interp1(t_sim_aug, X_sim_aug, t_true, 'pchip');
    end
    
    % Use augmented derivatives if available
    dXdt_aug_interp = [];
    if exist('dXdt_aug', 'var')
        % Interpolate augmented derivatives
        load('vdp_data.mat', 't_aug');
        dXdt_aug_interp = interp1(t_aug, dXdt_aug, t_true, 'pchip');
    end
    
    all_metrics.gp_sindy = compute_comprehensive_metrics(...
        X_true, dXdt_true, X_sim_aug_interp, dXdt_aug_interp, ...
        Xi_aug, Xi_true, library_names, 'GP+SINDy');
end

%% Display metrics summary
fprintf('\n=== Metrics Summary ===\n');
fprintf('\nTrajectory Metrics:\n');
fprintf('%-15s | %10s | %10s | %10s\n', 'Method', 'RMSE', 'MAE', 'R²');
fprintf('%s\n', repmat('-', 1, 50));

if phase2_available
    fprintf('%-15s | %10.4f | %10.4f | %10.4f\n', ...
        'GP-only', all_metrics.gp_only.rmse_traj, ...
        all_metrics.gp_only.mae_traj, all_metrics.gp_only.r2_traj);
end

if phase3_available
    fprintf('%-15s | %10.4f | %10.4f | %10.4f\n', ...
        'SINDy-only', all_metrics.sindy_only.rmse_traj, ...
        all_metrics.sindy_only.mae_traj, all_metrics.sindy_only.r2_traj);
end

if phase4_available
    fprintf('%-15s | %10.4f | %10.4f | %10.4f\n', ...
        'GP+SINDy', all_metrics.gp_sindy.rmse_traj, ...
        all_metrics.gp_sindy.mae_traj, all_metrics.gp_sindy.r2_traj);
end

fprintf('\nDerivative Metrics:\n');
fprintf('%-15s | %10s | %10s | %10s\n', 'Method', 'RMSE', 'MAE', 'Corr');
fprintf('%s\n', repmat('-', 1, 50));

if phase2_available && ~isnan(all_metrics.gp_only.rmse_deriv)
    fprintf('%-15s | %10.4f | %10.4f | %10.4f\n', ...
        'GP-only', all_metrics.gp_only.rmse_deriv, ...
        all_metrics.gp_only.mae_deriv, all_metrics.gp_only.corr_deriv);
end

if phase3_available && ~isnan(all_metrics.sindy_only.rmse_deriv)
    fprintf('%-15s | %10.4f | %10.4f | %10.4f\n', ...
        'SINDy-only', all_metrics.sindy_only.rmse_deriv, ...
        all_metrics.sindy_only.mae_deriv, all_metrics.sindy_only.corr_deriv);
end

if phase4_available && ~isnan(all_metrics.gp_sindy.rmse_deriv)
    fprintf('%-15s | %10.4f | %10.4f | %10.4f\n', ...
        'GP+SINDy', all_metrics.gp_sindy.rmse_deriv, ...
        all_metrics.gp_sindy.mae_deriv, all_metrics.gp_sindy.corr_deriv);
end

fprintf('\nCoefficient Metrics:\n');
fprintf('%-15s | %10s | %10s | %10s | %10s\n', 'Method', 'RMSE', 'MAE', 'Sparsity', 'ID Rate');
fprintf('%s\n', repmat('-', 1, 60));

if phase3_available && ~isnan(all_metrics.sindy_only.rmse_coeff)
    fprintf('%-15s | %10.4f | %10.4f | %10d | %10.2f%%\n', ...
        'SINDy-only', all_metrics.sindy_only.rmse_coeff, ...
        all_metrics.sindy_only.mae_coeff, all_metrics.sindy_only.sparsity_pred, ...
        all_metrics.sindy_only.coeff_identification_rate * 100);
end

if phase4_available && ~isnan(all_metrics.gp_sindy.rmse_coeff)
    fprintf('%-15s | %10.4f | %10.4f | %10d | %10.2f%%\n', ...
        'GP+SINDy', all_metrics.gp_sindy.rmse_coeff, ...
        all_metrics.gp_sindy.mae_coeff, all_metrics.gp_sindy.sparsity_pred, ...
        all_metrics.gp_sindy.coeff_identification_rate * 100);
end

%% Improvement analysis
if phase3_available && phase4_available
    fprintf('\n=== Improvement Analysis (GP+SINDy vs SINDy-only) ===\n');
    
    % Trajectory improvement
    traj_rmse_improvement = (1 - all_metrics.gp_sindy.rmse_traj / all_metrics.sindy_only.rmse_traj) * 100;
    traj_mae_improvement = (1 - all_metrics.gp_sindy.mae_traj / all_metrics.sindy_only.mae_traj) * 100;
    traj_r2_improvement = (all_metrics.gp_sindy.r2_traj - all_metrics.sindy_only.r2_traj) * 100;
    
    fprintf('Trajectory RMSE: %.1f%% improvement\n', traj_rmse_improvement);
    fprintf('Trajectory MAE:  %.1f%% improvement\n', traj_mae_improvement);
    fprintf('Trajectory R²:   %.1f%% improvement\n', traj_r2_improvement);
    
    % Coefficient improvement
    if ~isnan(all_metrics.sindy_only.rmse_coeff) && ~isnan(all_metrics.gp_sindy.rmse_coeff)
        coeff_rmse_improvement = (1 - all_metrics.gp_sindy.rmse_coeff / all_metrics.sindy_only.rmse_coeff) * 100;
        coeff_mae_improvement = (1 - all_metrics.gp_sindy.mae_coeff / all_metrics.sindy_only.mae_coeff) * 100;
        
        fprintf('\nCoefficient RMSE: %.1f%% improvement\n', coeff_rmse_improvement);
        fprintf('Coefficient MAE:  %.1f%% improvement\n', coeff_mae_improvement);
    end
end

%% Visualization: Comprehensive comparison
fprintf('\nCreating comprehensive comparison visualizations...\n');

figure('Color', 'w', 'Position', [100, 100, 1800, 1200]);

% Plot 1: Trajectory comparison - x(t)
subplot(3, 3, 1);
plot(t_true, X_true(:,1), '-', 'Color', [0, 0.4470, 0.7410], ...
    'LineWidth', 2.5, 'DisplayName', 'True'); hold on;
scatter(t_sparse, X_sparse(:,1), 100, [0.8500, 0.3250, 0.0980], ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', 'Data');
if phase2_available
    plot(t_dense, X_gp(:,1), '--', 'Color', [0.4660, 0.6740, 0.1880], ...
        'LineWidth', 2, 'DisplayName', 'GP-only');
end
if phase3_available
    plot(t_sim, X_sim(:,1), ':', 'Color', [0.4940, 0.1840, 0.5560], ...
        'LineWidth', 2, 'DisplayName', 'SINDy-only');
end
if phase4_available
    plot(t_sim_aug, X_sim_aug(:,1), '-.', 'Color', [0.9290, 0.6940, 0.1250], ...
        'LineWidth', 2, 'DisplayName', 'GP+SINDy');
end
xlabel('Time'); ylabel('x(t)');
title('Trajectory Comparison: x(t)');
legend('Location', 'best', 'FontSize', 9); grid on;

% Plot 2: Trajectory comparison - y(t)
subplot(3, 3, 2);
plot(t_true, X_true(:,2), '-', 'Color', [0.8500, 0.3250, 0.0980], ...
    'LineWidth', 2.5, 'DisplayName', 'True'); hold on;
scatter(t_sparse, X_sparse(:,2), 100, [0, 0.4470, 0.7410], ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', 'Data');
if phase2_available
    plot(t_dense, X_gp(:,2), '--', 'Color', [0.4660, 0.6740, 0.1880], ...
        'LineWidth', 2, 'DisplayName', 'GP-only');
end
if phase3_available
    plot(t_sim, X_sim(:,2), ':', 'Color', [0.4940, 0.1840, 0.5560], ...
        'LineWidth', 2, 'DisplayName', 'SINDy-only');
end
if phase4_available
    plot(t_sim_aug, X_sim_aug(:,2), '-.', 'Color', [0.9290, 0.6940, 0.1250], ...
        'LineWidth', 2, 'DisplayName', 'GP+SINDy');
end
xlabel('Time'); ylabel('y(t)');
title('Trajectory Comparison: y(t)');
legend('Location', 'best', 'FontSize', 9); grid on;

% Plot 3: Phase space
subplot(3, 3, 3);
plot(X_true(:,1), X_true(:,2), '-', 'Color', [0, 0.4470, 0.7410], ...
    'LineWidth', 2.5, 'DisplayName', 'True'); hold on;
scatter(X_sparse(:,1), X_sparse(:,2), 100, [0.8500, 0.3250, 0.0980], ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', 'Data');
if phase2_available
    plot(X_gp(:,1), X_gp(:,2), '--', 'Color', [0.4660, 0.6740, 0.1880], ...
        'LineWidth', 2, 'DisplayName', 'GP-only');
end
if phase3_available
    plot(X_sim(:,1), X_sim(:,2), ':', 'Color', [0.4940, 0.1840, 0.5560], ...
        'LineWidth', 2, 'DisplayName', 'SINDy-only');
end
if phase4_available
    plot(X_sim_aug(:,1), X_sim_aug(:,2), '-.', 'Color', [0.9290, 0.6940, 0.1250], ...
        'LineWidth', 2, 'DisplayName', 'GP+SINDy');
end
xlabel('x'); ylabel('y');
title('Phase Space Comparison');
legend('Location', 'best', 'FontSize', 9); grid on; axis equal;

% Plot 4: Trajectory RMSE comparison
subplot(3, 3, 4);
methods = {};
rmse_values = [];
if phase2_available
    methods{end+1} = 'GP-only';
    rmse_values(end+1) = all_metrics.gp_only.rmse_traj;
end
if phase3_available
    methods{end+1} = 'SINDy-only';
    rmse_values(end+1) = all_metrics.sindy_only.rmse_traj;
end
if phase4_available
    methods{end+1} = 'GP+SINDy';
    rmse_values(end+1) = all_metrics.gp_sindy.rmse_traj;
end
bar(rmse_values);
set(gca, 'XTickLabel', methods);
ylabel('RMSE');
title('Trajectory RMSE Comparison');
grid on;

% Plot 5: Trajectory R² comparison
subplot(3, 3, 5);
r2_values = [];
if phase2_available
    r2_values(end+1) = all_metrics.gp_only.r2_traj;
end
if phase3_available
    r2_values(end+1) = all_metrics.sindy_only.r2_traj;
end
if phase4_available
    r2_values(end+1) = all_metrics.gp_sindy.r2_traj;
end
bar(r2_values);
set(gca, 'XTickLabel', methods);
ylabel('R²');
title('Trajectory R² Comparison');
grid on;
ylim([0, 1.1]);

% Plot 6: Coefficient comparison (if available)
subplot(3, 3, 6);
if phase3_available && phase4_available
    % Compare coefficients for dx/dt
    bar([Xi(:,1), Xi_aug(:,1), Xi_true(:,1)]);
    set(gca, 'XTickLabel', library_names);
    xtickangle(45);
    ylabel('Coefficient Value');
    title('Coefficients: dx/dt');
    legend('SINDy-only', 'GP+SINDy', 'True', 'Location', 'best', 'FontSize', 8);
    grid on;
else
    text(0.5, 0.5, 'Run Phases 3-4 for coefficient comparison', ...
        'HorizontalAlignment', 'center', 'FontSize', 12);
    axis off;
end

% Plot 7: Coefficient comparison for dy/dt
subplot(3, 3, 7);
if phase3_available && phase4_available
    bar([Xi(:,2), Xi_aug(:,2), Xi_true(:,2)]);
    set(gca, 'XTickLabel', library_names);
    xtickangle(45);
    ylabel('Coefficient Value');
    title('Coefficients: dy/dt');
    legend('SINDy-only', 'GP+SINDy', 'True', 'Location', 'best', 'FontSize', 8);
    grid on;
else
    text(0.5, 0.5, 'Run Phases 3-4 for coefficient comparison', ...
        'HorizontalAlignment', 'center', 'FontSize', 12);
    axis off;
end

% Plot 8: Error evolution over time
subplot(3, 3, 8);
if phase2_available && ~isempty(X_gp_interp)
    traj_error_gp = X_gp_interp - X_true;
    plot(t_true, sqrt(sum(traj_error_gp.^2, 2)), '--', ...
        'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1.5, 'DisplayName', 'GP-only'); hold on;
end
if phase3_available
    traj_error_sindy = X_sim_interp - X_true;
    plot(t_true, sqrt(sum(traj_error_sindy.^2, 2)), ':', ...
        'Color', [0.4940, 0.1840, 0.5560], 'LineWidth', 1.5, 'DisplayName', 'SINDy-only');
end
if phase4_available
    traj_error_aug = X_sim_aug_interp - X_true;
    plot(t_true, sqrt(sum(traj_error_aug.^2, 2)), '-.', ...
        'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1.5, 'DisplayName', 'GP+SINDy');
end
xlabel('Time'); ylabel('Error Magnitude');
title('Error Evolution Over Time');
legend('Location', 'best', 'FontSize', 9); grid on;

% Plot 9: Summary metrics table (text)
subplot(3, 3, 9);
axis off;
text_str = sprintf('Summary Metrics\n\n');
if phase2_available
    text_str = [text_str sprintf('GP-only:\n  RMSE: %.4f\n  R²: %.4f\n\n', ...
        all_metrics.gp_only.rmse_traj, all_metrics.gp_only.r2_traj)];
end
if phase3_available
    text_str = [text_str sprintf('SINDy-only:\n  RMSE: %.4f\n  R²: %.4f\n  Sparsity: %d\n\n', ...
        all_metrics.sindy_only.rmse_traj, all_metrics.sindy_only.r2_traj, ...
        all_metrics.sindy_only.sparsity_pred)];
end
if phase4_available
    text_str = [text_str sprintf('GP+SINDy:\n  RMSE: %.4f\n  R²: %.4f\n  Sparsity: %d\n', ...
        all_metrics.gp_sindy.rmse_traj, all_metrics.gp_sindy.r2_traj, ...
        all_metrics.gp_sindy.sparsity_pred)];
end
text(0.1, 0.5, text_str, 'FontSize', 10, ...
    'VerticalAlignment', 'middle');

sgtitle(sprintf('Phase 5: Comprehensive Comparison - %d sparse points, %.0f%% noise', ...
    n_sparse_points, noise_fraction*100), 'FontSize', 14, 'FontWeight', 'bold');

%% Save Phase 5 results
fprintf('\nSaving Phase 5 results...\n');
save('vdp_data.mat', 'all_metrics', 'Xi_true', '-append');

fprintf('Results saved to vdp_data.mat\n');
fprintf('\nPhase 5 complete! All methods have been compared and analyzed.\n');

