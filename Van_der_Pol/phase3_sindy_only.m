%% Phase 3: SINDy-Only Implementation for Van der Pol Oscillator
% This script runs SINDy using finite difference derivatives from sparse data,
% discovers the model, and forward integrates it
%
% Requirements:
%   - Phase 1 must be run first (generates vdp_data.mat)

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

%% Phase 3: SINDy-Only Implementation
fprintf('\n=== Phase 3: SINDy-Only Implementation ===\n');

%% Step 1: Compute Finite Difference Derivatives
fprintf('\nStep 1: Computing finite difference derivatives...\n');
tic;
dXdt_fd = compute_finite_differences(t_sparse, X_sparse, 'central');
fd_time = toc;

fprintf('  Computed derivatives in %.4f seconds\n', fd_time);
fprintf('  Derivative range: dx/dt ∈ [%.4f, %.4f], dy/dt ∈ [%.4f, %.4f]\n', ...
    min(dXdt_fd(:,1)), max(dXdt_fd(:,1)), min(dXdt_fd(:,2)), max(dXdt_fd(:,2)));

% Compare to true derivatives at sparse points
dXdt_true_at_sparse = interp1(t_true, dXdt_true, t_sparse, 'pchip');
fd_error = dXdt_fd - dXdt_true_at_sparse;
rmse_fd = sqrt(mean(fd_error(:).^2));
fprintf('  FD error (vs true): RMSE = %.4f\n', rmse_fd);

%% Step 2: Build Library Matrix
fprintf('\nStep 2: Building library matrix...\n');
polyOrder = 3;  % Up to cubic terms to capture x²y term
[Theta, library_names] = build_library_vdp(X_sparse, polyOrder);

fprintf('  Library order: %d\n', polyOrder);
fprintf('  Number of terms: %d\n', length(library_names));
fprintf('  Terms: %s\n', strjoin(library_names, ', '));

%% Step 3: Run SINDy (STLS)
fprintf('\nStep 3: Running SINDy (STLS)...\n');

% Set sparsity threshold (data-adaptive)
lambda_base = 0.01;
lambda = lambda_base * std(dXdt_fd(:));
max_iterations = 10;

fprintf('  Sparsity threshold: λ = %.6f (%.4f * std(dXdt))\n', lambda, lambda_base);
fprintf('  Max iterations: %d\n', max_iterations);

tic;
[Xi, sindy_stats] = run_sindy_stls(Theta, dXdt_fd, lambda, max_iterations);
sindy_time = toc;

fprintf('  SINDy completed in %.4f seconds\n', sindy_time);
fprintf('  Final sparsity: %.1f%% (%.0f/%d nonzero terms)\n', ...
    sindy_stats.final_sparsity * 100, sindy_stats.num_nonzero, sindy_stats.num_terms);
fprintf('  Prediction RMSE: %.4f\n', sindy_stats.rmse);

% Display discovered coefficients
fprintf('\n  Discovered coefficients:\n');
fprintf('    Term        | dx/dt      | dy/dt\n');
fprintf('    ------------|------------|------------\n');
for i = 1:length(library_names)
    fprintf('    %-10s | %10.4f | %10.4f\n', library_names{i}, Xi(i,1), Xi(i,2));
end

% True coefficients for Van der Pol (for reference)
% dx/dt = y → coefficient of 'y' = 1
% dy/dt = μ(1-x²)y - x = μy - μx²y - x
% → coefficient of 'x' = -1, 'y' = μ, 'x²y' = -μ
fprintf('\n  True coefficients (for reference):\n');
fprintf('    dx/dt = y\n');
fprintf('    dy/dt = %.2f(1-x²)y - x = %.2fy - %.2fx²y - x\n', ...
    params.mu, params.mu, params.mu);

%% Step 4: Forward Integration
fprintf('\nStep 4: Forward integrating discovered model...\n');

% Use same time span as ground truth
t_span = [min(t_true), max(t_true)];
x0 = params.x0;

fprintf('  Time span: [%.2f, %.2f]\n', t_span(1), t_span(2));
fprintf('  Initial conditions: x0 = [%.2f; %.2f]\n', x0(1), x0(2));

tic;
[t_sim, X_sim] = forward_integrate_sindy(Xi, library_names, t_span, x0);
integration_time = toc;

fprintf('  Integration completed in %.4f seconds\n', integration_time);
fprintf('  Simulated %d time points\n', length(t_sim));

% Interpolate to true time grid for comparison
X_sim_interp = interp1(t_sim, X_sim, t_true, 'pchip');

% Compute trajectory errors
traj_error = X_sim_interp - X_true;
rmse_traj = sqrt(mean(traj_error(:).^2));
mae_traj = mean(abs(traj_error(:)));

fprintf('  Trajectory error (vs ground truth): RMSE = %.4f, MAE = %.4f\n', ...
    rmse_traj, mae_traj);

%% Visualization
fprintf('\nCreating visualizations...\n');

figure('Color', 'w', 'Position', [100, 100, 1400, 1000]);

% Plot 1: Finite difference derivatives vs true
subplot(3, 3, 1);
plot(t_sparse, dXdt_true_at_sparse(:,1), 'o', 'Color', [0, 0.4470, 0.7410], ...
    'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'True dx/dt'); hold on;
plot(t_sparse, dXdt_fd(:,1), 's', 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'FD dx/dt');
xlabel('Time'); ylabel('dx/dt');
title('Finite Difference: dx/dt');
legend('Location', 'best'); grid on;

subplot(3, 3, 2);
plot(t_sparse, dXdt_true_at_sparse(:,2), 'o', 'Color', [0, 0.4470, 0.7410], ...
    'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'True dy/dt'); hold on;
plot(t_sparse, dXdt_fd(:,2), 's', 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'FD dy/dt');
xlabel('Time'); ylabel('dy/dt');
title('Finite Difference: dy/dt');
legend('Location', 'best'); grid on;

% Plot 2: SINDy coefficients
subplot(3, 3, 3);
bar(Xi);
set(gca, 'XTickLabel', library_names);
xtickangle(45);
ylabel('Coefficient Value');
title('SINDy Discovered Coefficients');
legend('dx/dt', 'dy/dt', 'Location', 'best');
grid on;

% Plot 3: State x(t) - SINDy vs True
subplot(3, 3, 4);
plot(t_true, X_true(:,1), '-', 'Color', [0, 0.4470, 0.7410], ...
    'LineWidth', 2, 'DisplayName', 'True x(t)'); hold on;
scatter(t_sparse, X_sparse(:,1), 100, [0.8500, 0.3250, 0.0980], ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Sparse data');
plot(t_sim, X_sim(:,1), '--', 'Color', [0.9290, 0.6940, 0.1250], ...
    'LineWidth', 2, 'DisplayName', 'SINDy model');
xlabel('Time'); ylabel('x(t)');
title('SINDy Model: x(t)');
legend('Location', 'best'); grid on;

% Plot 4: State y(t) - SINDy vs True
subplot(3, 3, 5);
plot(t_true, X_true(:,2), '-', 'Color', [0.8500, 0.3250, 0.0980], ...
    'LineWidth', 2, 'DisplayName', 'True y(t)'); hold on;
scatter(t_sparse, X_sparse(:,2), 100, [0, 0.4470, 0.7410], ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Sparse data');
plot(t_sim, X_sim(:,2), '--', 'Color', [0.9290, 0.6940, 0.1250], ...
    'LineWidth', 2, 'DisplayName', 'SINDy model');
xlabel('Time'); ylabel('y(t)');
title('SINDy Model: y(t)');
legend('Location', 'best'); grid on;

% Plot 5: Phase space
subplot(3, 3, 6);
plot(X_true(:,1), X_true(:,2), '-', 'Color', [0, 0.4470, 0.7410], ...
    'LineWidth', 2, 'DisplayName', 'True trajectory'); hold on;
scatter(X_sparse(:,1), X_sparse(:,2), 100, [0.8500, 0.3250, 0.0980], ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Sparse data');
plot(X_sim(:,1), X_sim(:,2), '--', 'Color', [0.9290, 0.6940, 0.1250], ...
    'LineWidth', 2, 'DisplayName', 'SINDy model');
xlabel('x'); ylabel('y');
title('Phase Space');
legend('Location', 'best'); grid on; axis equal;

% Plot 6: Prediction errors
subplot(3, 3, 7);
plot(t_true, traj_error(:,1), '-', 'Color', [0, 0.4470, 0.7410], ...
    'LineWidth', 1.5, 'DisplayName', 'x error'); hold on;
plot(t_true, traj_error(:,2), '-', 'Color', [0.8500, 0.3250, 0.0980], ...
    'LineWidth', 1.5, 'DisplayName', 'y error');
xlabel('Time'); ylabel('Prediction Error');
title('SINDy Trajectory Errors');
legend('Location', 'best'); grid on;
yline(0, 'k--', 'LineWidth', 1);

% Plot 7: Derivative prediction
subplot(3, 3, 8);
dXdt_sindy = Theta * Xi;
plot(t_sparse, dXdt_true_at_sparse(:,1), 'o', 'Color', [0, 0.4470, 0.7410], ...
    'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'True dx/dt'); hold on;
plot(t_sparse, dXdt_sindy(:,1), 's', 'Color', [0.9290, 0.6940, 0.1250], ...
    'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'SINDy dx/dt');
xlabel('Time'); ylabel('dx/dt');
title('SINDy Derivative: dx/dt');
legend('Location', 'best'); grid on;

subplot(3, 3, 9);
plot(t_sparse, dXdt_true_at_sparse(:,2), 'o', 'Color', [0.8500, 0.3250, 0.0980], ...
    'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'True dy/dt'); hold on;
plot(t_sparse, dXdt_sindy(:,2), 's', 'Color', [0.9290, 0.6940, 0.1250], ...
    'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'SINDy dy/dt');
xlabel('Time'); ylabel('dy/dt');
title('SINDy Derivative: dy/dt');
legend('Location', 'best'); grid on;

sgtitle(sprintf('Phase 3: SINDy-Only - %d sparse points, %.0f%% noise', ...
    n_sparse_points, noise_fraction*100), 'FontSize', 14, 'FontWeight', 'bold');

%% Save Phase 3 results
fprintf('\nSaving Phase 3 results...\n');
save('vdp_data.mat', 'dXdt_fd', 'Theta', 'library_names', 'Xi', 'sindy_stats', ...
    't_sim', 'X_sim', 'X_sim_interp', 'rmse_traj', 'mae_traj', 'rmse_fd', ...
    'fd_time', 'sindy_time', 'integration_time', '-append');

fprintf('Results saved to vdp_data.mat\n');
fprintf('\nPhase 3 complete! Ready for Phase 4 (GP+SINDy Implementation).\n');

