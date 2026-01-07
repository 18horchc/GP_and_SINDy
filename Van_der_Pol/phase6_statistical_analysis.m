%% Phase 6: Statistical Robustness Analysis
% This script runs multiple trials with different noise realizations
% to assess statistical robustness of the methods
%
% Requirements:
%   - All phase functions must be available
%   - This does NOT modify previous phase code

clear; close all; clc;

fprintf('=== Phase 6: Statistical Robustness Analysis ===\n');

%% Parameters
num_trials = 50;  % Number of noise realizations
n_sparse_points = 9;  % Fixed number of sparse points
noise_fraction = 0.10;  % Fixed noise level (10%)

% Fixed parameters
mu = 1.0;
x0 = [2; 0];
t_start = 0;
t_end = 10;
dt_fine = 0.01;
t_span = t_start:dt_fine:t_end;
sampling_method = 'uniform';
polyOrder = 3;

fprintf('\nConfiguration:\n');
fprintf('  Number of trials: %d\n', num_trials);
fprintf('  Sparse points: %d\n', n_sparse_points);
fprintf('  Noise level: %.0f%%\n', noise_fraction*100);

%% Initialize results storage
all_metrics = struct();
all_metrics.gp_only = struct();
all_metrics.sindy_only = struct();
all_metrics.gp_sindy = struct();

% Initialize arrays for each metric
metrics_list = {'rmse_traj', 'mae_traj', 'r2_traj', 'rmse_deriv', 'rmse_coeff', ...
    'mae_coeff', 'coeff_identification_rate', 'sparsity_pred'};

for method = {'gp_only', 'sindy_only', 'gp_sindy'}
    for metric = metrics_list
        all_metrics.(method{1}).(metric{1}) = zeros(num_trials, 1);
    end
end

%% Compute true coefficients for comparison
[Xi_true, library_names] = compute_true_vdp_coefficients(mu, polyOrder);

%% Run multiple trials
fprintf('\nRunning %d trials with different noise realizations...\n', num_trials);

successful_trials = 0;
failed_trials = 0;

for trial = 1:num_trials
    fprintf('Trial %d/%d...', trial, num_trials);
    
    try
        % Set random seed for reproducibility (but different for each trial)
        rng(trial);
        
        %% Phase 1: Generate data
        [t_true, X_true, dXdt_true, params] = generate_vdp_ground_truth(t_span, x0, mu);
        [t_sparse, X_sparse, ~, ~] = generate_sparse_noisy_data(...
            t_true, X_true, n_sparse_points, noise_fraction, sampling_method);
        
        %% Phase 2: GP Implementation
        t_dense = t_true;
        gp_options = struct();
        gp_options.kernel = 'squaredexponential';
        gp_options.standardize = false;
        gp_options.basis = 'constant';
        
        [gp_models, X_gp, dXdt_gp, gp_stats] = fit_gp_vdp(t_sparse, X_sparse, t_dense, gp_options);
        
        %% Phase 3: SINDy-Only
        dXdt_fd = compute_finite_differences(t_sparse, X_sparse, 'central');
        [Theta, ~] = build_library_vdp(X_sparse, polyOrder);
        
        lambda_base = 0.01;
        lambda = lambda_base * std(dXdt_fd(:));
        max_iterations = 10;
        
        [Xi, sindy_stats] = run_sindy_stls(Theta, dXdt_fd, lambda, max_iterations);
        
        t_span_sim = [min(t_true), max(t_true)];
        [t_sim, X_sim] = forward_integrate_sindy(Xi, library_names, t_span_sim, x0);
        X_sim_interp = interp1(t_sim, X_sim, t_true, 'pchip');
        
        %% Phase 4: GP+SINDy
        num_gp_samples = 1;
        augment_options = struct();
        augment_options.combine_method = 'append';
        
        [t_aug, X_aug, dXdt_aug, ~] = augment_data_with_gp_samples(...
            t_sparse, X_sparse, gp_models, t_dense, num_gp_samples, augment_options);
        
        [Theta_aug, ~] = build_library_vdp(X_aug, polyOrder);
        lambda_aug = lambda_base * std(dXdt_aug(:));
        
        [Xi_aug, sindy_stats_aug] = run_sindy_stls(Theta_aug, dXdt_aug, lambda_aug, max_iterations);
        
        [t_sim_aug, X_sim_aug] = forward_integrate_sindy(Xi_aug, library_names, t_span_sim, x0);
        X_sim_aug_interp = interp1(t_sim_aug, X_sim_aug, t_true, 'pchip');
        
        %% Validate results (check for numerical instability)
        max_reasonable_value = 1e6;
        if any(abs(X_sim(:)) > max_reasonable_value) || any(isnan(X_sim(:))) || any(isinf(X_sim(:))) || ...
           any(abs(X_sim_aug(:)) > max_reasonable_value) || any(isnan(X_sim_aug(:))) || any(isinf(X_sim_aug(:)))
            error('Integration unstable');
        end
        
        %% Compute metrics
        X_gp_interp = interp1(t_dense, X_gp, t_true, 'pchip');
        dXdt_gp_interp = interp1(t_dense, dXdt_gp, t_true, 'pchip');
        
        metrics_gp = compute_comprehensive_metrics(...
            X_true, dXdt_true, X_gp_interp, dXdt_gp_interp, ...
            [], [], library_names, 'GP-only');
        
        metrics_sindy = compute_comprehensive_metrics(...
            X_true, dXdt_true, X_sim_interp, [], ...
            Xi, Xi_true, library_names, 'SINDy-only');
        
        metrics_gp_sindy = compute_comprehensive_metrics(...
            X_true, dXdt_true, X_sim_aug_interp, [], ...
            Xi_aug, Xi_true, library_names, 'GP+SINDy');
        
        % Validate metrics before storing
        max_reasonable_rmse = 1e3;
        if metrics_sindy.rmse_traj > max_reasonable_rmse || metrics_gp_sindy.rmse_traj > max_reasonable_rmse || ...
           isnan(metrics_sindy.rmse_traj) || isnan(metrics_gp_sindy.rmse_traj)
            error('Metrics indicate unstable solution');
        end
        
        % Store metrics
        for metric = metrics_list
            metric_name = metric{1};
            if isfield(metrics_gp, metric_name)
                all_metrics.gp_only.(metric_name)(trial) = metrics_gp.(metric_name);
            end
            if isfield(metrics_sindy, metric_name)
                all_metrics.sindy_only.(metric_name)(trial) = metrics_sindy.(metric_name);
            end
            if isfield(metrics_gp_sindy, metric_name)
                all_metrics.gp_sindy.(metric_name)(trial) = metrics_gp_sindy.(metric_name);
            end
        end
        
        successful_trials = successful_trials + 1;
        fprintf(' ✓\n');
        
    catch ME
        failed_trials = failed_trials + 1;
        fprintf(' ✗ Error: %s\n', ME.message);
        % Store NaN for failed trial
        for method = {'gp_only', 'sindy_only', 'gp_sindy'}
            for metric = metrics_list
                all_metrics.(method{1}).(metric{1})(trial) = NaN;
            end
        end
    end
end

fprintf('\nCompleted: %d successful, %d failed\n', successful_trials, failed_trials);

%% Compute statistics
fprintf('\n=== Statistical Summary ===\n');

% Remove NaN values for statistics
for method = {'gp_only', 'sindy_only', 'gp_sindy'}
    for metric = metrics_list
        values = all_metrics.(method{1}).(metric{1});
        valid_idx = ~isnan(values);
        if sum(valid_idx) > 0
            all_metrics.(method{1}).(metric{1}) = values(valid_idx);
        end
    end
end

% Display summary statistics
fprintf('\nTrajectory RMSE (mean ± std):\n');
fprintf('%-15s | %10s | %10s | %10s\n', 'Method', 'Mean', 'Std', 'N');
fprintf('%s\n', repmat('-', 1, 50));

for method = {'gp_only', 'sindy_only', 'gp_sindy'}
    method_name = strrep(method{1}, '_', '-');
    values = all_metrics.(method{1}).rmse_traj;
    valid_values = values(~isnan(values));
    if ~isempty(valid_values)
        mean_val = round(mean(valid_values), 4);
        std_val = round(std(valid_values), 4);
        fprintf('%-15s | %10.4f | %10.4f | %10d\n', ...
            method_name, mean_val, std_val, length(valid_values));
    end
end

fprintf('\nTrajectory R² (mean ± std):\n');
fprintf('%-15s | %10s | %10s | %10s\n', 'Method', 'Mean', 'Std', 'N');
fprintf('%s\n', repmat('-', 1, 50));

for method = {'gp_only', 'sindy_only', 'gp_sindy'}
    method_name = strrep(method{1}, '_', '-');
    values = all_metrics.(method{1}).r2_traj;
    valid_values = values(~isnan(values));
    if ~isempty(valid_values)
        mean_val = round(mean(valid_values), 4);
        std_val = round(std(valid_values), 4);
        fprintf('%-15s | %10.4f | %10.4f | %10d\n', ...
            method_name, mean_val, std_val, length(valid_values));
    end
end

fprintf('\nCoefficient RMSE (mean ± std):\n');
fprintf('%-15s | %10s | %10s | %10s\n', 'Method', 'Mean', 'Std', 'N');
fprintf('%s\n', repmat('-', 1, 50));

for method = {'sindy_only', 'gp_sindy'}
    method_name = strrep(method{1}, '_', '-');
    values = all_metrics.(method{1}).rmse_coeff;
    valid_values = values(~isnan(values));
    if ~isempty(valid_values)
        mean_val = round(mean(valid_values), 4);
        std_val = round(std(valid_values), 4);
        fprintf('%-15s | %10.4f | %10.4f | %10d\n', ...
            method_name, mean_val, std_val, length(valid_values));
    end
end

%% Improvement analysis
fprintf('\n=== Improvement Analysis ===\n');
fprintf('GP+SINDy vs SINDy-only:\n');

sindy_rmse = all_metrics.sindy_only.rmse_traj;
gp_sindy_rmse = all_metrics.gp_sindy.rmse_traj;

% Compute improvement for each trial
improvements = (1 - gp_sindy_rmse ./ sindy_rmse) * 100;
improvements = improvements(~isnan(improvements) & ~isinf(improvements));

fprintf('  RMSE improvement: %.2f%% ± %.2f%% (mean ± std)\n', ...
    round(mean(improvements), 2), round(std(improvements), 2));
fprintf('  Success rate: %.1f%% (%d/%d trials where GP+SINDy < SINDy-only)\n', ...
    round(sum(improvements > 0) / length(improvements) * 100, 1), ...
    sum(improvements > 0), length(improvements));

% Statistical test (paired t-test)
if length(improvements) > 1
    [h, p] = ttest(sindy_rmse, gp_sindy_rmse);
    fprintf('  Paired t-test: p = %.4f (h = %d, significant if p < 0.05)\n', round(p, 4), h);
end

%% Visualization
fprintf('\nCreating visualizations...\n');

figure('Color', 'w', 'Position', [100, 100, 1600, 1000]);

% Plot 1: RMSE distribution
subplot(2, 3, 1);
hold on;
histogram(all_metrics.gp_only.rmse_traj, 'FaceAlpha', 0.5, 'DisplayName', 'GP-only');
histogram(all_metrics.sindy_only.rmse_traj, 'FaceAlpha', 0.5, 'DisplayName', 'SINDy-only');
histogram(all_metrics.gp_sindy.rmse_traj, 'FaceAlpha', 0.5, 'DisplayName', 'GP+SINDy');
xlabel('Trajectory RMSE');
ylabel('Frequency');
title('RMSE Distribution');
legend('Location', 'best');
grid on;

% Plot 2: R² distribution
subplot(2, 3, 2);
hold on;
histogram(all_metrics.gp_only.r2_traj, 'FaceAlpha', 0.5, 'DisplayName', 'GP-only');
histogram(all_metrics.sindy_only.r2_traj, 'FaceAlpha', 0.5, 'DisplayName', 'SINDy-only');
histogram(all_metrics.gp_sindy.r2_traj, 'FaceAlpha', 0.5, 'DisplayName', 'GP+SINDy');
xlabel('Trajectory R²');
ylabel('Frequency');
title('R² Distribution');
legend('Location', 'best');
grid on;

% Plot 3: Box plot comparison
subplot(2, 3, 3);
rmse_data = [all_metrics.gp_only.rmse_traj; ...
             all_metrics.sindy_only.rmse_traj; ...
             all_metrics.gp_sindy.rmse_traj];
group = [ones(length(all_metrics.gp_only.rmse_traj), 1); ...
         2*ones(length(all_metrics.sindy_only.rmse_traj), 1); ...
         3*ones(length(all_metrics.gp_sindy.rmse_traj), 1)];
boxplot(rmse_data, group, 'Labels', {'GP-only', 'SINDy-only', 'GP+SINDy'});
ylabel('Trajectory RMSE');
title('RMSE Box Plot Comparison');
grid on;

% Plot 4: Improvement distribution
subplot(2, 3, 4);
histogram(improvements, 20, 'FaceColor', [0.9290, 0.6940, 0.1250]);
xlabel('Improvement (%)');
ylabel('Frequency');
title('Improvement Distribution (GP+SINDy vs SINDy-only)');
xline(0, 'r--', 'LineWidth', 2, 'DisplayName', 'No improvement');
grid on;
legend;

% Plot 5: Scatter plot: SINDy-only vs GP+SINDy RMSE
subplot(2, 3, 5);
scatter(sindy_rmse, gp_sindy_rmse, 50, 'filled', 'MarkerEdgeColor', 'k');
hold on;
plot([0, max([sindy_rmse; gp_sindy_rmse])], [0, max([sindy_rmse; gp_sindy_rmse])], ...
    'r--', 'LineWidth', 2, 'DisplayName', 'y=x (no improvement)');
xlabel('SINDy-only RMSE');
ylabel('GP+SINDy RMSE');
title('Method Comparison: RMSE');
legend('Location', 'best');
grid on;
axis equal;

% Plot 6: Coefficient RMSE comparison
subplot(2, 3, 6);
coeff_rmse_data = [all_metrics.sindy_only.rmse_coeff; ...
                   all_metrics.gp_sindy.rmse_coeff];
coeff_group = [ones(length(all_metrics.sindy_only.rmse_coeff), 1); ...
               2*ones(length(all_metrics.gp_sindy.rmse_coeff), 1)];
boxplot(coeff_rmse_data, coeff_group, 'Labels', {'SINDy-only', 'GP+SINDy'});
ylabel('Coefficient RMSE');
title('Coefficient RMSE Comparison');
grid on;

sgtitle(sprintf('Phase 6: Statistical Analysis (%d trials, %d points, %.0f%% noise)', ...
    num_trials, n_sparse_points, noise_fraction*100), 'FontSize', 14, 'FontWeight', 'bold');

%% Save results
fprintf('\nSaving results...\n');
save('phase6_statistical_analysis_results.mat', 'all_metrics', 'num_trials', ...
    'n_sparse_points', 'noise_fraction', 'improvements');
fprintf('Results saved to phase6_statistical_analysis_results.mat\n');

fprintf('\nPhase 6 Statistical Analysis complete!\n');

