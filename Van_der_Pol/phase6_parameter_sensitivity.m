%% Phase 6: Parameter Sensitivity Analysis
% This script tests how results vary with different parameters:
% - Number of sparse data points
% - Noise level
%
% Requirements:
%   - All phase functions must be available
%   - This does NOT modify previous phase code

clear; close all; clc;

fprintf('=== Phase 6: Parameter Sensitivity Analysis ===\n');

%% Parameter ranges to test
sparse_point_counts = [5, 7, 9, 11, 15];  % Number of sparse data points
noise_levels = [0.05, 0.10, 0.15, 0.20];  % Noise levels (5%, 10%, 15%, 20%)

% Fixed parameters
mu = 1.0;  % Van der Pol parameter
x0 = [2; 0];  % Initial conditions
t_start = 0;
t_end = 10;
dt_fine = 0.01;  % Fine resolution for ground truth
t_span = t_start:dt_fine:t_end;
sampling_method = 'uniform';

fprintf('\nParameter ranges:\n');
fprintf('  Sparse points: %s\n', mat2str(sparse_point_counts));
fprintf('  Noise levels: %s\n', mat2str(noise_levels * 100));
fprintf('  Total combinations: %d\n', length(sparse_point_counts) * length(noise_levels));

%% Initialize results storage
results = struct();
results.sparse_points = sparse_point_counts;
results.noise_levels = noise_levels;
results.metrics = cell(length(sparse_point_counts), length(noise_levels));

%% Run parameter sweep
fprintf('\nRunning parameter sweep...\n');
total_combinations = length(sparse_point_counts) * length(noise_levels);
current_combination = 0;

for i = 1:length(sparse_point_counts)
    n_sparse = sparse_point_counts(i);
    
    for j = 1:length(noise_levels)
        noise_frac = noise_levels(j);
        current_combination = current_combination + 1;
        
        fprintf('\n[%d/%d] Testing: %d points, %.0f%% noise\n', ...
            current_combination, total_combinations, n_sparse, noise_frac*100);
        
        try
            %% Phase 1: Generate data
            [t_true, X_true, dXdt_true, params] = generate_vdp_ground_truth(t_span, x0, mu);
            [t_sparse, X_sparse, X_true_at_sparse, noise_level] = ...
                generate_sparse_noisy_data(t_true, X_true, n_sparse, noise_frac, sampling_method);
            
            %% Phase 2: GP Implementation
            t_dense = t_true;
            gp_options = struct();
            gp_options.kernel = 'squaredexponential';
            gp_options.standardize = false;
            gp_options.basis = 'constant';
            
            [gp_models, X_gp, dXdt_gp, gp_stats] = fit_gp_vdp(t_sparse, X_sparse, t_dense, gp_options);
            
            %% Phase 3: SINDy-Only
            dXdt_fd = compute_finite_differences(t_sparse, X_sparse, 'central');
            polyOrder = 3;
            [Theta, library_names] = build_library_vdp(X_sparse, polyOrder);
            
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
            
            [t_aug, X_aug, dXdt_aug, aug_stats] = augment_data_with_gp_samples(...
                t_sparse, X_sparse, gp_models, t_dense, num_gp_samples, augment_options);
            
            [Theta_aug, ~] = build_library_vdp(X_aug, polyOrder);
            lambda_aug = lambda_base * std(dXdt_aug(:));
            
            [Xi_aug, sindy_stats_aug] = run_sindy_stls(Theta_aug, dXdt_aug, lambda_aug, max_iterations);
            
            [t_sim_aug, X_sim_aug] = forward_integrate_sindy(Xi_aug, library_names, t_span_sim, x0);
            X_sim_aug_interp = interp1(t_sim_aug, X_sim_aug, t_true, 'pchip');
            
            %% Validate results (check for numerical instability)
            % Check if integration produced reasonable values
            max_reasonable_value = 1e6;  % Maximum reasonable state value
            if any(abs(X_sim(:)) > max_reasonable_value) || any(isnan(X_sim(:))) || any(isinf(X_sim(:)))
                error('SINDy-only integration unstable');
            end
            if any(abs(X_sim_aug(:)) > max_reasonable_value) || any(isnan(X_sim_aug(:))) || any(isinf(X_sim_aug(:)))
                error('GP+SINDy integration unstable');
            end
            
            %% Compute metrics
            [Xi_true, ~] = compute_true_vdp_coefficients(mu, polyOrder);
            
            % Interpolate to common time grid
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
            
            % Validate metrics are reasonable
            max_reasonable_rmse = 1e3;  % Maximum reasonable RMSE
            if metrics_sindy.rmse_traj > max_reasonable_rmse || ...
               metrics_gp_sindy.rmse_traj > max_reasonable_rmse || ...
               isnan(metrics_sindy.rmse_traj) || isnan(metrics_gp_sindy.rmse_traj)
                error('Metrics indicate unstable solution');
            end
            
            % Store results
            results.metrics{i, j} = struct();
            results.metrics{i, j}.gp = metrics_gp;
            results.metrics{i, j}.sindy = metrics_sindy;
            results.metrics{i, j}.gp_sindy = metrics_gp_sindy;
            
            fprintf('  ✓ Completed successfully\n');
            
        catch ME
            fprintf('  ✗ Error: %s\n', ME.message);
            results.metrics{i, j} = struct();  % Empty struct for failed run
        end
    end
end

%% Display results summary
fprintf('\n=== Parameter Sensitivity Results Summary ===\n');

% Create summary tables
fprintf('\nTrajectory RMSE by Parameter Combination:\n');
fprintf('\n%-10s |', 'Noise %');
for j = 1:length(noise_levels)
    fprintf(' %6.0f%% |', noise_levels(j)*100);
end
fprintf('\n');
fprintf('%s\n', repmat('-', 1, 10 + 9*length(noise_levels)));

for i = 1:length(sparse_point_counts)
    fprintf('%-10d |', sparse_point_counts(i));
    for j = 1:length(noise_levels)
        if ~isempty(results.metrics{i, j}) && isfield(results.metrics{i, j}, 'gp_sindy')
            rmse = results.metrics{i, j}.gp_sindy.rmse_traj;
            if ~isnan(rmse) && ~isinf(rmse) && rmse < 1e3
                rmse = round(rmse, 4);
                fprintf(' %7.4f |', rmse);
            else
                fprintf(' %7s |', 'Unstable');
            end
        else
            fprintf(' %7s |', 'N/A');
        end
    end
    fprintf('\n');
end

fprintf('\nImprovement (GP+SINDy vs SINDy-only) by Parameter Combination:\n');
fprintf('\n%-10s |', 'Noise %');
for j = 1:length(noise_levels)
    fprintf(' %6.0f%% |', noise_levels(j)*100);
end
fprintf('\n');
fprintf('%s\n', repmat('-', 1, 10 + 9*length(noise_levels)));

for i = 1:length(sparse_point_counts)
    fprintf('%-10d |', sparse_point_counts(i));
    for j = 1:length(noise_levels)
        if ~isempty(results.metrics{i, j}) && isfield(results.metrics{i, j}, 'gp_sindy') && ...
           isfield(results.metrics{i, j}, 'sindy')
            sindy_rmse = results.metrics{i, j}.sindy.rmse_traj;
            gp_sindy_rmse = results.metrics{i, j}.gp_sindy.rmse_traj;
            if sindy_rmse > 0 && ~isnan(sindy_rmse) && ~isnan(gp_sindy_rmse) && ...
               ~isinf(sindy_rmse) && ~isinf(gp_sindy_rmse) && ...
               sindy_rmse < 1e3 && gp_sindy_rmse < 1e3
                improvement = round((1 - gp_sindy_rmse / sindy_rmse) * 100, 1);
                % Cap improvement display at reasonable range
                if abs(improvement) > 1000
                    fprintf(' %7s |', 'Unstable');
                else
                    fprintf(' %6.1f%% |', improvement);
                end
            else
                fprintf(' %7s |', 'N/A');
            end
        else
            fprintf(' %7s |', 'N/A');
        end
    end
    fprintf('\n');
end

%% Visualization
fprintf('\nCreating visualizations...\n');

figure('Color', 'w', 'Position', [100, 100, 1600, 1000]);

% Extract data for plotting
rmse_gp = zeros(length(sparse_point_counts), length(noise_levels));
rmse_sindy = zeros(length(sparse_point_counts), length(noise_levels));
rmse_gp_sindy = zeros(length(sparse_point_counts), length(noise_levels));
improvement = zeros(length(sparse_point_counts), length(noise_levels));

for i = 1:length(sparse_point_counts)
    for j = 1:length(noise_levels)
        if ~isempty(results.metrics{i, j}) && isfield(results.metrics{i, j}, 'gp_sindy')
            rmse_gp(i, j) = results.metrics{i, j}.gp.rmse_traj;
            rmse_sindy(i, j) = results.metrics{i, j}.sindy.rmse_traj;
            rmse_gp_sindy(i, j) = results.metrics{i, j}.gp_sindy.rmse_traj;
            if rmse_sindy(i, j) > 0
                improvement(i, j) = (1 - rmse_gp_sindy(i, j) / rmse_sindy(i, j)) * 100;
            end
        else
            rmse_gp(i, j) = NaN;
            rmse_sindy(i, j) = NaN;
            rmse_gp_sindy(i, j) = NaN;
            improvement(i, j) = NaN;
        end
    end
end

% Plot 1: RMSE heatmap for GP+SINDy
subplot(2, 3, 1);
imagesc(noise_levels*100, sparse_point_counts, rmse_gp_sindy);
colorbar;
xlabel('Noise Level (%)');
ylabel('Number of Sparse Points');
title('GP+SINDy: Trajectory RMSE');
set(gca, 'YDir', 'normal');
colormap(gca, 'hot');

% Plot 2: RMSE heatmap for SINDy-only
subplot(2, 3, 2);
imagesc(noise_levels*100, sparse_point_counts, rmse_sindy);
colorbar;
xlabel('Noise Level (%)');
ylabel('Number of Sparse Points');
title('SINDy-only: Trajectory RMSE');
set(gca, 'YDir', 'normal');
colormap(gca, 'hot');

% Plot 3: Improvement heatmap
subplot(2, 3, 3);
imagesc(noise_levels*100, sparse_point_counts, improvement);
colorbar;
xlabel('Noise Level (%)');
ylabel('Number of Sparse Points');
title('Improvement: GP+SINDy vs SINDy-only (%)');
set(gca, 'YDir', 'normal');
colormap(gca, 'cool');

% Plot 4: RMSE vs sparse points (different noise levels)
subplot(2, 3, 4);
hold on;
for j = 1:length(noise_levels)
    plot(sparse_point_counts, rmse_gp_sindy(:, j), '-o', 'LineWidth', 2, ...
        'DisplayName', sprintf('GP+SINDy, %.0f%% noise', noise_levels(j)*100));
    plot(sparse_point_counts, rmse_sindy(:, j), '--s', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('SINDy-only, %.0f%% noise', noise_levels(j)*100));
end
xlabel('Number of Sparse Points');
ylabel('Trajectory RMSE');
title('RMSE vs Sparse Points');
legend('Location', 'best', 'FontSize', 8);
grid on;

% Plot 5: RMSE vs noise level (different sparse point counts)
subplot(2, 3, 5);
hold on;
for i = 1:length(sparse_point_counts)
    plot(noise_levels*100, rmse_gp_sindy(i, :), '-o', 'LineWidth', 2, ...
        'DisplayName', sprintf('GP+SINDy, %d points', sparse_point_counts(i)));
    plot(noise_levels*100, rmse_sindy(i, :), '--s', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('SINDy-only, %d points', sparse_point_counts(i)));
end
xlabel('Noise Level (%)');
ylabel('Trajectory RMSE');
title('RMSE vs Noise Level');
legend('Location', 'best', 'FontSize', 8);
grid on;

% Plot 6: Improvement vs parameters
subplot(2, 3, 6);
imagesc(noise_levels*100, sparse_point_counts, improvement);
colorbar;
xlabel('Noise Level (%)');
ylabel('Number of Sparse Points');
title('Improvement (%)');
set(gca, 'YDir', 'normal');
colormap(gca, 'cool');
% Add text annotations
for i = 1:length(sparse_point_counts)
    for j = 1:length(noise_levels)
        if ~isnan(improvement(i, j))
            text(noise_levels(j)*100, sparse_point_counts(i), ...
                sprintf('%.1f%%', improvement(i, j)), ...
                'HorizontalAlignment', 'center', 'Color', 'white', 'FontWeight', 'bold');
        end
    end
end

sgtitle('Phase 6: Parameter Sensitivity Analysis', 'FontSize', 14, 'FontWeight', 'bold');

%% Save results
fprintf('\nSaving results...\n');
save('phase6_parameter_sensitivity_results.mat', 'results', 'sparse_point_counts', 'noise_levels');
fprintf('Results saved to phase6_parameter_sensitivity_results.mat\n');

fprintf('\nPhase 6 Parameter Sensitivity Analysis complete!\n');

