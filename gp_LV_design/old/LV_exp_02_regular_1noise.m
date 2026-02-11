%% GP Lotka-Volterra Experiment 02: Regular sampling, 1% noise, varying N and kernel
%
% Same structure as LV_exp_01 but with 1% Gaussian noise on sampled prey and predator.
% Uses same parameters as lotka_volterra_model.m / all_toy_problems.m (defaults: alpha=1, beta=0.2, delta=0.5, gamma=0.2, x0=1, y0=2, tspan=[0 50]).
%
% Scope:
%   - Ground truth: 500 pts on [0, 50] from LV (two states: prey, predator)
%   - Sampling: REGULAR only, 1% Gaussian noise (sigma = 0.01 * std per state)
%   - Sparsity N: 5, 10, 25, 50
%   - Kernels: 5 built-in + custom periodic (6 total)
%   - One GP per state (prey and predator); same metrics as logistic design
%
% Output: Three metric tables (Kernel, N, State, ...); RMSE vs N figure; one figure per kernel
%         with 2 rows (prey, predator) x 4 subplots (N=5,10,25,50): points, ground truth, GP mean + 95%% band.

clear; clc; close all;

% This folder (lotka_volterra_model, ground_truth_LV) and parent (metric_helpers)
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..'));

%% 1. Ground truth (500 points, [0, 50]) — same setup as lotka_volterra_model defaults
[t_gt, prey_gt, pred_gt] = ground_truth_LV(500);
fprintf('Ground truth LV: %d points on [0, 50]. Prey [%.2f, %.2f], Predator [%.2f, %.2f]\n', ...
    length(t_gt), min(prey_gt), max(prey_gt), min(pred_gt), max(pred_gt));

%% 2. Design (5 built-in + custom periodic kernel)
N_list = [5, 10, 25, 50];
kernel_list = {'squaredexponential', 'exponential', 'matern32', 'matern52', 'rationalquadratic', @periodicKernel};
kernel_labels = {'SqExp', 'Matern 1/2', 'Matern 3/2', 'Matern 5/2', 'RatQuad', 'Periodic'};
n_anticipated_cycles = 3;  % for periodic kernel initial period

% 1% noise: sigma = 0.01 * std(ground truth) per state, per design plan
noise_pct = 0.01;
sigma_noise_prey = noise_pct * std(prey_gt);
sigma_noise_pred = noise_pct * std(pred_gt);

rng(42);  % seed for reproducibility

%% 3. Loop: each N → sample regularly, 1% noise; each kernel → fit GP for prey and predator, compute metrics
results_point = [];
results_prob = [];
results_calib = [];
pred_ymu_prey  = cell(length(kernel_list), length(N_list));
pred_ystd_prey = cell(length(kernel_list), length(N_list));
pred_ymu_pred  = cell(length(kernel_list), length(N_list));
pred_ystd_pred = cell(length(kernel_list), length(N_list));
t_obs_by_N   = cell(1, length(N_list));
prey_obs_by_N = cell(1, length(N_list));
pred_obs_by_N = cell(1, length(N_list));

for in = 1:length(N_list)
    N = N_list(in);
    t_obs = linspace(0, 50, N)';
    prey_obs = interp1(t_gt, prey_gt, t_obs, 'linear', 'extrap');
    pred_obs = interp1(t_gt, pred_gt, t_obs, 'linear', 'extrap');
    % Add 1% Gaussian noise
    prey_obs = prey_obs + sigma_noise_prey * randn(size(prey_obs));
    pred_obs = pred_obs + sigma_noise_pred * randn(size(pred_obs));
    t_obs_by_N{in} = t_obs;
    prey_obs_by_N{in} = prey_obs;
    pred_obs_by_N{in} = pred_obs;

    % Initial parameters for periodic kernel (same period/lengthScale for both states)
    t_range = max(t_obs) - min(t_obs);
    period0 = t_range / n_anticipated_cycles;
    lengthScale0 = period0 / 4;

    for ik = 1:length(kernel_list)
        kern = kernel_list{ik};
        is_periodic = isa(kern, 'function_handle');
        try
            if is_periodic
                theta0_prey = [std(prey_obs), period0, lengthScale0];
                theta0_pred = [std(pred_obs), period0, lengthScale0];
            end
            % GP for prey
            if is_periodic
                gpr_prey = fitrgp(t_obs, prey_obs, ...
                    'KernelFunction', @periodicKernel, 'KernelParameters', theta0_prey, ...
                    'BasisFunction', 'constant', 'FitMethod', 'exact', 'PredictMethod', 'exact', 'Standardize', true);
            else
                gpr_prey = fitrgp(t_obs, prey_obs, ...
                    'KernelFunction', kern, 'BasisFunction', 'constant', ...
                    'FitMethod', 'exact', 'PredictMethod', 'exact', 'Standardize', true);
            end
            [ymu_prey, ystd_prey, ~] = predict(gpr_prey, t_gt);
            pred_ymu_prey{ik, in} = ymu_prey;
            pred_ystd_prey{ik, in} = ystd_prey;

            [rmse_p, mae_p, r2_p] = metric_helpers('point_estimate', prey_gt, ymu_prey);
            [nlpd_p, msll_p, crps_p] = metric_helpers('probabilistic', prey_gt, ymu_prey, ystd_prey, prey_obs);
            [smse_p, cov_p, nlml_p] = metric_helpers('calibration', prey_gt, ymu_prey, ystd_prey, prey_obs, gpr_prey);

            results_point = [results_point; {kernel_labels{ik}, N, 'Prey', rmse_p, mae_p, r2_p}];
            results_prob  = [results_prob;  {kernel_labels{ik}, N, 'Prey', nlpd_p, msll_p, crps_p}];
            results_calib = [results_calib; {kernel_labels{ik}, N, 'Prey', smse_p, cov_p, nlml_p}];

            % GP for predator
            if is_periodic
                gpr_pred = fitrgp(t_obs, pred_obs, ...
                    'KernelFunction', @periodicKernel, 'KernelParameters', theta0_pred, ...
                    'BasisFunction', 'constant', 'FitMethod', 'exact', 'PredictMethod', 'exact', 'Standardize', true);
            else
                gpr_pred = fitrgp(t_obs, pred_obs, ...
                    'KernelFunction', kern, 'BasisFunction', 'constant', ...
                    'FitMethod', 'exact', 'PredictMethod', 'exact', 'Standardize', true);
            end
            [ymu_pred, ystd_pred, ~] = predict(gpr_pred, t_gt);
            pred_ymu_pred{ik, in} = ymu_pred;
            pred_ystd_pred{ik, in} = ystd_pred;

            [rmse_y, mae_y, r2_y] = metric_helpers('point_estimate', pred_gt, ymu_pred);
            [nlpd_y, msll_y, crps_y] = metric_helpers('probabilistic', pred_gt, ymu_pred, ystd_pred, pred_obs);
            [smse_y, cov_y, nlml_y] = metric_helpers('calibration', pred_gt, ymu_pred, ystd_pred, pred_obs, gpr_pred);

            results_point = [results_point; {kernel_labels{ik}, N, 'Predator', rmse_y, mae_y, r2_y}];
            results_prob  = [results_prob;  {kernel_labels{ik}, N, 'Predator', nlpd_y, msll_y, crps_y}];
            results_calib = [results_calib; {kernel_labels{ik}, N, 'Predator', smse_y, cov_y, nlml_y}];
        catch me
            warning('Failed: N=%d, kernel=%s. %s', N, kernel_labels{ik}, me.message);
            results_point = [results_point; {kernel_labels{ik}, N, 'Prey', NaN, NaN, NaN}; {kernel_labels{ik}, N, 'Predator', NaN, NaN, NaN}];
            results_prob  = [results_prob;  {kernel_labels{ik}, N, 'Prey', NaN, NaN, NaN}; {kernel_labels{ik}, N, 'Predator', NaN, NaN, NaN}];
            results_calib = [results_calib; {kernel_labels{ik}, N, 'Prey', NaN, NaN, NaN}; {kernel_labels{ik}, N, 'Predator', NaN, NaN, NaN}];
            pred_ymu_prey{ik, in} = nan(size(t_gt)); pred_ystd_prey{ik, in} = nan(size(t_gt));
            pred_ymu_pred{ik, in} = nan(size(t_gt)); pred_ystd_pred{ik, in} = nan(size(t_gt));
        end
    end
end

%% 4. Tables (include State: Prey / Predator)
T_point = cell2table(results_point, 'VariableNames', {'Kernel', 'N', 'State', 'RMSE', 'MAE', 'R2'});
T_prob  = cell2table(results_prob,  'VariableNames', {'Kernel', 'N', 'State', 'NLPD', 'MSLL', 'CRPS'});
T_calib = cell2table(results_calib, 'VariableNames', {'Kernel', 'N', 'State', 'sMSE', 'Coverage', 'NLML'});

fprintf('--- Point-Estimate Metrics (Accuracy) ---\n');
disp(T_point);
fprintf('--- Probabilistic Metrics (Uncertainty Quality) ---\n');
disp(T_prob);
fprintf('--- Calibration & Robustness Metrics ---\n');
disp(T_calib);

%% 5. RMSE vs N: one line per kernel, separate by state (or average); show Prey and Predator in one figure with two subplots
figure;
subplot(1,2,1);
hold on;
for ik = 1:length(kernel_labels)
    idx = strcmp(T_point.Kernel, kernel_labels{ik}) & strcmp(T_point.State, 'Prey');
    plot(T_point.N(idx), T_point.RMSE(idx), '-o', 'LineWidth', 1.5, 'DisplayName', kernel_labels{ik});
end
xlabel('N'); ylabel('RMSE'); title('LV Exp 02: Prey — Regular, 1% noise');
legend('Location', 'best'); grid on; set(gca, 'XTick', N_list);
subplot(1,2,2);
hold on;
for ik = 1:length(kernel_labels)
    idx = strcmp(T_point.Kernel, kernel_labels{ik}) & strcmp(T_point.State, 'Predator');
    plot(T_point.N(idx), T_point.RMSE(idx), '-o', 'LineWidth', 1.5, 'DisplayName', kernel_labels{ik});
end
xlabel('N'); ylabel('RMSE'); title('LV Exp 02: Predator — Regular, 1% noise');
legend('Location', 'best'); grid on; set(gca, 'XTick', N_list);
sgtitle('GP LV Exp 02: Regular sampling, 1% noise');

%% 6. One figure per kernel: 2 rows (prey, predator) x 4 cols (N = 5, 10, 25, 50)
for ik = 1:length(kernel_labels)
    figure('Name', sprintf('LV Exp 02: %s', kernel_labels{ik}));
    for in = 1:length(N_list)
        % Row 1: Prey
        subplot(2, 4, in);
        t_obs = t_obs_by_N{in};
        prey_obs = prey_obs_by_N{in};
        ymu = pred_ymu_prey{ik, in};
        ystd = pred_ystd_prey{ik, in};
        plot(t_gt, prey_gt, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Ground truth');
        hold on;
        plot(t_obs, prey_obs, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
        plot(t_gt, ymu, 'b-', 'LineWidth', 1.2, 'DisplayName', 'GP mean');
        if ~all(isnan(ystd))
            ylo = ymu - 1.96*ystd; yhi = ymu + 1.96*ystd;
            fill([t_gt; flip(t_gt)], [ylo; flip(yhi)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        hold off;
        xlabel('Time'); ylabel('Prey');
        title(sprintf('Prey, N = %d', N_list(in)));
        legend('Location', 'best', 'FontSize', 8);
        grid on; xlim([0 50]);

        % Row 2: Predator
        subplot(2, 4, 4 + in);
        pred_obs = pred_obs_by_N{in};
        ymu = pred_ymu_pred{ik, in};
        ystd = pred_ystd_pred{ik, in};
        plot(t_gt, pred_gt, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Ground truth');
        hold on;
        plot(t_obs, pred_obs, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
        plot(t_gt, ymu, 'r-', 'LineWidth', 1.2, 'DisplayName', 'GP mean');
        if ~all(isnan(ystd))
            ylo = ymu - 1.96*ystd; yhi = ymu + 1.96*ystd;
            fill([t_gt; flip(t_gt)], [ylo; flip(yhi)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        hold off;
        xlabel('Time'); ylabel('Predator');
        title(sprintf('Predator, N = %d', N_list(in)));
        legend('Location', 'best', 'FontSize', 8);
        grid on; xlim([0 50]);
    end
    sgtitle(sprintf('GP LV Exp 02: %s — Regular, 1%% noise', kernel_labels{ik}));
end

%% 7. Save (optional)
save(fullfile(script_dir, 'results_exp_02_regular_1noise.mat'), 'T_point', 'T_prob', 'T_calib', 't_gt', 'prey_gt', 'pred_gt');
writetable(T_point, fullfile(script_dir, 'results_exp_02_point_metrics.csv'));
writetable(T_prob, fullfile(script_dir, 'results_exp_02_prob_metrics.csv'));
writetable(T_calib, fullfile(script_dir, 'results_exp_02_calib_metrics.csv'));

fprintf('Done. Total rows: %d (Prey + Predator per Kernel,N)\n', height(T_point));
