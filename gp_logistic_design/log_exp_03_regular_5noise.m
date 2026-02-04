%% GP Logistic Experiment 03: Regular sampling, 5% noise, varying N and kernel
%
% Part of the GP-on-logistic-growth experimental design (see GP_Logistic_Experimental_Design_Plan.md).
%
% Scope:
%   - Ground truth: 500 pts on [0, 50] from logistic growth
%   - Sampling: REGULAR only
%   - Noise: 5% Gaussian (sigma = 0.05 * std(y_gt))
%   - Sparsity N: 5, 10, 25, 50
%   - Kernels: all 5 (Squared Exp, Matern 1/2, 3/2, 5/2, Rational Quadratic)
%   - Optimization: default (no OptimizeHyperparameters)
%
% Naming: log_exp_NN_<sampling>_<noise>.m — NN = experiment number, build on as we add steps.
%
% Output: Three metric tables (Point-Estimate, Probabilistic, Calibration);
%         figure RMSE vs N by kernel;
%         one figure per kernel with 4 subplots (N=5,10,25,50): points, ground truth, GP mean + 95%% band.

clear; clc; close all;

% Ensure this folder (logistic_growth, ground_truth_logistic) and parent (metric_helpers) are on the path
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..'));

%% 1. Ground truth (500 points, [0, 50])
[t_gt, y_gt] = ground_truth_logistic(500);  % 500 x 1 each
fprintf('Ground truth: %d points on [0, 50]. Range y = [%.2f, %.2f]\n', ...
    length(t_gt), min(y_gt), max(y_gt));

%% 2. Design for this experiment: N and Kernel only
N_list = [5, 10, 25, 50];
kernel_list = {'squaredexponential', 'exponential', 'matern32', 'matern52', 'rationalquadratic'};
kernel_labels = {'SqExp', 'Matern 1/2', 'Matern 3/2', 'Matern 5/2', 'RatQuad'};

% 5% noise: sigma = 0.05 * std(ground truth), per design plan
noise_pct = 0.05;
sigma_noise = noise_pct * std(y_gt);

%% 3. Loop: for each N, sample regularly with 5% noise; for each kernel, fit GP and compute all metrics
% Store predictions and observed points for later plotting (GP curves vs ground truth)
results_point = [];
results_prob = [];
results_calib = [];
pred_ymu = cell(length(kernel_list), length(N_list));   % pred_ymu{ik, in}
pred_ystd = cell(length(kernel_list), length(N_list));
t_obs_by_N = cell(1, length(N_list));
y_obs_by_N = cell(1, length(N_list));

rng(42); %setting seed for reproducibility
for in = 1:length(N_list)
    N = N_list(in);
    % Regular sampling: N even intervals including 0 and 50
    t_obs = linspace(0, 50, N)';
    y_obs = interp1(t_gt, y_gt, t_obs, 'linear', 'extrap');
    % Add 5% Gaussian noise
    y_obs = y_obs + sigma_noise * randn(size(y_obs));
    t_obs_by_N{in} = t_obs;
    y_obs_by_N{in} = y_obs;

    for ik = 1:length(kernel_list)
        kern = kernel_list{ik};
        try
            gpr = fitrgp(t_obs, y_obs, ...
                'KernelFunction', kern, ...
                'BasisFunction', 'constant', ...
                'FitMethod', 'exact', ...
                'PredictMethod', 'exact', ...
                'Standardize', true); %Think about if I need/want this

            [ymu, ystd, ~] = predict(gpr, t_gt);
            pred_ymu{ik, in} = ymu;
            pred_ystd{ik, in} = ystd;
            % ystd values returned by predict are the standard deviations of
            % the predicted response distribution at the new points, which
            % accounts for both the mean prediction uncertainty and the
            % observation noise

            % Point-estimate metrics: RMSE, MAE, R²
            [rmse, mae, r2] = metric_helpers('point_estimate', y_gt, ymu);
            results_point = [results_point; {kernel_labels{ik}, N, rmse, mae, r2}];

            % Probabilistic metrics: NLPD, MSLL, CRPS
            [nlpd, msll, crps] = metric_helpers('probabilistic', y_gt, ymu, ystd, y_obs);
            results_prob = [results_prob; {kernel_labels{ik}, N, nlpd, msll, crps}];

            % Calibration & robustness: sMSE, Coverage, NLML
            [smse, coverage, nlml] = metric_helpers('calibration', y_gt, ymu, ystd, y_obs, gpr);
            results_calib = [results_calib; {kernel_labels{ik}, N, smse, coverage, nlml}];
        catch me
            warning('Failed: N=%d, kernel=%s. %s', N, kern, me.message);
            results_point = [results_point; {kernel_labels{ik}, N, NaN, NaN, NaN}];
            results_prob = [results_prob; {kernel_labels{ik}, N, NaN, NaN, NaN}];
            results_calib = [results_calib; {kernel_labels{ik}, N, NaN, NaN, NaN}];
            pred_ymu{ik, in} = nan(size(t_gt));
            pred_ystd{ik, in} = nan(size(t_gt));
        end
    end
end

%% 4. Build and display three metric tables
T_point = cell2table(results_point, 'VariableNames', {'Kernel', 'N', 'RMSE', 'MAE', 'R2'});
T_prob = cell2table(results_prob, 'VariableNames', {'Kernel', 'N', 'NLPD', 'MSLL', 'CRPS'});
T_calib = cell2table(results_calib, 'VariableNames', {'Kernel', 'N', 'sMSE', 'Coverage', 'NLML'});

fprintf('--- Point-Estimate Metrics (Accuracy) ---\n');
disp(T_point);
fprintf('--- Probabilistic Metrics (Uncertainty Quality) ---\n');
disp(T_prob);
fprintf('--- Calibration & Robustness Metrics ---\n');
disp(T_calib);

%% 5. Simple plot: RMSE vs N, one line per kernel
figure;
hold on;
for ik = 1:length(kernel_labels)
    idx = strcmp(T_point.Kernel, kernel_labels{ik});
    plot(T_point.N(idx), T_point.RMSE(idx), '-o', 'LineWidth', 1.5, 'DisplayName', kernel_labels{ik});
end
xlabel('Number of points (N)');
ylabel('RMSE');
title('GP Logistic Exp 03: Regular sampling, 5% noise');
legend('Location', 'best');
grid on;
set(gca, 'XTick', N_list);

%% 6. One figure per kernel: 4 subplots (N = 5, 10, 25, 50), each with points, ground truth, GP mean and 95% band
for ik = 1:length(kernel_labels)
    figure('Name', sprintf('Exp 03: %s', kernel_labels{ik}));
    for in = 1:length(N_list)
        subplot(2, 2, in);
        N = N_list(in);
        t_obs = t_obs_by_N{in};
        y_obs = y_obs_by_N{in};
        ymu = pred_ymu{ik, in};
        ystd = pred_ystd{ik, in};

        % Ground truth curve
        plot(t_gt, y_gt, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Ground truth');
        hold on;
        % Observed points
        plot(t_obs, y_obs, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
        % GP mean
        plot(t_gt, ymu, 'b-', 'LineWidth', 1.2, 'DisplayName', 'GP mean');
        % 95% confidence band (mean +/- 1.96*std)
        if ~all(isnan(ystd))
            ylo = ymu - 1.96 * ystd;
            yhi = ymu + 1.96 * ystd;
            fill([t_gt; flip(t_gt)], [ylo; flip(yhi)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        hold off;
        xlabel('Time');
        ylabel('Population');
        title(sprintf('N = %d', N));
        legend('Location', 'best', 'FontSize', 8);
        grid on;
        xlim([0 50]);
    end
    sgtitle(sprintf('GP Logistic Exp 03: %s — Regular sampling, 5%% noise', kernel_labels{ik}));
end

%% 7. Save (optional)
% save(fullfile(script_dir, 'results_exp_03_regular_5noise.mat'), 'T_point', 'T_prob', 'T_calib', 't_gt', 'y_gt');
% writetable(T_point, fullfile(script_dir, 'results_exp_03_point_metrics.csv'));
% writetable(T_prob, fullfile(script_dir, 'results_exp_03_prob_metrics.csv'));
% writetable(T_calib, fullfile(script_dir, 'results_exp_03_calib_metrics.csv'));

fprintf('Done. Total runs: %d\n', height(T_point));
