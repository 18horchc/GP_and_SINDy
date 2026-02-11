%% GP Logistic Experiment 16: Irregular sampling, 0% noise, Auto hyperparameter optimization
% Same as log_exp_06 but with OptimizeHyperparameters = 'auto'. Saves results_exp_16_irregular_0noise_auto.mat.
clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..'));

[t_gt, y_gt] = ground_truth_logistic(500);
fprintf('Ground truth: %d points on [0, 50]. Range y = [%.2f, %.2f]\n', length(t_gt), min(y_gt), max(y_gt));

N_list = [5, 10, 25, 50];
kernel_list = {'squaredexponential', 'exponential', 'matern32', 'matern52', 'rationalquadratic'};
kernel_labels = {'SqExp', 'Matern 1/2', 'Matern 3/2', 'Matern 5/2', 'RatQuad'};
rng(42);

results_point = []; results_prob = []; results_calib = [];
pred_ymu = cell(length(kernel_list), length(N_list));
pred_ystd = cell(length(kernel_list), length(N_list));
t_obs_by_N = cell(1, length(N_list));
y_obs_by_N = cell(1, length(N_list));

for in = 1:length(N_list)
    N = N_list(in);
    t_obs = sort(randperm(51, N) - 1)';
    y_obs = interp1(t_gt, y_gt, t_obs, 'linear', 'extrap');
    t_obs_by_N{in} = t_obs;
    y_obs_by_N{in} = y_obs;

    for ik = 1:length(kernel_list)
        kern = kernel_list{ik};
        try
            gpr = fitrgp(t_obs, y_obs, ...
                'KernelFunction', kern, 'BasisFunction', 'constant', ...
                'FitMethod', 'exact', 'PredictMethod', 'exact', 'Standardize', true, ...
                'OptimizeHyperparameters', 'auto', ...
                'HyperparameterOptimizationOptions', struct('ShowPlots', false, 'Verbose', 0));
            [ymu, ystd, ~] = predict(gpr, t_gt);
            pred_ymu{ik, in} = ymu;
            pred_ystd{ik, in} = ystd;
            [rmse, mae, r2] = metric_helpers('point_estimate', y_gt, ymu);
            results_point = [results_point; {kernel_labels{ik}, N, rmse, mae, r2}];
            [nlpd, msll, crps] = metric_helpers('probabilistic', y_gt, ymu, ystd, y_obs);
            results_prob = [results_prob; {kernel_labels{ik}, N, nlpd, msll, crps}];
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

T_point = cell2table(results_point, 'VariableNames', {'Kernel', 'N', 'RMSE', 'MAE', 'R2'});
T_prob = cell2table(results_prob, 'VariableNames', {'Kernel', 'N', 'NLPD', 'MSLL', 'CRPS'});
T_calib = cell2table(results_calib, 'VariableNames', {'Kernel', 'N', 'sMSE', 'Coverage', 'NLML'});
fprintf('--- Point-Estimate Metrics ---\n'); disp(T_point);
fprintf('--- Probabilistic Metrics ---\n'); disp(T_prob);
fprintf('--- Calibration Metrics ---\n'); disp(T_calib);

figure;
hold on;
for ik = 1:length(kernel_labels)
    idx = strcmp(T_point.Kernel, kernel_labels{ik});
    plot(T_point.N(idx), T_point.RMSE(idx), '-o', 'LineWidth', 1.5, 'DisplayName', kernel_labels{ik});
end
xlabel('N'); ylabel('RMSE'); title('GP Logistic Exp 16: Irregular, 0% noise, Auto');
legend('Location', 'best'); grid on; set(gca, 'XTick', N_list);

for ik = 1:length(kernel_labels)
    figure('Name', sprintf('Exp 16: %s', kernel_labels{ik}));
    for in = 1:length(N_list)
        subplot(2, 2, in);
        t_obs = t_obs_by_N{in}; y_obs = y_obs_by_N{in};
        ymu = pred_ymu{ik, in}; ystd = pred_ystd{ik, in};
        plot(t_gt, y_gt, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Ground truth'); hold on;
        plot(t_obs, y_obs, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
        plot(t_gt, ymu, 'b-', 'LineWidth', 1.2, 'DisplayName', 'GP mean');
        if ~all(isnan(ystd))
            ylo = ymu - 1.96*ystd; yhi = ymu + 1.96*ystd;
            fill([t_gt; flip(t_gt)], [ylo; flip(yhi)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        hold off; xlabel('Time'); ylabel('Population'); title(sprintf('N = %d', N_list(in)));
        legend('Location', 'best', 'FontSize', 8); grid on; xlim([0 50]);
    end
    sgtitle(sprintf('GP Logistic Exp 16: %s — Irregular, 0%% noise, Auto', kernel_labels{ik}));
end

save(fullfile(script_dir, 'results_exp_16_irregular_0noise_auto.mat'), 'T_point', 'T_prob', 'T_calib', 't_gt', 'y_gt');
writetable(T_point, fullfile(script_dir, 'results_exp_16_point_metrics.csv'));
writetable(T_prob, fullfile(script_dir, 'results_exp_16_prob_metrics.csv'));
writetable(T_calib, fullfile(script_dir, 'results_exp_16_calib_metrics.csv'));
fprintf('Done. Total runs: %d\n', height(T_point));
