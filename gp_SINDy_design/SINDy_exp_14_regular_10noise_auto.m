%% GP SINDy Experiment 14: Regular sampling, 10% noise, Auto hyperparameter optimization
% Same as SINDy_exp_04 but with OptimizeHyperparameters = 'auto'. Saves results_SINDy_exp_14_regular_110noise_auto.mat.
clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..'));

[t_gt, M1_gt, M2_gt] = ground_truth_SINDy(500);
fprintf('Ground truth SINDy: %d points on [0, 50]. M1 [%.2f, %.2f], M2 [%.2f, %.2f]\n', ...
    length(t_gt), min(M1_gt), max(M1_gt), min(M2_gt), max(M2_gt));

N_list = [5, 10, 25, 50];
kernel_list = {'squaredexponential', 'exponential', 'matern32', 'matern52', 'rationalquadratic'};
kernel_labels = {'SqExp', 'Matern 1/2', 'Matern 3/2', 'Matern 5/2', 'RatQuad'};
state_labels = {'M1', 'M2'};
hopts = struct('ShowPlots', false, 'Verbose', 0);

results_point = []; results_prob = []; results_calib = [];
pred_ymu_M1 = cell(length(kernel_list), length(N_list));
pred_ystd_M1 = cell(length(kernel_list), length(N_list));
pred_ymu_M2 = cell(length(kernel_list), length(N_list));
pred_ystd_M2 = cell(length(kernel_list), length(N_list));
t_obs_by_N   = cell(1, length(N_list));
M1_obs_by_N  = cell(1, length(N_list));
M2_obs_by_N  = cell(1, length(N_list));

for in = 1:length(N_list)
    N = N_list(in);
    t_obs = linspace(0, 50, N)';
    M1_obs = interp1(t_gt, M1_gt, t_obs, 'linear', 'extrap');
    M2_obs = interp1(t_gt, M2_gt, t_obs, 'linear', 'extrap');
    M1_obs = M1_obs + sigma_noise_M1 * randn(size(M1_obs));
    M2_obs = M2_obs + sigma_noise_M2 * randn(size(M2_obs));
    t_obs_by_N{in} = t_obs; M1_obs_by_N{in} = M1_obs; M2_obs_by_N{in} = M2_obs;

    for ik = 1:length(kernel_list)
        kern = kernel_list{ik};
        try
            gpr_M1 = fitrgp(t_obs, M1_obs, 'KernelFunction', kern, 'BasisFunction', 'constant', ...
                'FitMethod', 'exact', 'PredictMethod', 'exact', 'Standardize', true, ...
                'OptimizeHyperparameters', 'auto', 'HyperparameterOptimizationOptions', hopts);
            [ymu_M1, ystd_M1, ~] = predict(gpr_M1, t_gt);
            pred_ymu_M1{ik, in} = ymu_M1; pred_ystd_M1{ik, in} = ystd_M1;
            [rmse_1, mae_1, r2_1] = metric_helpers('point_estimate', M1_gt, ymu_M1);
            [nlpd_1, msll_1, crps_1] = metric_helpers('probabilistic', M1_gt, ymu_M1, ystd_M1, M1_obs);
            [smse_1, cov_1, nlml_1] = metric_helpers('calibration', M1_gt, ymu_M1, ystd_M1, M1_obs, gpr_M1);
            results_point = [results_point; {kernel_labels{ik}, N, 'M1', rmse_1, mae_1, r2_1}];
            results_prob  = [results_prob;  {kernel_labels{ik}, N, 'M1', nlpd_1, msll_1, crps_1}];
            results_calib = [results_calib; {kernel_labels{ik}, N, 'M1', smse_1, cov_1, nlml_1}];

            gpr_M2 = fitrgp(t_obs, M2_obs, 'KernelFunction', kern, 'BasisFunction', 'constant', ...
                'FitMethod', 'exact', 'PredictMethod', 'exact', 'Standardize', true, ...
                'OptimizeHyperparameters', 'auto', 'HyperparameterOptimizationOptions', hopts);
            [ymu_M2, ystd_M2, ~] = predict(gpr_M2, t_gt);
            pred_ymu_M2{ik, in} = ymu_M2; pred_ystd_M2{ik, in} = ystd_M2;
            [rmse_2, mae_2, r2_2] = metric_helpers('point_estimate', M2_gt, ymu_M2);
            [nlpd_2, msll_2, crps_2] = metric_helpers('probabilistic', M2_gt, ymu_M2, ystd_M2, M2_obs);
            [smse_2, cov_2, nlml_2] = metric_helpers('calibration', M2_gt, ymu_M2, ystd_M2, M2_obs, gpr_M2);
            results_point = [results_point; {kernel_labels{ik}, N, 'M2', rmse_2, mae_2, r2_2}];
            results_prob  = [results_prob;  {kernel_labels{ik}, N, 'M2', nlpd_2, msll_2, crps_2}];
            results_calib = [results_calib; {kernel_labels{ik}, N, 'M2', smse_2, cov_2, nlml_2}];
        catch me
            warning('Failed: N=%d, kernel=%s. %s', N, kernel_labels{ik}, me.message);
            results_point = [results_point; {kernel_labels{ik}, N, 'M1', NaN, NaN, NaN}; {kernel_labels{ik}, N, 'M2', NaN, NaN, NaN}];
            results_prob  = [results_prob;  {kernel_labels{ik}, N, 'M1', NaN, NaN, NaN}; {kernel_labels{ik}, N, 'M2', NaN, NaN, NaN}];
            results_calib = [results_calib; {kernel_labels{ik}, N, 'M1', NaN, NaN, NaN}; {kernel_labels{ik}, N, 'M2', NaN, NaN, NaN}];
            pred_ymu_M1{ik, in} = nan(size(t_gt)); pred_ystd_M1{ik, in} = nan(size(t_gt));
            pred_ymu_M2{ik, in} = nan(size(t_gt)); pred_ystd_M2{ik, in} = nan(size(t_gt));
        end
    end
end

T_point = cell2table(results_point, 'VariableNames', {'Kernel', 'N', 'State', 'RMSE', 'MAE', 'R2'});
T_prob  = cell2table(results_prob,  'VariableNames', {'Kernel', 'N', 'State', 'NLPD', 'MSLL', 'CRPS'});
T_calib = cell2table(results_calib, 'VariableNames', {'Kernel', 'N', 'State', 'sMSE', 'Coverage', 'NLML'});
fprintf('--- Point-Estimate Metrics (Accuracy) ---\n'); disp(T_point);
fprintf('--- Probabilistic Metrics (Uncertainty Quality) ---\n'); disp(T_prob);
fprintf('--- Calibration & Robustness Metrics ---\n'); disp(T_calib);

figure;
subplot(1, 2, 1); hold on;
for ik = 1:length(kernel_labels)
    idx = strcmp(T_point.Kernel, kernel_labels{ik}) & strcmp(T_point.State, 'M1');
    plot(T_point.N(idx), T_point.RMSE(idx), '-o', 'LineWidth', 1.5, 'DisplayName', kernel_labels{ik});
end
xlabel('N'); ylabel('RMSE'); title('SINDy Exp 11: M1 — Regular, 10% noise, Auto');
legend('Location', 'best'); grid on; set(gca, 'XTick', N_list);
subplot(1, 2, 2); hold on;
for ik = 1:length(kernel_labels)
    idx = strcmp(T_point.Kernel, kernel_labels{ik}) & strcmp(T_point.State, 'M2');
    plot(T_point.N(idx), T_point.RMSE(idx), '-o', 'LineWidth', 1.5, 'DisplayName', kernel_labels{ik});
end
xlabel('N'); ylabel('RMSE'); title('SINDy Exp 11: M2 — Regular, 10% noise, Auto');
legend('Location', 'best'); grid on; set(gca, 'XTick', N_list);
sgtitle('GP SINDy Exp 11: Regular sampling, 0% noise, Auto');

for ik = 1:length(kernel_labels)
    figure('Name', sprintf('SINDy Exp 11: %s', kernel_labels{ik}));
    for in = 1:length(N_list)
        subplot(2, 4, in);
        t_obs = t_obs_by_N{in}; M1_obs = M1_obs_by_N{in};
        ymu = pred_ymu_M1{ik, in}; ystd = pred_ystd_M1{ik, in};
        plot(t_gt, M1_gt, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Ground truth'); hold on;
        plot(t_obs, M1_obs, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
        plot(t_gt, ymu, 'b-', 'LineWidth', 1.2, 'DisplayName', 'GP mean');
        if ~all(isnan(ystd))
            ylo = ymu - 1.96*ystd; yhi = ymu + 1.96*ystd;
            fill([t_gt; flip(t_gt)], [ylo; flip(yhi)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        hold off; xlabel('Time'); ylabel('M1 (cells/mm^2)'); title(sprintf('M1, N = %d', N_list(in)));
        legend('Location', 'best', 'FontSize', 8); grid on; xlim([0 50]);
        subplot(2, 4, 4 + in);
        M2_obs = M2_obs_by_N{in}; ymu = pred_ymu_M2{ik, in}; ystd = pred_ystd_M2{ik, in};
        plot(t_gt, M2_gt, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Ground truth'); hold on;
        plot(t_obs, M2_obs, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
        plot(t_gt, ymu, 'r-', 'LineWidth', 1.2, 'DisplayName', 'GP mean');
        if ~all(isnan(ystd))
            ylo = ymu - 1.96*ystd; yhi = ymu + 1.96*ystd;
            fill([t_gt; flip(t_gt)], [ylo; flip(yhi)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        hold off; xlabel('Time'); ylabel('M2 (cells/mm^2)'); title(sprintf('M2, N = %d', N_list(in)));
        legend('Location', 'best', 'FontSize', 8); grid on; xlim([0 50]);
    end
    sgtitle(sprintf('GP SINDy Exp 11: %s — Regular, 0%% noise, Auto', kernel_labels{ik}));
end

save(fullfile(script_dir, 'results_SINDy_exp_14_regular_10noise_auto.mat'), 'T_point', 'T_prob', 'T_calib', 't_gt', 'M1_gt', 'M2_gt');
writetable(T_point, fullfile(script_dir, 'results_SINDy_exp_11_point_metrics.csv'));
writetable(T_prob, fullfile(script_dir, 'results_SINDy_exp_11_prob_metrics.csv'));
writetable(T_calib, fullfile(script_dir, 'results_SINDy_exp_11_calib_metrics.csv'));
fprintf('Done. Total rows: %d (M1 + M2 per Kernel,N)\n', height(T_point));
