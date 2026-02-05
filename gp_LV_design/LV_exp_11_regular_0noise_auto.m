%% GP Lotka-Volterra Experiment 11: Regular sampling, 0% noise, Auto hyperparameter optimization
% Same as LV_exp_01 but with OptimizeHyperparameters = 'auto' for built-in kernels (Periodic unchanged). Saves results_exp_11_regular_0noise_auto.mat.
clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..'));

[t_gt, prey_gt, pred_gt] = ground_truth_LV(500);
fprintf('Ground truth LV: %d points on [0, 50]. Prey [%.2f, %.2f], Predator [%.2f, %.2f]\n', ...
    length(t_gt), min(prey_gt), max(prey_gt), min(pred_gt), max(pred_gt));

N_list = [5, 10, 25, 50];
kernel_list = {'squaredexponential', 'exponential', 'matern32', 'matern52', 'rationalquadratic', @periodicKernel};
kernel_labels = {'SqExp', 'Matern 1/2', 'Matern 3/2', 'Matern 5/2', 'RatQuad', 'Periodic'};
n_anticipated_cycles = 3;

results_point = []; results_prob = []; results_calib = [];
pred_ymu_prey  = cell(length(kernel_list), length(N_list));
pred_ystd_prey = cell(length(kernel_list), length(N_list));
pred_ymu_pred  = cell(length(kernel_list), length(N_list));
pred_ystd_pred = cell(length(kernel_list), length(N_list));
t_obs_by_N   = cell(1, length(N_list));
prey_obs_by_N = cell(1, length(N_list));
pred_obs_by_N = cell(1, length(N_list));

hopts = struct('ShowPlots', false, 'Verbose', 0);

for in = 1:length(N_list)
    N = N_list(in);
    t_obs = linspace(0, 50, N)';
    prey_obs = interp1(t_gt, prey_gt, t_obs, 'linear', 'extrap');
    pred_obs = interp1(t_gt, pred_gt, t_obs, 'linear', 'extrap');
    t_obs_by_N{in} = t_obs;
    prey_obs_by_N{in} = prey_obs;
    pred_obs_by_N{in} = pred_obs;

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
            if is_periodic
                gpr_prey = fitrgp(t_obs, prey_obs, ...
                    'KernelFunction', @periodicKernel, 'KernelParameters', theta0_prey, ...
                    'BasisFunction', 'constant', 'FitMethod', 'exact', 'PredictMethod', 'exact', 'Standardize', true);
            else
                gpr_prey = fitrgp(t_obs, prey_obs, ...
                    'KernelFunction', kern, 'BasisFunction', 'constant', ...
                    'FitMethod', 'exact', 'PredictMethod', 'exact', 'Standardize', true, ...
                    'OptimizeHyperparameters', 'auto', 'HyperparameterOptimizationOptions', hopts);
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

            if is_periodic
                gpr_pred = fitrgp(t_obs, pred_obs, ...
                    'KernelFunction', @periodicKernel, 'KernelParameters', theta0_pred, ...
                    'BasisFunction', 'constant', 'FitMethod', 'exact', 'PredictMethod', 'exact', 'Standardize', true);
            else
                gpr_pred = fitrgp(t_obs, pred_obs, ...
                    'KernelFunction', kern, 'BasisFunction', 'constant', ...
                    'FitMethod', 'exact', 'PredictMethod', 'exact', 'Standardize', true, ...
                    'OptimizeHyperparameters', 'auto', 'HyperparameterOptimizationOptions', hopts);
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

T_point = cell2table(results_point, 'VariableNames', {'Kernel', 'N', 'State', 'RMSE', 'MAE', 'R2'});
T_prob  = cell2table(results_prob,  'VariableNames', {'Kernel', 'N', 'State', 'NLPD', 'MSLL', 'CRPS'});
T_calib = cell2table(results_calib, 'VariableNames', {'Kernel', 'N', 'State', 'sMSE', 'Coverage', 'NLML'});
fprintf('--- Point-Estimate Metrics (Accuracy) ---\n'); disp(T_point);
fprintf('--- Probabilistic Metrics (Uncertainty Quality) ---\n'); disp(T_prob);
fprintf('--- Calibration & Robustness Metrics ---\n'); disp(T_calib);

figure;
subplot(1,2,1);
hold on;
for ik = 1:length(kernel_labels)
    idx = strcmp(T_point.Kernel, kernel_labels{ik}) & strcmp(T_point.State, 'Prey');
    plot(T_point.N(idx), T_point.RMSE(idx), '-o', 'LineWidth', 1.5, 'DisplayName', kernel_labels{ik});
end
xlabel('N'); ylabel('RMSE'); title('LV Exp 11: Prey — Regular, 0% noise, Auto');
legend('Location', 'best'); grid on; set(gca, 'XTick', N_list);
subplot(1,2,2);
hold on;
for ik = 1:length(kernel_labels)
    idx = strcmp(T_point.Kernel, kernel_labels{ik}) & strcmp(T_point.State, 'Predator');
    plot(T_point.N(idx), T_point.RMSE(idx), '-o', 'LineWidth', 1.5, 'DisplayName', kernel_labels{ik});
end
xlabel('N'); ylabel('RMSE'); title('LV Exp 11: Predator — Regular, 0% noise, Auto');
legend('Location', 'best'); grid on; set(gca, 'XTick', N_list);
sgtitle('GP LV Exp 11: Regular sampling, 0% noise, Auto');

for ik = 1:length(kernel_labels)
    figure('Name', sprintf('LV Exp 11: %s', kernel_labels{ik}));
    for in = 1:length(N_list)
        subplot(2, 4, in);
        t_obs = t_obs_by_N{in}; prey_obs = prey_obs_by_N{in};
        ymu = pred_ymu_prey{ik, in}; ystd = pred_ystd_prey{ik, in};
        plot(t_gt, prey_gt, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Ground truth'); hold on;
        plot(t_obs, prey_obs, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
        plot(t_gt, ymu, 'b-', 'LineWidth', 1.2, 'DisplayName', 'GP mean');
        if ~all(isnan(ystd))
            ylo = ymu - 1.96*ystd; yhi = ymu + 1.96*ystd;
            fill([t_gt; flip(t_gt)], [ylo; flip(yhi)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        hold off; xlabel('Time'); ylabel('Prey'); title(sprintf('Prey, N = %d', N_list(in)));
        legend('Location', 'best', 'FontSize', 8); grid on; xlim([0 50]);
        subplot(2, 4, 4 + in);
        pred_obs = pred_obs_by_N{in}; ymu = pred_ymu_pred{ik, in}; ystd = pred_ystd_pred{ik, in};
        plot(t_gt, pred_gt, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Ground truth'); hold on;
        plot(t_obs, pred_obs, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
        plot(t_gt, ymu, 'r-', 'LineWidth', 1.2, 'DisplayName', 'GP mean');
        if ~all(isnan(ystd))
            ylo = ymu - 1.96*ystd; yhi = ymu + 1.96*ystd;
            fill([t_gt; flip(t_gt)], [ylo; flip(yhi)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        hold off; xlabel('Time'); ylabel('Predator'); title(sprintf('Predator, N = %d', N_list(in)));
        legend('Location', 'best', 'FontSize', 8); grid on; xlim([0 50]);
    end
    sgtitle(sprintf('GP LV Exp 11: %s — Regular, 0%% noise, Auto', kernel_labels{ik}));
end

save(fullfile(script_dir, 'results_exp_11_regular_0noise_auto.mat'), 'T_point', 'T_prob', 'T_calib', 't_gt', 'prey_gt', 'pred_gt');
writetable(T_point, fullfile(script_dir, 'results_exp_11_point_metrics.csv'));
writetable(T_prob, fullfile(script_dir, 'results_exp_11_prob_metrics.csv'));
writetable(T_calib, fullfile(script_dir, 'results_exp_11_calib_metrics.csv'));
fprintf('Done. Total rows: %d (Prey + Predator per Kernel,N)\n', height(T_point));
