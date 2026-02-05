%% GP Kinetics Experiment 08: Irregular sampling, 5% noise, varying N and kernel
%
% Same as kin_exp_06 but with 5% Gaussian noise on sampled X, Y, Z.
% Irregular sampling: N distinct whole-number times in [0, 50], one fixed scheme per N (rng(42)).
%
% Scope: IRREGULAR sampling; 5% noise (sigma = 0.05*std per state); N = 5, 10, 25, 50; 5 kernels; one GP per state (X, Y, Z).
% Output: Same tables and figures as kin_exp_06 (titles/saves say 08, 5% noise).

clear; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..'));

[t_gt, x_gt, y_gt, z_gt] = ground_truth_kinetics(500);
fprintf('Ground truth kinetics: %d points on [0, 50]. X [%.4f, %.4f], Y [%.2e, %.2e], Z [%.4f, %.4f]\n', ...
    length(t_gt), min(x_gt), max(x_gt), min(y_gt), max(y_gt), min(z_gt), max(z_gt));

N_list = [5, 10, 25, 50];
kernel_list = {'squaredexponential', 'exponential', 'matern32', 'matern52', 'rationalquadratic'};
kernel_labels = {'SqExp', 'Matern 1/2', 'Matern 3/2', 'Matern 5/2', 'RatQuad'};
state_labels = {'X', 'Y', 'Z'};

noise_pct = 0.05;
sigma_noise_x = noise_pct * std(x_gt);
sigma_noise_y = noise_pct * std(y_gt);
sigma_noise_z = noise_pct * std(z_gt);
rng(42);

results_point = []; results_prob = []; results_calib = [];
pred_ymu_x = cell(length(kernel_list), length(N_list)); pred_ystd_x = pred_ymu_x;
pred_ymu_y = pred_ymu_x; pred_ystd_y = pred_ymu_x; pred_ymu_z = pred_ymu_x; pred_ystd_z = pred_ymu_x;
t_obs_by_N = cell(1, length(N_list)); x_obs_by_N = t_obs_by_N; y_obs_by_N = t_obs_by_N; z_obs_by_N = t_obs_by_N;

for in = 1:length(N_list)
    N = N_list(in);
    t_obs = sort(randperm(51, N) - 1)';
    x_obs = interp1(t_gt, x_gt, t_obs, 'linear', 'extrap');
    y_obs = interp1(t_gt, y_gt, t_obs, 'linear', 'extrap');
    z_obs = interp1(t_gt, z_gt, t_obs, 'linear', 'extrap');
    x_obs = x_obs + sigma_noise_x * randn(size(x_obs));
    y_obs = y_obs + sigma_noise_y * randn(size(y_obs));
    z_obs = z_obs + sigma_noise_z * randn(size(z_obs));
    t_obs_by_N{in} = t_obs; x_obs_by_N{in} = x_obs; y_obs_by_N{in} = y_obs; z_obs_by_N{in} = z_obs;

    for ik = 1:length(kernel_list)
        kern = kernel_list{ik};
        try
            gpr_x = fitrgp(t_obs, x_obs, 'KernelFunction', kern, 'BasisFunction', 'constant', 'FitMethod', 'exact', 'PredictMethod', 'exact', 'Standardize', true);
            [ymu_x, ystd_x, ~] = predict(gpr_x, t_gt);
            pred_ymu_x{ik, in} = ymu_x; pred_ystd_x{ik, in} = ystd_x;
            [rmse_x, mae_x, r2_x] = metric_helpers('point_estimate', x_gt, ymu_x);
            [nlpd_x, msll_x, crps_x] = metric_helpers('probabilistic', x_gt, ymu_x, ystd_x, x_obs);
            [smse_x, cov_x, nlml_x] = metric_helpers('calibration', x_gt, ymu_x, ystd_x, x_obs, gpr_x);
            results_point = [results_point; {kernel_labels{ik}, N, 'X', rmse_x, mae_x, r2_x}];
            results_prob  = [results_prob;  {kernel_labels{ik}, N, 'X', nlpd_x, msll_x, crps_x}];
            results_calib = [results_calib; {kernel_labels{ik}, N, 'X', smse_x, cov_x, nlml_x}];

            gpr_y = fitrgp(t_obs, y_obs, 'KernelFunction', kern, 'BasisFunction', 'constant', 'FitMethod', 'exact', 'PredictMethod', 'exact', 'Standardize', true);
            [ymu_y, ystd_y, ~] = predict(gpr_y, t_gt);
            pred_ymu_y{ik, in} = ymu_y; pred_ystd_y{ik, in} = ystd_y;
            [rmse_y, mae_y, r2_y] = metric_helpers('point_estimate', y_gt, ymu_y);
            [nlpd_y, msll_y, crps_y] = metric_helpers('probabilistic', y_gt, ymu_y, ystd_y, y_obs);
            [smse_y, cov_y, nlml_y] = metric_helpers('calibration', y_gt, ymu_y, ystd_y, y_obs, gpr_y);
            results_point = [results_point; {kernel_labels{ik}, N, 'Y', rmse_y, mae_y, r2_y}];
            results_prob  = [results_prob;  {kernel_labels{ik}, N, 'Y', nlpd_y, msll_y, crps_y}];
            results_calib = [results_calib; {kernel_labels{ik}, N, 'Y', smse_y, cov_y, nlml_y}];

            gpr_z = fitrgp(t_obs, z_obs, 'KernelFunction', kern, 'BasisFunction', 'constant', 'FitMethod', 'exact', 'PredictMethod', 'exact', 'Standardize', true);
            [ymu_z, ystd_z, ~] = predict(gpr_z, t_gt);
            pred_ymu_z{ik, in} = ymu_z; pred_ystd_z{ik, in} = ystd_z;
            [rmse_z, mae_z, r2_z] = metric_helpers('point_estimate', z_gt, ymu_z);
            [nlpd_z, msll_z, crps_z] = metric_helpers('probabilistic', z_gt, ymu_z, ystd_z, z_obs);
            [smse_z, cov_z, nlml_z] = metric_helpers('calibration', z_gt, ymu_z, ystd_z, z_obs, gpr_z);
            results_point = [results_point; {kernel_labels{ik}, N, 'Z', rmse_z, mae_z, r2_z}];
            results_prob  = [results_prob;  {kernel_labels{ik}, N, 'Z', nlpd_z, msll_z, crps_z}];
            results_calib = [results_calib; {kernel_labels{ik}, N, 'Z', smse_z, cov_z, nlml_z}];
        catch me
            warning('Failed: N=%d, kernel=%s. %s', N, kernel_labels{ik}, me.message);
            for s = 1:3
                st = state_labels{s};
                results_point = [results_point; {kernel_labels{ik}, N, st, NaN, NaN, NaN}];
                results_prob  = [results_prob;  {kernel_labels{ik}, N, st, NaN, NaN, NaN}];
                results_calib = [results_calib; {kernel_labels{ik}, N, st, NaN, NaN, NaN}];
            end
            pred_ymu_x{ik, in} = nan(size(t_gt)); pred_ystd_x{ik, in} = nan(size(t_gt));
            pred_ymu_y{ik, in} = nan(size(t_gt)); pred_ystd_y{ik, in} = nan(size(t_gt));
            pred_ymu_z{ik, in} = nan(size(t_gt)); pred_ystd_z{ik, in} = nan(size(t_gt));
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
for is = 1:3
    subplot(1, 3, is); hold on;
    for ik = 1:length(kernel_labels)
        idx = strcmp(T_point.Kernel, kernel_labels{ik}) & strcmp(T_point.State, state_labels{is});
        plot(T_point.N(idx), T_point.RMSE(idx), '-o', 'LineWidth', 1.5, 'DisplayName', kernel_labels{ik});
    end
    xlabel('N'); ylabel('RMSE'); title(sprintf('Kin Exp 08: %s — Irregular, 5%% noise', state_labels{is}));
    legend('Location', 'best'); grid on; set(gca, 'XTick', N_list);
end
sgtitle('GP Kinetics Exp 08: Irregular sampling, 5% noise');

gt_cell = {x_gt, y_gt, z_gt}; obs_cell = {x_obs_by_N, y_obs_by_N, z_obs_by_N};
ymu_cell = {pred_ymu_x, pred_ymu_y, pred_ymu_z}; ystd_cell = {pred_ystd_x, pred_ystd_y, pred_ystd_z};
colors = {'b', 'g', 'r'};
for ik = 1:length(kernel_labels)
    figure('Name', sprintf('Kin Exp 08: %s', kernel_labels{ik}));
    for is = 1:3
        y_gt_s = gt_cell{is};
        for in = 1:length(N_list)
            subplot(3, 4, (is-1)*4 + in);
            t_obs = t_obs_by_N{in}; y_obs_s = obs_cell{is}{in}; ymu = ymu_cell{is}{ik, in}; ystd = ystd_cell{is}{ik, in};
            plot(t_gt, y_gt_s, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Ground truth'); hold on;
            plot(t_obs, y_obs_s, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
            plot(t_gt, ymu, '-', 'Color', colors{is}, 'LineWidth', 1.2, 'DisplayName', 'GP mean');
            if ~all(isnan(ystd))
                ylo = ymu - 1.96*ystd; yhi = ymu + 1.96*ystd;
                fill([t_gt; flip(t_gt)], [ylo; flip(yhi)], colors{is}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            end
            hold off; xlabel('Time'); ylabel(state_labels{is}); title(sprintf('%s, N = %d', state_labels{is}, N_list(in)));
            legend('Location', 'best', 'FontSize', 8); grid on; xlim([0 50]);
        end
    end
    sgtitle(sprintf('GP Kinetics Exp 08: %s — Irregular, 5%% noise', kernel_labels{ik}));
end

save(fullfile(script_dir, 'results_kin_exp_08_irregular_5noise.mat'), 'T_point', 'T_prob', 'T_calib', 't_gt', 'x_gt', 'y_gt', 'z_gt');
writetable(T_point, fullfile(script_dir, 'results_kin_exp_08_point_metrics.csv'));
writetable(T_prob, fullfile(script_dir, 'results_kin_exp_08_prob_metrics.csv'));
writetable(T_calib, fullfile(script_dir, 'results_kin_exp_08_calib_metrics.csv'));
fprintf('Done. Total rows: %d (X, Y, Z per Kernel,N)\n', height(T_point));
