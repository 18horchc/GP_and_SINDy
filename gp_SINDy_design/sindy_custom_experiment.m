%% SINDY_CUSTOM_EXPERIMENT  Standalone GP experiment on SINDy [0,14] with custom sampling.
%
%  Completely separate from run_SINDy_experiment. Does NOT save to any metrics files.
%
%  Custom sampling:
%    M1: 2 @ t=0, 2 @ t=1, 2 @ t=2, 5 @ t=3, 1 @ t=5, 4 @ t=7, 3 @ t=14  (19 pts)
%    M2: 2 @ t=0, 3 @ t=1, 2 @ t=2, 5 @ t=3, 1 @ t=5, 4 @ t=7, 3 @ t=14  (20 pts)
%
%  SET NOISE LEVEL HERE (fraction, e.g. 0.05 = 5%, 0 = no noise):
noise_pct = 0.2;   % <-- CHANGE THIS to set noise level

%% Setup
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..'));

%% Ground truth on [0, 14]
params = struct('tspan', [0 14]);
[t_raw, M1_raw, M2_raw, ~] = sindy_mcmc_microglia(params, false);
t_gt = linspace(0, 14, 500)';
M1_gt = interp1(t_raw, M1_raw, t_gt, 'linear', 'extrap');
M2_gt = interp1(t_raw, M2_raw, t_gt, 'linear', 'extrap');

%% Custom sampling (different replicate counts per time point)
% M1: [t, n_replicates]
schedule_M1 = [0 2; 1 2; 2 2; 3 5; 5 1; 7 4; 14 3];
% M2: [t, n_replicates]
schedule_M2 = [0 2; 1 3; 2 2; 3 5; 5 1; 7 4; 14 3];

[t_obs_M1, M1_obs] = build_obs(M1_gt, t_gt, schedule_M1, noise_pct);
[t_obs_M2, M2_obs] = build_obs(M2_gt, t_gt, schedule_M2, noise_pct);

%% GP setup (same 5 kernels as main experiments)
kernel_list = {'squaredexponential', 'exponential', 'matern32', 'matern52', 'rationalquadratic'};
kernel_labels = {'SqExp', 'Matern 1/2', 'Matern 3/2', 'Matern 5/2', 'RatQuad'};
base_opts = {'BasisFunction', 'constant', 'FitMethod', 'exact', ...
    'PredictMethod', 'exact', 'Standardize', true};

fprintf('\n=== SINDy Custom Experiment [0, 14] ===\n');
fprintf('Noise: %.0f%%\n', noise_pct * 100);
fprintf('M1: %d observations | M2: %d observations\n\n', numel(M1_obs), numel(M2_obs));

%% Fit GP and compute metrics for each kernel
pred_M1 = cell(1, length(kernel_list));
pred_M2 = cell(1, length(kernel_list));
% Collect metrics for tables
point_data = [];
prob_data = [];
calib_data = [];
for ik = 1:length(kernel_list)
    kern = kernel_list{ik};
    fit_opts = [{'KernelFunction', kern}, base_opts];
    try
        % M1
        gpr_M1 = fitrgp(t_obs_M1, M1_obs, fit_opts{:});
        [ymu_M1, ystd_M1] = predict(gpr_M1, t_gt);
        [rmse_1, mae_1, r2_1] = metric_helpers('point_estimate', M1_gt, ymu_M1);
        [nlpd_1, msll_1, crps_1] = metric_helpers('probabilistic', M1_gt, ymu_M1, ystd_M1, M1_obs);
        [smse_1, cov_1, nlml_1] = metric_helpers('calibration', M1_gt, ymu_M1, ystd_M1, M1_obs, gpr_M1);

        % M2
        gpr_M2 = fitrgp(t_obs_M2, M2_obs, fit_opts{:});
        [ymu_M2, ystd_M2] = predict(gpr_M2, t_gt);
        [rmse_2, mae_2, r2_2] = metric_helpers('point_estimate', M2_gt, ymu_M2);
        [nlpd_2, msll_2, crps_2] = metric_helpers('probabilistic', M2_gt, ymu_M2, ystd_M2, M2_obs);
        [smse_2, cov_2, nlml_2] = metric_helpers('calibration', M2_gt, ymu_M2, ystd_M2, M2_obs, gpr_M2);

        point_data = [point_data; {kernel_labels{ik}, 'M1', rmse_1, mae_1, r2_1}; {kernel_labels{ik}, 'M2', rmse_2, mae_2, r2_2}];
        prob_data  = [prob_data;  {kernel_labels{ik}, 'M1', nlpd_1, msll_1, crps_1}; {kernel_labels{ik}, 'M2', nlpd_2, msll_2, crps_2}];
        calib_data = [calib_data; {kernel_labels{ik}, 'M1', smse_1, cov_1, nlml_1}; {kernel_labels{ik}, 'M2', smse_2, cov_2, nlml_2}];

        pred_M1{ik} = struct('ymu', ymu_M1, 'ystd', ystd_M1);
        pred_M2{ik} = struct('ymu', ymu_M2, 'ystd', ystd_M2);
    catch me
        warning('Kernel %s failed: %s', kernel_labels{ik}, me.message);
        point_data = [point_data; {kernel_labels{ik}, 'M1', NaN, NaN, NaN}; {kernel_labels{ik}, 'M2', NaN, NaN, NaN}];
        prob_data  = [prob_data;  {kernel_labels{ik}, 'M1', NaN, NaN, NaN}; {kernel_labels{ik}, 'M2', NaN, NaN, NaN}];
        calib_data = [calib_data; {kernel_labels{ik}, 'M1', NaN, NaN, NaN}; {kernel_labels{ik}, 'M2', NaN, NaN, NaN}];
        pred_M1{ik} = struct('ymu', nan(size(t_gt)), 'ystd', nan(size(t_gt)));
        pred_M2{ik} = pred_M1{ik};
    end
end

%% Display metrics in three tables
T_point = cell2table(point_data, 'VariableNames', {'Kernel', 'State', 'RMSE', 'MAE', 'R2'});
T_prob  = cell2table(prob_data,  'VariableNames', {'Kernel', 'State', 'NLPD', 'MSLL', 'CRPS'});
T_calib = cell2table(calib_data, 'VariableNames', {'Kernel', 'State', 'sMSE', 'Coverage', 'NLML'});

fprintf('\n--- Point estimate metrics ---\n');
disp(T_point);
fprintf('\n--- Probabilistic metrics ---\n');
disp(T_prob);
fprintf('\n--- Calibration and robustness metrics ---\n');
disp(T_calib);

%% Plot (one figure per kernel, same layout as main experiments)
for ik = 1:length(kernel_list)
    figure('Name', sprintf('SINDy Custom [0,14]: %s', kernel_labels{ik}));
    subplot(1, 2, 1);
    plot(t_gt, M1_gt, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Ground truth'); hold on;
    plot(t_obs_M1, M1_obs, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
    ymu = pred_M1{ik}.ymu; ystd = pred_M1{ik}.ystd;
    plot(t_gt, ymu, 'b-', 'LineWidth', 1.2, 'DisplayName', 'GP mean');
    if ~all(isnan(ystd))
        ylo = ymu - 1.96*ystd; yhi = ymu + 1.96*ystd;
        fill([t_gt; flip(t_gt)], [ylo; flip(yhi)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    hold off; xlabel('Time'); ylabel('M1'); title(sprintf('M1 — %s', kernel_labels{ik}));
    legend('Location', 'best'); grid on; xlim([0 14]);

    subplot(1, 2, 2);
    plot(t_gt, M2_gt, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Ground truth'); hold on;
    plot(t_obs_M2, M2_obs, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
    ymu = pred_M2{ik}.ymu; ystd = pred_M2{ik}.ystd;
    plot(t_gt, ymu, 'r-', 'LineWidth', 1.2, 'DisplayName', 'GP mean');
    if ~all(isnan(ystd))
        ylo = ymu - 1.96*ystd; yhi = ymu + 1.96*ystd;
        fill([t_gt; flip(t_gt)], [ylo; flip(yhi)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    hold off; xlabel('Time'); ylabel('M2'); title(sprintf('M2 — %s', kernel_labels{ik}));
    legend('Location', 'best'); grid on; xlim([0 14]);

    sgtitle(sprintf('SINDy Custom [0,14], %.0f%% noise — %s', noise_pct*100, kernel_labels{ik}));
end

fprintf('Done. Figures show GP fits for each kernel.\n');

%% Helper: build (t_obs, y_obs) from schedule with optional noise
function [t_obs, y_obs] = build_obs(y_gt, t_gt, schedule, noise_pct)
    t_obs = [];
    y_obs = [];
    sigma_floor = 1e-10 * std(y_gt);
    for row = 1:size(schedule, 1)
        t_val = schedule(row, 1);
        n_rep = schedule(row, 2);
        y_val = interp1(t_gt, y_gt, t_val, 'linear', 'extrap');
        t_obs = [t_obs; repmat(t_val, n_rep, 1)];
        if noise_pct > 0
            sigma_i = max(noise_pct * abs(y_val), sigma_floor);
            y_vals = y_val + sigma_i * randn(n_rep, 1);
        else
            y_vals = repmat(y_val, n_rep, 1);
        end
        y_obs = [y_obs; y_vals];
    end
    t_obs = t_obs(:);
    y_obs = y_obs(:);
end