%% SINDY_CUSTOM_SIMULATION  1000-run simulation of custom experiment at 20% noise.
%
%  Runs the SINDy custom experiment 1000 times, accumulates metrics, and
%  displays the averaged metrics in three tables (same format as sindy_custom_experiment).
%  No plotting.
%
%  Custom sampling (same as sindy_custom_experiment):
%    M1: 2 @ t=0, 2 @ t=1, 2 @ t=2, 5 @ t=3, 1 @ t=5, 4 @ t=7, 3 @ t=14  (19 pts)
%    M2: 2 @ t=0, 3 @ t=1, 2 @ t=2, 5 @ t=3, 1 @ t=5, 4 @ t=7, 3 @ t=14  (20 pts)

noise_pct = 0.20;   % 20% noise (fixed for this simulation)
n_runs = 1000;

%% Setup
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..'));

%% Ground truth on [0, 14] (same each run)
params = struct('tspan', [0 14]);
[t_raw, M1_raw, M2_raw, ~] = sindy_mcmc_microglia(params, false);
t_gt = linspace(0, 14, 500)';
M1_gt = interp1(t_raw, M1_raw, t_gt, 'linear', 'extrap');
M2_gt = interp1(t_raw, M2_raw, t_gt, 'linear', 'extrap');

%% Custom sampling
schedule_M1 = [0 2; 1 2; 2 2; 3 5; 5 1; 7 4; 14 3];
schedule_M2 = [0 2; 1 3; 2 2; 3 5; 5 1; 7 4; 14 3];

%% GP setup
kernel_list = {'squaredexponential', 'exponential', 'matern32', 'matern52', 'rationalquadratic'};
kernel_labels = {'SqExp', 'Matern 1/2', 'Matern 3/2', 'Matern 5/2', 'RatQuad'};
base_opts = {'BasisFunction', 'constant', 'FitMethod', 'exact', ...
    'PredictMethod', 'exact', 'Standardize', true};

%% Accumulators for averaged metrics (10 rows: 5 kernels x 2 states)
n_rows = length(kernel_list) * 2;  % 10
point_sum = zeros(n_rows, 3);   % RMSE, MAE, R2
prob_sum  = zeros(n_rows, 3);   % NLPD, MSLL, CRPS
calib_sum = zeros(n_rows, 3);   % sMSE, Coverage, NLML
row_count = zeros(n_rows, 1);   % count successful runs per row (for averaging if some fail)

fprintf('\n=== SINDy Custom Simulation [0, 14] ===\n');
fprintf('Noise: %.0f%% | Runs: %d\n', noise_pct * 100, n_runs);
fprintf('Running...\n');

for ir = 1:n_runs
    if mod(ir, 100) == 0
        fprintf('  Run %d/%d\n', ir, n_runs);
    end

    [t_obs_M1, M1_obs] = build_obs(M1_gt, t_gt, schedule_M1, noise_pct);
    [t_obs_M2, M2_obs] = build_obs(M2_gt, t_gt, schedule_M2, noise_pct);

    for ik = 1:length(kernel_list)
        kern = kernel_list{ik};
        fit_opts = [{'KernelFunction', kern}, base_opts];
        row_m1 = (ik - 1) * 2 + 1;
        row_m2 = (ik - 1) * 2 + 2;
        try
            % M1
            gpr_M1 = fitrgp(t_obs_M1, M1_obs, fit_opts{:});
            [ymu_M1, ystd_M1] = predict(gpr_M1, t_gt);
            [rmse_1, mae_1, r2_1] = metric_helpers('point_estimate', M1_gt, ymu_M1);
            [nlpd_1, msll_1, crps_1] = metric_helpers('probabilistic', M1_gt, ymu_M1, ystd_M1, M1_obs);
            [smse_1, cov_1, nlml_1] = metric_helpers('calibration', M1_gt, ymu_M1, ystd_M1, M1_obs, gpr_M1);

            point_sum(row_m1, :) = point_sum(row_m1, :) + [rmse_1, mae_1, r2_1];
            prob_sum(row_m1, :)  = prob_sum(row_m1, :)  + [nlpd_1, msll_1, crps_1];
            calib_sum(row_m1,:) = calib_sum(row_m1,:)  + [smse_1, cov_1, nlml_1];
            row_count(row_m1) = row_count(row_m1) + 1;

            % M2
            gpr_M2 = fitrgp(t_obs_M2, M2_obs, fit_opts{:});
            [ymu_M2, ystd_M2] = predict(gpr_M2, t_gt);
            [rmse_2, mae_2, r2_2] = metric_helpers('point_estimate', M2_gt, ymu_M2);
            [nlpd_2, msll_2, crps_2] = metric_helpers('probabilistic', M2_gt, ymu_M2, ystd_M2, M2_obs);
            [smse_2, cov_2, nlml_2] = metric_helpers('calibration', M2_gt, ymu_M2, ystd_M2, M2_obs, gpr_M2);

            point_sum(row_m2, :) = point_sum(row_m2, :) + [rmse_2, mae_2, r2_2];
            prob_sum(row_m2, :)  = prob_sum(row_m2, :)  + [nlpd_2, msll_2, crps_2];
            calib_sum(row_m2,:)  = calib_sum(row_m2,:)  + [smse_2, cov_2, nlml_2];
            row_count(row_m2) = row_count(row_m2) + 1;
        catch
            % Skip failed runs for this kernel-state
        end
    end
end

%% Average (handle possible failures)
for r = 1:n_rows
    if row_count(r) > 0
        point_sum(r, :) = point_sum(r, :) / row_count(r);
        prob_sum(r, :)  = prob_sum(r, :)  / row_count(r);
        calib_sum(r,:) = calib_sum(r,:)  / row_count(r);
    else
        point_sum(r, :) = NaN;
        prob_sum(r, :)  = NaN;
        calib_sum(r,:)  = NaN;
    end
end

%% Build tables (same structure as sindy_custom_experiment)
kern_col = repelem(kernel_labels(:), 2);
state_col = repmat({'M1'; 'M2'}, length(kernel_list), 1);

T_point = table(kern_col, state_col, point_sum(:,1), point_sum(:,2), point_sum(:,3), ...
    'VariableNames', {'Kernel', 'State', 'RMSE', 'MAE', 'R2'});
T_prob  = table(kern_col, state_col, prob_sum(:,1), prob_sum(:,2), prob_sum(:,3), ...
    'VariableNames', {'Kernel', 'State', 'NLPD', 'MSLL', 'CRPS'});
T_calib = table(kern_col, state_col, calib_sum(:,1), calib_sum(:,2), calib_sum(:,3), ...
    'VariableNames', {'Kernel', 'State', 'sMSE', 'Coverage', 'NLML'});

%% Display
fprintf('\n--- Point estimate metrics (mean over %d runs) ---\n', n_runs);
disp(T_point);
fprintf('\n--- Probabilistic metrics (mean over %d runs) ---\n', n_runs);
disp(T_prob);
fprintf('\n--- Calibration and robustness metrics (mean over %d runs) ---\n', n_runs);
disp(T_calib);
fprintf('Done.\n');

%% Helper
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
