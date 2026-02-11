function [T_point, T_prob, T_calib] = run_LV_experiment(cfg, opts)
%% RUN_LV_EXPERIMENT Run one GP Lotka-Volterra experiment with given configuration.
%
%   [T_point, T_prob, T_calib] = run_LV_experiment(cfg)
%   [T_point, T_prob, T_calib] = run_LV_experiment(cfg, opts)
%
%   cfg: struct with fields:
%        .exp_id     - experiment number (1-20), used for output filenames
%        .noise_pct  - noise fraction (0, 0.01, 0.05, 0.10, 0.20)
%        .is_regular - true = linspace sampling, false = randperm sampling
%        .is_auto    - true = OptimizeHyperparameters='auto' for built-in kernels
%                      (Periodic kernel unchanged; uses KernelParameters)
%
%   opts (optional): struct with fields:
%        .out_dir    - where to save (default: folder containing this file)
%        .make_plots - true/false (default true)
%        .rng_seed   - for reproducibility (default 42)
%
%   Output: T_point, T_prob, T_calib (tables with Kernel, N, State). Saves .mat
%   and one merged .csv with all metrics per experiment. Two states: Prey,
%   Predator. Six kernels: 5 built-in + custom periodicKernel.

    if nargin < 2, opts = struct(); end
    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    addpath(fullfile(script_dir, '..'));

    out_dir = get_opt(opts, 'out_dir', script_dir);
    make_plots = get_opt(opts, 'make_plots', true);
    rng_seed = get_opt(opts, 'rng_seed', 42);

    rng(rng_seed);

    %% Ground truth (two states: prey, predator)
    [t_gt, prey_gt, pred_gt] = ground_truth_LV(500);
    sigma_noise_prey = cfg.noise_pct * std(prey_gt);
    sigma_noise_pred = cfg.noise_pct * std(pred_gt);

    %% Design: 5 built-in + custom periodic kernel
    N_list = [5, 10, 25, 50];
    kernel_list = {'squaredexponential', 'exponential', 'matern32', 'matern52', 'rationalquadratic', @periodicKernel};
    kernel_labels = {'SqExp', 'Matern 1/2', 'Matern 3/2', 'Matern 5/2', 'RatQuad', 'Periodic'};
    n_anticipated_cycles = 3;

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

    hopts = struct('ShowPlots', false, 'Verbose', 0);

    %% Fit GP loop
    for in = 1:length(N_list)
        N = N_list(in);
        if cfg.is_regular
            t_obs = linspace(0, 50, N)';
        else
            t_obs = sort(randperm(51, N) - 1)';
        end
        prey_obs = interp1(t_gt, prey_gt, t_obs, 'linear', 'extrap');
        pred_obs = interp1(t_gt, pred_gt, t_obs, 'linear', 'extrap');
        if cfg.noise_pct > 0
            prey_obs = prey_obs + sigma_noise_prey * randn(size(prey_obs));
            pred_obs = pred_obs + sigma_noise_pred * randn(size(pred_obs));
        end
        t_obs_by_N{in} = t_obs;
        prey_obs_by_N{in} = prey_obs;
        pred_obs_by_N{in} = pred_obs;

        t_range = max(t_obs) - min(t_obs);
        period0 = t_range / n_anticipated_cycles;
        lengthScale0 = period0 / 4;
        theta0_prey = [std(prey_obs), period0, lengthScale0];
        theta0_pred = [std(pred_obs), period0, lengthScale0];

        for ik = 1:length(kernel_list)
            kern = kernel_list{ik};
            is_periodic = isa(kern, 'function_handle');
            try
                % GP for prey
                if is_periodic
                    gpr_prey = fitrgp(t_obs, prey_obs, ...
                        'KernelFunction', @periodicKernel, 'KernelParameters', theta0_prey, ...
                        'BasisFunction', 'constant', 'FitMethod', 'exact', 'PredictMethod', 'exact', 'Standardize', true);
                else
                    base_prey = {'KernelFunction', kern, 'BasisFunction', 'constant', ...
                        'FitMethod', 'exact', 'PredictMethod', 'exact', 'Standardize', true};
                    if cfg.is_auto
                        base_prey = [base_prey, {'OptimizeHyperparameters', 'auto'}, ...
                            {'HyperparameterOptimizationOptions', hopts}];
                    end
                    gpr_prey = fitrgp(t_obs, prey_obs, base_prey{:});
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
                    base_pred = {'KernelFunction', kern, 'BasisFunction', 'constant', ...
                        'FitMethod', 'exact', 'PredictMethod', 'exact', 'Standardize', true};
                    if cfg.is_auto
                        base_pred = [base_pred, {'OptimizeHyperparameters', 'auto'}, ...
                            {'HyperparameterOptimizationOptions', hopts}];
                    end
                    gpr_pred = fitrgp(t_obs, pred_obs, base_pred{:});
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

    %% Plots
    if make_plots
        sam_str = iif(cfg.is_regular, 'Regular', 'Irregular');
        noise_str = sprintf('%.0f%%', cfg.noise_pct * 100);
        tit_suffix = iif(cfg.is_auto, ', Auto', '');
        tit_base = sprintf('GP LV Exp %02d: %s sampling, %s noise%s', cfg.exp_id, sam_str, noise_str, tit_suffix);

        figure;
        subplot(1, 2, 1);
        hold on;
        for ik = 1:length(kernel_labels)
            idx = strcmp(T_point.Kernel, kernel_labels{ik}) & strcmp(T_point.State, 'Prey');
            plot(T_point.N(idx), T_point.RMSE(idx), '-o', 'LineWidth', 1.5, 'DisplayName', kernel_labels{ik});
        end
        xlabel('N'); ylabel('RMSE'); title(sprintf('Exp %02d: Prey — %s, %s%s', cfg.exp_id, sam_str, noise_str, tit_suffix));
        legend('Location', 'best'); grid on; set(gca, 'XTick', N_list);
        subplot(1, 2, 2);
        hold on;
        for ik = 1:length(kernel_labels)
            idx = strcmp(T_point.Kernel, kernel_labels{ik}) & strcmp(T_point.State, 'Predator');
            plot(T_point.N(idx), T_point.RMSE(idx), '-o', 'LineWidth', 1.5, 'DisplayName', kernel_labels{ik});
        end
        xlabel('N'); ylabel('RMSE'); title(sprintf('Exp %02d: Predator — %s, %s%s', cfg.exp_id, sam_str, noise_str, tit_suffix));
        legend('Location', 'best'); grid on; set(gca, 'XTick', N_list);
        sgtitle(tit_base);

        for ik = 1:length(kernel_labels)
            figure('Name', sprintf('LV Exp %02d: %s', cfg.exp_id, kernel_labels{ik}));
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
            sgtitle(sprintf('GP LV Exp %02d: %s — %s, %s%s', cfg.exp_id, kernel_labels{ik}, sam_str, noise_str, tit_suffix));
        end
    end

    %% Save (MAT uses full name; CSVs use exp_id only)
    [mat_fstem, csv_fstem] = get_output_fstems(cfg);
    save(fullfile(out_dir, [mat_fstem '.mat']), 'T_point', 'T_prob', 'T_calib', 't_gt', 'prey_gt', 'pred_gt');
    writetable(T_point, fullfile(out_dir, [csv_fstem '_point_metrics.csv']));
    writetable(T_prob, fullfile(out_dir, [csv_fstem '_prob_metrics.csv']));
    writetable(T_calib, fullfile(out_dir, [csv_fstem '_calib_metrics.csv']));
end

function v = get_opt(opts, name, default)
    if isfield(opts, name), v = opts.(name); else, v = default; end
end

function v = iif(cond, a, b)
    if cond, v = a; else, v = b; end
end

function fstem = get_output_fstem(cfg)
    % Same stem for .mat and .csv (e.g. results_exp_01_regular_0noise)
    sam = iif(cfg.is_regular, 'regular', 'irregular');
    noise = cfg.noise_pct;
    if noise == 0, ns = '0noise'; elseif noise == 0.01, ns = '1noise';
    elseif noise == 0.05, ns = '5noise'; elseif noise == 0.10, ns = '10noise';
    elseif noise == 0.20, ns = '20noise'; else, ns = sprintf('%.0fnoise', noise*100); end
    auto = iif(cfg.is_auto, '_auto', '');
    fstem = sprintf('results_exp_%02d_%s_%s%s', cfg.exp_id, sam, ns, auto);
end
