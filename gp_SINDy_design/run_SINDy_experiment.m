function [T_point, T_prob, T_calib] = run_SINDy_experiment(cfg, opts)
%% RUN_SINDY_EXPERIMENT Run one GP SINDy experiment with given configuration.
%
%   [T_point, T_prob, T_calib] = run_SINDy_experiment(cfg)
%   [T_point, T_prob, T_calib] = run_SINDy_experiment(cfg, opts)
%
%   cfg: struct with fields:
%        .exp_id      - experiment number (1-60), used for output filenames
%        .noise_pct   - noise fraction (0, 0.01, 0.05, 0.10, 0.20)
%        .is_regular  - true = linspace sampling, false = randperm sampling
%        .is_auto     - true = OptimizeHyperparameters='auto', false = default
%        .n_replicates - 1, 3, or 8; observations per time point per state (default 1)
%
%   opts (optional): struct with fields:
%        .out_dir    - where to save (default: folder containing this file)
%        .make_plots - true/false (default true)
%        .rng_seed   - for reproducibility (default 42)
%
%   Output: T_point, T_prob, T_calib (tables with Kernel, N, State). Saves .mat
%   and one merged .csv with all metrics per experiment. Two states:
%   M1, M2 (SINDy+MCMC microglia). Five built-in kernels only.

    if nargin < 2, opts = struct(); end
    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    addpath(fullfile(script_dir, '..'));

    out_dir = get_opt(opts, 'out_dir', script_dir);
    make_plots = get_opt(opts, 'make_plots', true);
    rng_seed = get_opt(opts, 'rng_seed', 42);

    rng(rng_seed);

    %% Ground truth (two states: M1, M2 - SINDy+MCMC microglia)
    [t_gt, M1_gt, M2_gt] = ground_truth_SINDy(500);
    n_rep = get_opt(cfg, 'n_replicates', 1);

    %% Design: N time points, n_replicates per state; 5 built-in kernels
    N_list = [5, 10, 25, 50];
    kernel_list = {'squaredexponential', 'exponential', 'matern32', 'matern52', 'rationalquadratic'};
    kernel_labels = {'SqExp', 'Matern 1/2', 'Matern 3/2', 'Matern 5/2', 'RatQuad'};
    state_labels = {'M1', 'M2'};

    results_point = [];
    results_prob = [];
    results_calib = [];
    pred_ymu_M1  = cell(length(kernel_list), length(N_list));
    pred_ystd_M1 = pred_ymu_M1;
    pred_ymu_M2  = pred_ymu_M1;
    pred_ystd_M2 = pred_ymu_M1;
    t_obs_by_N  = cell(1, length(N_list));
    M1_obs_by_N = t_obs_by_N;
    M2_obs_by_N = t_obs_by_N;

    base_opts = {'BasisFunction', 'constant', 'FitMethod', 'exact', ...
        'PredictMethod', 'exact', 'Standardize', true};
    if cfg.is_auto
        base_opts = [base_opts, {'OptimizeHyperparameters', 'auto'}, ...
            {'HyperparameterOptimizationOptions', struct('ShowPlots', false, 'Verbose', 0)}];
    end

    %% Fit GP loop
    for in = 1:length(N_list)
        N = N_list(in);
        if cfg.is_regular
            t_points = linspace(0, 50, N)';
        else
            t_points = sort(randperm(51, N) - 1)';
        end
        [t_obs, M1_obs] = generate_replicates(t_gt, M1_gt, t_points, n_rep, cfg.noise_pct);
        [~, M2_obs] = generate_replicates(t_gt, M2_gt, t_points, n_rep, cfg.noise_pct);
        t_obs_by_N{in} = t_obs;
        M1_obs_by_N{in} = M1_obs;
        M2_obs_by_N{in} = M2_obs;

        for ik = 1:length(kernel_list)
            kern = kernel_list{ik};
            fit_opts = [{'KernelFunction', kern}, base_opts];
            try
                gpr_M1 = fitrgp(t_obs, M1_obs, fit_opts{:});
                [ymu_M1, ystd_M1, ~] = predict(gpr_M1, t_gt);
                pred_ymu_M1{ik, in} = ymu_M1;
                pred_ystd_M1{ik, in} = ystd_M1;
                [rmse_1, mae_1, r2_1] = metric_helpers('point_estimate', M1_gt, ymu_M1);
                [nlpd_1, msll_1, crps_1] = metric_helpers('probabilistic', M1_gt, ymu_M1, ystd_M1, M1_obs);
                [smse_1, cov_1, nlml_1] = metric_helpers('calibration', M1_gt, ymu_M1, ystd_M1, M1_obs, gpr_M1);
                results_point = [results_point; {kernel_labels{ik}, N, n_rep, 'M1', rmse_1, mae_1, r2_1}];
                results_prob  = [results_prob;  {kernel_labels{ik}, N, n_rep, 'M1', nlpd_1, msll_1, crps_1}];
                results_calib = [results_calib; {kernel_labels{ik}, N, n_rep, 'M1', smse_1, cov_1, nlml_1}];

                gpr_M2 = fitrgp(t_obs, M2_obs, fit_opts{:});
                [ymu_M2, ystd_M2, ~] = predict(gpr_M2, t_gt);
                pred_ymu_M2{ik, in} = ymu_M2;
                pred_ystd_M2{ik, in} = ystd_M2;
                [rmse_2, mae_2, r2_2] = metric_helpers('point_estimate', M2_gt, ymu_M2);
                [nlpd_2, msll_2, crps_2] = metric_helpers('probabilistic', M2_gt, ymu_M2, ystd_M2, M2_obs);
                [smse_2, cov_2, nlml_2] = metric_helpers('calibration', M2_gt, ymu_M2, ystd_M2, M2_obs, gpr_M2);
                results_point = [results_point; {kernel_labels{ik}, N, n_rep, 'M2', rmse_2, mae_2, r2_2}];
                results_prob  = [results_prob;  {kernel_labels{ik}, N, n_rep, 'M2', nlpd_2, msll_2, crps_2}];
                results_calib = [results_calib; {kernel_labels{ik}, N, n_rep, 'M2', smse_2, cov_2, nlml_2}];
            catch me
                warning('Failed: N=%d, kernel=%s. %s', N, kernel_labels{ik}, me.message);
                results_point = [results_point; {kernel_labels{ik}, N, n_rep, 'M1', NaN, NaN, NaN}; {kernel_labels{ik}, N, n_rep, 'M2', NaN, NaN, NaN}];
                results_prob  = [results_prob;  {kernel_labels{ik}, N, n_rep, 'M1', NaN, NaN, NaN}; {kernel_labels{ik}, N, n_rep, 'M2', NaN, NaN, NaN}];
                results_calib = [results_calib; {kernel_labels{ik}, N, n_rep, 'M1', NaN, NaN, NaN}; {kernel_labels{ik}, N, n_rep, 'M2', NaN, NaN, NaN}];
                pred_ymu_M1{ik, in} = nan(size(t_gt)); pred_ystd_M1{ik, in} = nan(size(t_gt));
                pred_ymu_M2{ik, in} = nan(size(t_gt)); pred_ystd_M2{ik, in} = nan(size(t_gt));
            end
        end
    end

    T_point = cell2table(results_point, 'VariableNames', {'Kernel', 'N', 'N_replicates', 'State', 'RMSE', 'MAE', 'R2'});
    T_prob  = cell2table(results_prob,  'VariableNames', {'Kernel', 'N', 'N_replicates', 'State', 'NLPD', 'MSLL', 'CRPS'});
    T_calib = cell2table(results_calib, 'VariableNames', {'Kernel', 'N', 'N_replicates', 'State', 'sMSE', 'Coverage', 'NLML'});

    %% Plots
    if make_plots
        sam_str = iif(cfg.is_regular, 'Regular', 'Irregular');
        noise_str = sprintf('%.0f%%', cfg.noise_pct * 100);
        tit_suffix = iif(cfg.is_auto, ', Auto', '');
        rep_str = iif(n_rep > 1, sprintf(' — Replicates, %d', n_rep), '');
        tit_base = sprintf('GP SINDy Exp %02d: %s sampling, %s noise%s%s', cfg.exp_id, sam_str, noise_str, tit_suffix, rep_str);
        N_total_list = N_list * n_rep;

        figure;
        subplot(1, 2, 1);
        hold on;
        for ik = 1:length(kernel_labels)
            idx = strcmp(T_point.Kernel, kernel_labels{ik}) & strcmp(T_point.State, 'M1');
            plot(N_total_list, T_point.RMSE(idx), '-o', 'LineWidth', 1.5, 'DisplayName', kernel_labels{ik});
        end
        xlabel('N'); ylabel('RMSE'); title(sprintf('Exp %02d: M1 — %s, %s%s%s', cfg.exp_id, sam_str, noise_str, tit_suffix, rep_str));
        legend('Location', 'best'); grid on; set(gca, 'XTick', N_total_list);
        subplot(1, 2, 2);
        hold on;
        for ik = 1:length(kernel_labels)
            idx = strcmp(T_point.Kernel, kernel_labels{ik}) & strcmp(T_point.State, 'M2');
            plot(N_total_list, T_point.RMSE(idx), '-o', 'LineWidth', 1.5, 'DisplayName', kernel_labels{ik});
        end
        xlabel('N'); ylabel('RMSE'); title(sprintf('Exp %02d: M2 — %s, %s%s%s', cfg.exp_id, sam_str, noise_str, tit_suffix, rep_str));
        legend('Location', 'best'); grid on; set(gca, 'XTick', N_total_list);
        sgtitle(tit_base);

        for ik = 1:length(kernel_labels)
            figure('Name', sprintf('SINDy Exp %02d: %s', cfg.exp_id, kernel_labels{ik}));
            for in = 1:length(N_list)
                subplot(2, 4, in);
                t_obs = t_obs_by_N{in};
                M1_obs = M1_obs_by_N{in};
                ymu = pred_ymu_M1{ik, in};
                ystd = pred_ystd_M1{ik, in};
                plot(t_gt, M1_gt, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Ground truth'); hold on;
                plot(t_obs, M1_obs, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
                plot(t_gt, ymu, 'b-', 'LineWidth', 1.2, 'DisplayName', 'GP mean');
                if ~all(isnan(ystd))
                    ylo = ymu - 1.96*ystd; yhi = ymu + 1.96*ystd;
                    fill([t_gt; flip(t_gt)], [ylo; flip(yhi)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                end
                hold off; xlabel('Time'); ylabel('M1 (cells/mm^2)'); title(sprintf('M1, N = %d', N_total_list(in)));
                legend('Location', 'best', 'FontSize', 8); grid on; xlim([0 50]);

                subplot(2, 4, 4 + in);
                M2_obs = M2_obs_by_N{in};
                ymu = pred_ymu_M2{ik, in};
                ystd = pred_ystd_M2{ik, in};
                plot(t_gt, M2_gt, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Ground truth'); hold on;
                plot(t_obs, M2_obs, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
                plot(t_gt, ymu, 'r-', 'LineWidth', 1.2, 'DisplayName', 'GP mean');
                if ~all(isnan(ystd))
                    ylo = ymu - 1.96*ystd; yhi = ymu + 1.96*ystd;
                    fill([t_gt; flip(t_gt)], [ylo; flip(yhi)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                end
                hold off; xlabel('Time'); ylabel('M2 (cells/mm^2)'); title(sprintf('M2, N = %d', N_total_list(in)));
                legend('Location', 'best', 'FontSize', 8); grid on; xlim([0 50]);
            end
            sgtitle(sprintf('GP SINDy Exp %02d: %s — %s, %s%s%s', cfg.exp_id, kernel_labels{ik}, sam_str, noise_str, tit_suffix, rep_str));
        end
    end

    %% Save (MAT + single CSV with all metrics per experiment; same stem for both)
    fstem = get_output_fstem(cfg);
    save(fullfile(out_dir, [fstem '.mat']), 'T_point', 'T_prob', 'T_calib', 't_gt', 'M1_gt', 'M2_gt');
    T_merged = join(T_point, T_prob, 'Keys', {'Kernel', 'N', 'N_replicates', 'State'});
    T_merged = join(T_merged, T_calib, 'Keys', {'Kernel', 'N', 'N_replicates', 'State'});
    writetable(T_merged, fullfile(out_dir, [fstem '.csv']));
end

function v = get_opt(opts, name, default)
    if isfield(opts, name), v = opts.(name); else, v = default; end
end

function v = iif(cond, a, b)
    if cond, v = a; else, v = b; end
end

function [t_obs, y_obs] = generate_replicates(t_gt, y_gt, t_points, n_rep, noise_pct)
    y_true_at_t = interp1(t_gt, y_gt, t_points, 'linear', 'extrap');
    t_obs = repelem(t_points, n_rep);
    y_obs = repelem(y_true_at_t, n_rep);
    if noise_pct > 0
        sigma_floor = 1e-10 * std(y_gt);
        for i = 1:length(t_points)
            y_true_i = y_true_at_t(i);
            sigma_i = max(noise_pct * abs(y_true_i), sigma_floor);
            idx = (i-1)*n_rep + (1:n_rep);
            y_obs(idx) = y_true_i + sigma_i * randn(n_rep, 1);
        end
    end
    t_obs = t_obs(:);
    y_obs = y_obs(:);
end

function fstem = get_output_fstem(cfg)
    sam = iif(cfg.is_regular, 'regular', 'irregular');
    noise = cfg.noise_pct;
    if noise == 0, ns = '0noise'; elseif noise == 0.01, ns = '1noise';
    elseif noise == 0.05, ns = '5noise'; elseif noise == 0.10, ns = '10noise';
    elseif noise == 0.20, ns = '20noise'; else, ns = sprintf('%.0fnoise', noise*100); end
    auto = iif(cfg.is_auto, '_auto', '');
    n_rep = get_opt(cfg, 'n_replicates', 1);
    rep_suf = iif(n_rep > 1, sprintf('_%drep', n_rep), '');
    fstem = sprintf('results_SINDy_exp_%02d_%s_%s%s%s', cfg.exp_id, sam, ns, auto, rep_suf);
end
