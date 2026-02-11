function [T_point, T_prob, T_calib] = run_kinetics_experiment(cfg, opts)
%% RUN_KINETICS_EXPERIMENT Run one GP kinetics experiment with given configuration.
%
%   [T_point, T_prob, T_calib] = run_kinetics_experiment(cfg)
%   [T_point, T_prob, T_calib] = run_kinetics_experiment(cfg, opts)
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
%   and one merged .csv with all metrics per experiment. Three states:
%   X, Y, Z (Robertson kinetics). Five built-in kernels only (no periodic).

    if nargin < 2, opts = struct(); end
    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    addpath(fullfile(script_dir, '..'));

    out_dir = get_opt(opts, 'out_dir', script_dir);
    make_plots = get_opt(opts, 'make_plots', true);
    rng_seed = get_opt(opts, 'rng_seed', 42);

    rng(rng_seed);

    %% Ground truth (three states: X, Y, Z)
    [t_gt, x_gt, y_gt, z_gt] = ground_truth_kinetics(500);
    n_rep = get_opt(cfg, 'n_replicates', 1);

    %% Design: N time points, n_replicates per state
    N_list = [5, 10, 25, 50];
    kernel_list = {'squaredexponential', 'exponential', 'matern32', 'matern52', 'rationalquadratic'};
    kernel_labels = {'SqExp', 'Matern 1/2', 'Matern 3/2', 'Matern 5/2', 'RatQuad'};
    state_labels = {'X', 'Y', 'Z'};

    results_point = [];
    results_prob = [];
    results_calib = [];
    pred_ymu_x = cell(length(kernel_list), length(N_list));
    pred_ystd_x = pred_ymu_x;
    pred_ymu_y = pred_ymu_x;
    pred_ystd_y = pred_ymu_x;
    pred_ymu_z = pred_ymu_x;
    pred_ystd_z = pred_ymu_x;
    t_obs_by_N = cell(1, length(N_list));
    x_obs_by_N = t_obs_by_N;
    y_obs_by_N = t_obs_by_N;
    z_obs_by_N = t_obs_by_N;

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
        [t_obs, x_obs] = generate_replicates(t_gt, x_gt, t_points, n_rep, cfg.noise_pct);
        [~, y_obs] = generate_replicates(t_gt, y_gt, t_points, n_rep, cfg.noise_pct);
        [~, z_obs] = generate_replicates(t_gt, z_gt, t_points, n_rep, cfg.noise_pct);
        t_obs_by_N{in} = t_obs;
        x_obs_by_N{in} = x_obs;
        y_obs_by_N{in} = y_obs;
        z_obs_by_N{in} = z_obs;

        for ik = 1:length(kernel_list)
            kern = kernel_list{ik};
            fit_opts = [{'KernelFunction', kern}, base_opts];
            try
                gpr_x = fitrgp(t_obs, x_obs, fit_opts{:});
                [ymu_x, ystd_x, ~] = predict(gpr_x, t_gt);
                pred_ymu_x{ik, in} = ymu_x;
                pred_ystd_x{ik, in} = ystd_x;
                [rmse_x, mae_x, r2_x] = metric_helpers('point_estimate', x_gt, ymu_x);
                [nlpd_x, msll_x, crps_x] = metric_helpers('probabilistic', x_gt, ymu_x, ystd_x, x_obs);
                [smse_x, cov_x, nlml_x] = metric_helpers('calibration', x_gt, ymu_x, ystd_x, x_obs, gpr_x);
                results_point = [results_point; {kernel_labels{ik}, N, n_rep, 'X', rmse_x, mae_x, r2_x}];
                results_prob  = [results_prob;  {kernel_labels{ik}, N, n_rep, 'X', nlpd_x, msll_x, crps_x}];
                results_calib = [results_calib; {kernel_labels{ik}, N, n_rep, 'X', smse_x, cov_x, nlml_x}];

                gpr_y = fitrgp(t_obs, y_obs, fit_opts{:});
                [ymu_y, ystd_y, ~] = predict(gpr_y, t_gt);
                pred_ymu_y{ik, in} = ymu_y;
                pred_ystd_y{ik, in} = ystd_y;
                [rmse_y, mae_y, r2_y] = metric_helpers('point_estimate', y_gt, ymu_y);
                [nlpd_y, msll_y, crps_y] = metric_helpers('probabilistic', y_gt, ymu_y, ystd_y, y_obs);
                [smse_y, cov_y, nlml_y] = metric_helpers('calibration', y_gt, ymu_y, ystd_y, y_obs, gpr_y);
                results_point = [results_point; {kernel_labels{ik}, N, n_rep, 'Y', rmse_y, mae_y, r2_y}];
                results_prob  = [results_prob;  {kernel_labels{ik}, N, n_rep, 'Y', nlpd_y, msll_y, crps_y}];
                results_calib = [results_calib; {kernel_labels{ik}, N, n_rep, 'Y', smse_y, cov_y, nlml_y}];

                gpr_z = fitrgp(t_obs, z_obs, fit_opts{:});
                [ymu_z, ystd_z, ~] = predict(gpr_z, t_gt);
                pred_ymu_z{ik, in} = ymu_z;
                pred_ystd_z{ik, in} = ystd_z;
                [rmse_z, mae_z, r2_z] = metric_helpers('point_estimate', z_gt, ymu_z);
                [nlpd_z, msll_z, crps_z] = metric_helpers('probabilistic', z_gt, ymu_z, ystd_z, z_obs);
                [smse_z, cov_z, nlml_z] = metric_helpers('calibration', z_gt, ymu_z, ystd_z, z_obs, gpr_z);
                results_point = [results_point; {kernel_labels{ik}, N, n_rep, 'Z', rmse_z, mae_z, r2_z}];
                results_prob  = [results_prob;  {kernel_labels{ik}, N, n_rep, 'Z', nlpd_z, msll_z, crps_z}];
                results_calib = [results_calib; {kernel_labels{ik}, N, n_rep, 'Z', smse_z, cov_z, nlml_z}];
            catch me
                warning('Failed: N=%d, kernel=%s. %s', N, kernel_labels{ik}, me.message);
                for s = 1:3
                    st = state_labels{s};
                    results_point = [results_point; {kernel_labels{ik}, N, n_rep, st, NaN, NaN, NaN}];
                    results_prob  = [results_prob;  {kernel_labels{ik}, N, n_rep, st, NaN, NaN, NaN}];
                    results_calib = [results_calib; {kernel_labels{ik}, N, n_rep, st, NaN, NaN, NaN}];
                end
                pred_ymu_x{ik, in} = nan(size(t_gt)); pred_ystd_x{ik, in} = nan(size(t_gt));
                pred_ymu_y{ik, in} = nan(size(t_gt)); pred_ystd_y{ik, in} = nan(size(t_gt));
                pred_ymu_z{ik, in} = nan(size(t_gt)); pred_ystd_z{ik, in} = nan(size(t_gt));
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
        tit_base = sprintf('GP Kinetics Exp %02d: %s sampling, %s noise%s%s', cfg.exp_id, sam_str, noise_str, tit_suffix, rep_str);
        N_total_list = N_list * n_rep;

        figure;
        for is = 1:3
            subplot(1, 3, is);
            hold on;
            for ik = 1:length(kernel_labels)
                idx = strcmp(T_point.Kernel, kernel_labels{ik}) & strcmp(T_point.State, state_labels{is});
                plot(N_total_list, T_point.RMSE(idx), '-o', 'LineWidth', 1.5, 'DisplayName', kernel_labels{ik});
            end
            xlabel('N'); ylabel('RMSE'); title(sprintf('Exp %02d: %s — %s, %s%s%s', cfg.exp_id, state_labels{is}, sam_str, noise_str, tit_suffix, rep_str));
            legend('Location', 'best'); grid on; set(gca, 'XTick', N_total_list);
        end
        sgtitle(tit_base);

        gt_cell = {x_gt, y_gt, z_gt};
        obs_cell = {x_obs_by_N, y_obs_by_N, z_obs_by_N};
        ymu_cell = {pred_ymu_x, pred_ymu_y, pred_ymu_z};
        ystd_cell = {pred_ystd_x, pred_ystd_y, pred_ystd_z};
        colors = {'b', 'g', 'r'};

        for ik = 1:length(kernel_labels)
            figure('Name', sprintf('Kin Exp %02d: %s', cfg.exp_id, kernel_labels{ik}));
            for is = 1:3
                y_gt_s = gt_cell{is};
                for in = 1:length(N_list)
                    subplot(3, 4, (is-1)*4 + in);
                    t_obs = t_obs_by_N{in};
                    y_obs_s = obs_cell{is}{in};
                    ymu = ymu_cell{is}{ik, in};
                    ystd = ystd_cell{is}{ik, in};
                    plot(t_gt, y_gt_s, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Ground truth'); hold on;
                    plot(t_obs, y_obs_s, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
                    plot(t_gt, ymu, '-', 'Color', colors{is}, 'LineWidth', 1.2, 'DisplayName', 'GP mean');
                    if ~all(isnan(ystd))
                        ylo = ymu - 1.96*ystd; yhi = ymu + 1.96*ystd;
                        fill([t_gt; flip(t_gt)], [ylo; flip(yhi)], colors{is}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                    end
                    hold off;
                    xlabel('Time'); ylabel(state_labels{is}); title(sprintf('%s, N = %d', state_labels{is}, N_total_list(in)));
                    legend('Location', 'best', 'FontSize', 8); grid on; xlim([0 50]);
                end
            end
            sgtitle(sprintf('GP Kinetics Exp %02d: %s — %s, %s%s%s', cfg.exp_id, kernel_labels{ik}, sam_str, noise_str, tit_suffix, rep_str));
        end
    end

    %% Save (MAT + single CSV with all metrics per experiment; same stem for both)
    fstem = get_output_fstem(cfg);
    save(fullfile(out_dir, [fstem '.mat']), 'T_point', 'T_prob', 'T_calib', 't_gt', 'x_gt', 'y_gt', 'z_gt');
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
    fstem = sprintf('results_kin_exp_%02d_%s_%s%s%s', cfg.exp_id, sam, ns, auto, rep_suf);
end
