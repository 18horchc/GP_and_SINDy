function [T_point, T_prob, T_calib] = run_logistic_experiment(cfg, opts)
%% RUN_LOGISTIC_EXPERIMENT Run one GP logistic experiment with given configuration.
%
%   [T_point, T_prob, T_calib] = run_logistic_experiment(cfg)
%   [T_point, T_prob, T_calib] = run_logistic_experiment(cfg, opts)
%
%   cfg: struct with fields:
%        .exp_id      - experiment number (1-60), used for output filenames
%        .noise_pct   - noise fraction (0, 0.01, 0.05, 0.10, 0.20)
%        .is_regular  - true = linspace sampling, false = randperm sampling
%        .is_auto     - true = OptimizeHyperparameters='auto', false = default
%        .n_replicates - 1, 3, or 8; observations per time point (default 1)
%
%   opts (optional): struct with fields:
%        .out_dir    - where to save (default: folder containing this file)
%        .make_plots - true/false (default true)
%        .rng_seed   - for reproducibility (default 42)
%
%   Output: T_point, T_prob, T_calib (tables). Also saves .mat and one merged
%   .csv with all metrics per experiment.
%
%   N_list: [5, 10, 25, 50] time points. With n_replicates>1, each time
%   point has n_replicates observations (slightly different via proportional
%   replicate noise: sigma = noise_pct * y_true(t), realistic for growth data).

    if nargin < 2, opts = struct(); end
    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    addpath(fullfile(script_dir, '..'));

    out_dir = get_opt(opts, 'out_dir', script_dir);
    make_plots = get_opt(opts, 'make_plots', true);
    rng_seed = get_opt(opts, 'rng_seed', 42);

    rng(rng_seed);

    %% Ground truth
    [t_gt, y_gt] = ground_truth_logistic(500);
    t_max = 50;
    n_rep = get_opt(cfg, 'n_replicates', 1);

    %% Design: N time points, n_replicates observations at each
    N_list = [5, 10, 25, 50];
    kernel_list = {'squaredexponential', 'exponential', 'matern32', 'matern52', 'rationalquadratic'};
    kernel_labels = {'SqExp', 'Matern 1/2', 'Matern 3/2', 'Matern 5/2', 'RatQuad'};

    results_point = [];
    results_prob = [];
    results_calib = [];
    pred_ymu = cell(length(kernel_list), length(N_list));
    pred_ystd = cell(length(kernel_list), length(N_list));
    t_obs_by_N = cell(1, length(N_list));
    y_obs_by_N = cell(1, length(N_list));

    %% Fit GP loop
    for in = 1:length(N_list)
        N = N_list(in);
        if cfg.is_regular
            t_points = linspace(0, t_max, N)';
        else
            t_points = sort(randperm(51, N) - 1)';
        end
        % Generate n_replicates at each time point (systematic, realistic)
        [t_obs, y_obs] = generate_replicates(t_gt, y_gt, t_points, n_rep, cfg.noise_pct);
        t_obs_by_N{in} = t_obs;
        y_obs_by_N{in} = y_obs;

        base_opts = {'BasisFunction', 'constant', 'FitMethod', 'exact', ...
            'PredictMethod', 'exact', 'Standardize', true};
        if cfg.is_auto
            base_opts = [base_opts, {'OptimizeHyperparameters', 'auto'}, ...
                {'HyperparameterOptimizationOptions', struct('ShowPlots', false, 'Verbose', 0)}];
        end

        for ik = 1:length(kernel_list)
            kern = kernel_list{ik};
            fit_opts = [{'KernelFunction', kern}, base_opts];
            try
                gpr = fitrgp(t_obs, y_obs, fit_opts{:});
                [ymu, ystd, ~] = predict(gpr, t_gt);
                pred_ymu{ik, in} = ymu;
                pred_ystd{ik, in} = ystd;

                [rmse, mae, r2] = metric_helpers('point_estimate', y_gt, ymu);
                results_point = [results_point; {kernel_labels{ik}, N, n_rep, rmse, mae, r2}];

                [nlpd, msll, crps] = metric_helpers('probabilistic', y_gt, ymu, ystd, y_obs);
                results_prob = [results_prob; {kernel_labels{ik}, N, n_rep, nlpd, msll, crps}];

                [smse, coverage, nlml] = metric_helpers('calibration', y_gt, ymu, ystd, y_obs, gpr);
                results_calib = [results_calib; {kernel_labels{ik}, N, n_rep, smse, coverage, nlml}];
            catch me
                warning('Failed: N=%d, kernel=%s. %s', N, kern, me.message);
                results_point = [results_point; {kernel_labels{ik}, N, n_rep, NaN, NaN, NaN}];
                results_prob = [results_prob; {kernel_labels{ik}, N, n_rep, NaN, NaN, NaN}];
                results_calib = [results_calib; {kernel_labels{ik}, N, n_rep, NaN, NaN, NaN}];
                pred_ymu{ik, in} = nan(size(t_gt));
                pred_ystd{ik, in} = nan(size(t_gt));
            end
        end
    end

    T_point = cell2table(results_point, 'VariableNames', {'Kernel', 'N', 'N_replicates', 'RMSE', 'MAE', 'R2'});
    T_prob = cell2table(results_prob, 'VariableNames', {'Kernel', 'N', 'N_replicates', 'NLPD', 'MSLL', 'CRPS'});
    T_calib = cell2table(results_calib, 'VariableNames', {'Kernel', 'N', 'N_replicates', 'sMSE', 'Coverage', 'NLML'});

    %% Plots
    if make_plots
        sam_str = iif(cfg.is_regular, 'Regular', 'Irregular');
        noise_str = sprintf('%.0f%%', cfg.noise_pct * 100);
        auto_str = iif(cfg.is_auto, 'Auto', 'Default');
        tit_base = sprintf('GP Logistic Exp %02d: %s sampling, %s noise%s', ...
            cfg.exp_id, sam_str, noise_str, iif(cfg.is_auto, ', Auto', ''));

        figure;
        hold on;
        for ik = 1:length(kernel_labels)
            idx = strcmp(T_point.Kernel, kernel_labels{ik});
            plot(T_point.N(idx), T_point.RMSE(idx), '-o', 'LineWidth', 1.5, 'DisplayName', kernel_labels{ik});
        end
        xlabel('Number of points (N)'); ylabel('RMSE'); title(tit_base);
        legend('Location', 'best'); grid on; set(gca, 'XTick', N_list);

        for ik = 1:length(kernel_labels)
            figure('Name', sprintf('Exp %02d: %s', cfg.exp_id, kernel_labels{ik}));
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
            sgtitle(sprintf('Exp %02d: %s — %s, %s%s', cfg.exp_id, kernel_labels{ik}, sam_str, noise_str, iif(cfg.is_auto, ', Auto', '')));
        end
    end

    %% Save (MAT + single CSV with all metrics per experiment; same stem for both)
    fstem = get_output_fstem(cfg);
    save(fullfile(out_dir, [fstem '.mat']), 'T_point', 'T_prob', 'T_calib', 't_gt', 'y_gt');
    T_merged = join(T_point, T_prob, 'Keys', {'Kernel', 'N', 'N_replicates'});
    T_merged = join(T_merged, T_calib, 'Keys', {'Kernel', 'N', 'N_replicates'});
    writetable(T_merged, fullfile(out_dir, [fstem '.csv']));
end

function v = get_opt(opts, name, default)
    if isfield(opts, name), v = opts.(name); else, v = default; end
end

function v = iif(cond, a, b)
    if cond, v = a; else, v = b; end
end

function [t_obs, y_obs] = generate_replicates(t_gt, y_gt, t_points, n_rep, noise_pct)
% Generate n_rep observations at each time point. Replicates differ via
% proportional noise (sigma = noise_pct * y_true) — realistic for growth
% data where variance often scales with the mean. Systematic (reproducible
% with fixed rng). When noise_pct=0, replicates are identical.
    y_true_at_t = interp1(t_gt, y_gt, t_points, 'linear', 'extrap');
    t_obs = repelem(t_points, n_rep);
    y_obs = repelem(y_true_at_t, n_rep);
    if noise_pct > 0 && n_rep > 1
        % Proportional replicate noise: sigma_i = noise_pct * y_true(t_i).
        % Realistic for growth data (variance scales with mean). Floor for stability.
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
    % Same stem for .mat and .csv (e.g. results_exp_01_regular_0noise or ..._3rep)
    sam = iif(cfg.is_regular, 'regular', 'irregular');
    noise = cfg.noise_pct;
    if noise == 0, ns = '0noise'; elseif noise == 0.01, ns = '1noise';
    elseif noise == 0.05, ns = '5noise'; elseif noise == 0.10, ns = '10noise';
    elseif noise == 0.20, ns = '20noise'; else, ns = sprintf('%.0fnoise', noise*100); end
    auto = iif(cfg.is_auto, '_auto', '');
    n_rep = get_opt(cfg, 'n_replicates', 1);
    rep_suf = iif(n_rep > 1, sprintf('_%drep', n_rep), '');
    fstem = sprintf('results_exp_%02d_%s_%s%s%s', cfg.exp_id, sam, ns, auto, rep_suf);
end
