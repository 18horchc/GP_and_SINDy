function [T_point, T_prob, T_calib] = run_logistic_experiment(cfg, opts)
%% RUN_LOGISTIC_EXPERIMENT Run one GP logistic experiment with given configuration.
%
%   [T_point, T_prob, T_calib] = run_logistic_experiment(cfg)
%   [T_point, T_prob, T_calib] = run_logistic_experiment(cfg, opts)
%
%   cfg: struct with fields:
%        .exp_id     - experiment number (1-20), used for output filenames
%        .noise_pct  - noise fraction (0, 0.01, 0.05, 0.10, 0.20)
%        .is_regular - true = linspace sampling, false = randperm sampling
%        .is_auto    - true = OptimizeHyperparameters='auto', false = default
%
%   opts (optional): struct with fields:
%        .out_dir    - where to save (default: folder containing this file)
%        .make_plots - true/false (default true)
%        .rng_seed   - for reproducibility when noise or irregular (default 42)
%
%   Output: T_point, T_prob, T_calib (tables). Also saves .mat and .csv files
%   with names matching aggregate_logistic_metrics expectations.

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
    sigma_noise = cfg.noise_pct * std(y_gt);

    %% Design
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
            t_obs = linspace(0, 50, N)';
        else
            t_obs = sort(randperm(51, N) - 1)';
        end
        y_obs = interp1(t_gt, y_gt, t_obs, 'linear', 'extrap');
        if cfg.noise_pct > 0
            y_obs = y_obs + sigma_noise * randn(size(y_obs));
        end
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

    %% Save (MAT uses full name; CSVs use exp_id only, matching original convention)
    [mat_fstem, csv_fstem] = get_output_fstems(cfg);
    save(fullfile(out_dir, [mat_fstem '.mat']), 'T_point', 'T_prob', 'T_calib', 't_gt', 'y_gt');
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

function [mat_fstem, csv_fstem] = get_output_fstems(cfg)
    % MAT: full name for aggregate_logistic_metrics.m; CSV: exp_id only
    sam = iif(cfg.is_regular, 'regular', 'irregular');
    noise = cfg.noise_pct;
    if noise == 0, ns = '0noise'; elseif noise == 0.01, ns = '1noise';
    elseif noise == 0.05, ns = '5noise'; elseif noise == 0.10, ns = '10noise';
    elseif noise == 0.20, ns = '20noise'; else, ns = sprintf('%.0fnoise', noise*100); end
    auto = iif(cfg.is_auto, '_auto', '');
    mat_fstem = sprintf('results_exp_%02d_%s_%s%s', cfg.exp_id, sam, ns, auto);
    csv_fstem = sprintf('results_exp_%02d', cfg.exp_id);
end
