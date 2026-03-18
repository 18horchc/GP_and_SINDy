function [T_compare, gpr_naive, gpr_baseline] = run_naive_vs_baseline_mean(opts)
%% RUN_NAIVE_VS_BASELINE_MEAN  Compare naive GP vs baseline-mean encoding on acute transient.
%
%   [T_compare, gpr_naive, gpr_baseline] = run_naive_vs_baseline_mean()
%   [T_compare, gpr_naive, gpr_baseline] = run_naive_vs_baseline_mean(opts)
%
%   Data: t = 0,1,2,3,5,7,14, 2 rep, 5% noise. Domain [0,15].
%   Naive: SE kernel, constant mean (fitted).
%   Baseline-encoded: Same + pseudo-observations at (0, baseline) and (15, baseline)
%   to encode prior belief that the curve returns to baseline at boundaries.
%
%   opts: .baseline (5), .rng_seed (42), .make_plots (true), .out_dir (results folder)
%
%   Returns: T_compare (metrics table with Model column), both fitted GPs.
%
%   See also: run_naive_acute_transient, ground_truth_bio, generate_bio_data

    if nargin < 1, opts = struct(); end

    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    addpath(fullfile(script_dir, '..'));

    baseline  = get_opt(opts, 'baseline', 5);
    rng_seed  = get_opt(opts, 'rng_seed', 42);
    make_plots = get_opt(opts, 'make_plots', true);
    out_dir   = get_opt(opts, 'out_dir', fullfile(script_dir, 'results'));

    rng(rng_seed);

    %% Data: t = 0,1,2,3,5,7,14, 2 rep, 5% noise, domain [0,15]
    data_cfg = struct('t_points', [0; 1; 2; 3; 5; 7; 14], 'n_replicates', 2, ...
        'noise_pct', 0.05, 't_domain', [0, 15]);

    params = struct('t_domain', [0, 15]);
    [t_gt, y_gt, ~] = ground_truth_bio('acute_transient', 500, params);
    [t_obs, y_obs] = generate_bio_data(t_gt, y_gt, data_cfg, rng_seed);

    %% Naive GP
    rng(rng_seed);
    gpr_naive = fitrgp(t_obs, y_obs, 'KernelFunction', 'squaredexponential', ...
        'BasisFunction', 'constant', 'FitMethod', 'exact', 'PredictMethod', 'exact', ...
        'Standardize', true, 'OptimizeHyperparameters', 'auto', ...
        'HyperparameterOptimizationOptions', struct('ShowPlots', false, 'Verbose', 0, 'UseParallel', false));

    [ymu_naive, ystd_naive, ~] = predict(gpr_naive, t_gt);

    %% Baseline-mean encoding: add pseudo-obs at boundaries (0, baseline) and (15, baseline)
    t_max = 15;
    t_fit = [t_obs; 0; t_max];
    y_fit = [y_obs; baseline; baseline];

    rng(rng_seed);
    gpr_baseline = fitrgp(t_fit, y_fit, 'KernelFunction', 'squaredexponential', ...
        'BasisFunction', 'constant', 'FitMethod', 'exact', 'PredictMethod', 'exact', ...
        'Standardize', true, 'OptimizeHyperparameters', 'auto', ...
        'HyperparameterOptimizationOptions', struct('ShowPlots', false, 'Verbose', 0, 'UseParallel', false));

    [ymu_baseline, ystd_baseline, ~] = predict(gpr_baseline, t_gt);

    %% Metrics for both
    [rmse_n, mae_n, r2_n] = metric_helpers('point_estimate', y_gt, ymu_naive);
    [nlpd_n, msll_n, crps_n] = metric_helpers('probabilistic', y_gt, ymu_naive, ystd_naive, y_obs);
    [smse_n, cov_n, nlml_n] = metric_helpers('calibration', y_gt, ymu_naive, ystd_naive, y_obs, gpr_naive);

    [rmse_b, mae_b, r2_b] = metric_helpers('point_estimate', y_gt, ymu_baseline);
    [nlpd_b, msll_b, crps_b] = metric_helpers('probabilistic', y_gt, ymu_baseline, ystd_baseline, y_obs);
    [smse_b, cov_b, nlml_b] = metric_helpers('calibration', y_gt, ymu_baseline, ystd_baseline, y_obs, gpr_baseline);

    T_compare = table(...
        {'naive'; 'baseline_mean'}, ...
        [rmse_n; rmse_b], [mae_n; mae_b], [r2_n; r2_b], ...
        [nlpd_n; nlpd_b], [msll_n; msll_b], [crps_n; crps_b], ...
        [smse_n; smse_b], [cov_n; cov_b], [nlml_n; nlml_b], ...
        'VariableNames', {'Model','RMSE','MAE','R2','NLPD','MSLL','CRPS','sMSE','Coverage','NLML'});

    fprintf('\nNaive vs Baseline-Mean (baseline=%.1f):\n', baseline);
    disp(T_compare);

    %% Save
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    fstem = sprintf('results_naive_vs_baseline_mean_b%d', round(baseline));
    save(fullfile(out_dir, [fstem '.mat']), 'T_compare', 't_gt', 'y_gt', 't_obs', 'y_obs', ...
        'ymu_naive', 'ystd_naive', 'ymu_baseline', 'ystd_baseline', ...
        'gpr_naive', 'gpr_baseline', 'baseline', 'data_cfg');
    writetable(T_compare, fullfile(out_dir, [fstem '.csv']));

    %% Plot: both models with uncertainty bounds
    if make_plots
        figure('Name', 'Naive vs Baseline-Mean Encoding');
        plot(t_gt, y_gt, 'k-', 'LineWidth', 2, 'DisplayName', 'Ground truth');
        hold on;
        plot(t_obs, y_obs, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');

        % Naive
        plot(t_gt, ymu_naive, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Naive');
        ylo_n = ymu_naive - 1.96 * ystd_naive;
        yhi_n = ymu_naive + 1.96 * ystd_naive;
        fill([t_gt; flip(t_gt)], [ylo_n; flip(yhi_n)], 'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

        % Baseline-encoded
        plot(t_gt, ymu_baseline, 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('Baseline mean (b=%.0f)', baseline));
        ylo_b = ymu_baseline - 1.96 * ystd_baseline;
        yhi_b = ymu_baseline + 1.96 * ystd_baseline;
        fill([t_gt; flip(t_gt)], [ylo_b; flip(yhi_b)], 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

        hold off;
        xlabel('t');
        ylabel('f(t)');
        title(sprintf('Acute Transient: Naive vs Baseline-Mean (b=%.0f), t=0,1,2,3,5,7,14, 2 rep, 5%% noise', baseline));
        legend('Location', 'best');
        grid on;
        xlim([0 15]);
    end
end

function v = get_opt(opts, name, default)
    if isempty(opts) || ~isfield(opts, name)
        v = default;
    else
        v = opts.(name);
    end
end
