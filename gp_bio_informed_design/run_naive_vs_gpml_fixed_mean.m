function [T_compare, hyp_naive, hyp_fixed] = run_naive_vs_gpml_fixed_mean(opts)
%% RUN_NAIVE_VS_GPML_FIXED_MEAN  Compare naive GP vs GPML fixed prior mean (both use GPML).
%
%   [T_compare, hyp_naive, hyp_fixed] = run_naive_vs_gpml_fixed_mean()
%   [T_compare, hyp_naive, hyp_fixed] = run_naive_vs_gpml_fixed_mean(opts)
%
%   Data: t = 0,1,2,3,5,7,14, 2 rep, 5% noise. Domain [0,15].
%   Naive: GPML meanConst + covSEiso, all hyperparameters optimized.
%   Fixed mean: GPML meanConst(5) fixed, cov and lik optimized.
%
%   opts: .baseline (5), .rng_seed (42), .make_plots (true), .out_dir, .gpml_path
%
%   Requires: GPML toolbox. Pass opts.gpml_path if not on path.
%
%   See also: fit_gp_gpml_naive, fit_gp_gpml_fixed_mean, run_naive_acute_transient

    if nargin < 1, opts = struct(); end

    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    addpath(fullfile(script_dir, '..'));

    baseline  = get_opt(opts, 'baseline', 5);
    rng_seed  = get_opt(opts, 'rng_seed', 42);
    make_plots = get_opt(opts, 'make_plots', true);
    out_dir   = get_opt(opts, 'out_dir', fullfile(script_dir, 'results'));
    gpml_path = get_opt(opts, 'gpml_path', []);

    rng(rng_seed);

    %% Data
    data_cfg = struct('t_points', [0; 1; 2; 3; 5; 7; 14], 'n_replicates', 2, ...
        'noise_pct', 0.05, 't_domain', [0, 15]);

    params = struct('t_domain', [0, 15]);
    [t_gt, y_gt, ~] = ground_truth_bio('acute_transient', 500, params);
    [t_obs, y_obs] = generate_bio_data(t_gt, y_gt, data_cfg, rng_seed);

    %% Both models use GPML
    try
        rng(rng_seed);
        [ymu_naive, ys2_naive, hyp_naive] = fit_gp_gpml_naive(t_obs, y_obs, t_gt, gpml_path);
        ystd_naive = sqrt(ys2_naive);

        [ymu_gpml, ys2_gpml, hyp_fixed] = fit_gp_gpml_fixed_mean(t_obs, y_obs, t_gt, baseline, gpml_path);
        ystd_gpml = sqrt(ys2_gpml);
        gpml_ok = true;
    catch me
        warning('GPML failed: %s. Install from https://gitlab.com/hnickisch/gpml-matlab', me.message);
        ymu_naive = nan(size(t_gt)); ystd_naive = nan(size(t_gt)); hyp_naive = [];
        ymu_gpml = nan(size(t_gt)); ystd_gpml = nan(size(t_gt)); hyp_fixed = [];
        gpml_ok = false;
    end

    %% Metrics
    [rmse_n, mae_n, r2_n] = metric_helpers('point_estimate', y_gt, ymu_naive);
    [nlpd_n, msll_n, crps_n] = metric_helpers('probabilistic', y_gt, ymu_naive, ystd_naive, y_obs);
    [smse_n, cov_n, ~] = metric_helpers('calibration', y_gt, ymu_naive, ystd_naive, y_obs, []);
    nlml_n = gpml_nlz_if_ok(gpml_ok, hyp_naive, t_obs, y_obs);

    if gpml_ok
        [rmse_g, mae_g, r2_g] = metric_helpers('point_estimate', y_gt, ymu_gpml);
        [nlpd_g, msll_g, crps_g] = metric_helpers('probabilistic', y_gt, ymu_gpml, ystd_gpml, y_obs);
        [smse_g, cov_g, ~] = metric_helpers('calibration', y_gt, ymu_gpml, ystd_gpml, y_obs, []);
        nlml_g = gpml_nlz_if_ok(true, hyp_fixed, t_obs, y_obs);
    else
        [rmse_g, mae_g, r2_g, nlpd_g, msll_g, crps_g, smse_g, cov_g, nlml_g] = deal(NaN);
    end

    T_compare = table(...
        {'naive'; 'gpml_fixed_mean'}, ...
        [rmse_n; rmse_g], [mae_n; mae_g], [r2_n; r2_g], ...
        [nlpd_n; nlpd_g], [msll_n; msll_g], [crps_n; crps_g], ...
        [smse_n; smse_g], [cov_n; cov_g], [nlml_n; nlml_g], ...
        'VariableNames', {'Model','RMSE','MAE','R2','NLPD','MSLL','CRPS','sMSE','Coverage','NLML'});

    fprintf('\nNaive vs GPML Fixed Mean (baseline=%.1f):\n', baseline);
    disp(T_compare);

    %% Save
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    fstem = sprintf('results_naive_vs_gpml_fixed_mean_b%d', round(baseline));
    save(fullfile(out_dir, [fstem '.mat']), 'T_compare', 't_gt', 'y_gt', 't_obs', 'y_obs', ...
        'ymu_naive', 'ystd_naive', 'ymu_gpml', 'ystd_gpml', ...
        'hyp_naive', 'hyp_fixed', 'baseline', 'data_cfg');
    writetable(T_compare, fullfile(out_dir, [fstem '.csv']));

    %% Plot
    if make_plots
        figure('Name', 'Naive vs GPML Fixed Prior Mean');
        plot(t_gt, y_gt, 'k-', 'LineWidth', 2, 'DisplayName', 'Ground truth');
        hold on;
        plot(t_obs, y_obs, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');

        if gpml_ok
            plot(t_gt, ymu_naive, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Naive (GPML)');
            ylo_n = ymu_naive - 1.96 * ystd_naive;
            yhi_n = ymu_naive + 1.96 * ystd_naive;
            fill([t_gt; flip(t_gt)], [ylo_n; flip(yhi_n)], 'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

            plot(t_gt, ymu_gpml, 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('Fixed mean (b=%.0f)', baseline));
            ylo_g = ymu_gpml - 1.96 * ystd_gpml;
            yhi_g = ymu_gpml + 1.96 * ystd_gpml;
            fill([t_gt; flip(t_gt)], [ylo_g; flip(yhi_g)], 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end

        hold off;
        xlabel('t');
        ylabel('f(t)');
        title(sprintf('Acute Transient: Naive vs Fixed Mean (both GPML), b=%.0f', baseline));
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

function nlml = gpml_nlz_if_ok(ok, hyp, x, y)
    if ~ok || isempty(hyp)
        nlml = NaN;
        return;
    end
    try
        [nlZ, ~] = gp(hyp, @infGaussLik, @meanConst, @covSEiso, @likGauss, x, y);
        nlml = nlZ;
    catch
        nlml = NaN;
    end
end
