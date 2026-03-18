function [T_compare, hyp_naive, hyp_pseudo] = run_naive_vs_pseudo_obs(opts)
%% RUN_NAIVE_VS_PSEUDO_OBS  Compare naive GP vs pseudo-observations (baseline around day 11).
%
%   [T_compare, hyp_naive, hyp_pseudo] = run_naive_vs_pseudo_obs()
%   [T_compare, hyp_naive, hyp_pseudo] = run_naive_vs_pseudo_obs(opts)
%
%   Data: t = 0,1,2,3,5,7,14, 2 rep, 5% noise. Domain [0,15].
%   Naive: GPML on (t_obs, y_obs) only.
%   Pseudo-obs: Same GPML, but augment with (10, b), (11, b), (12, b) to encode
%   "trajectory returns to baseline around day 11".
%
%   opts: .baseline (5), .rng_seed (42), .make_plots (true), .out_dir, .gpml_path
%
%   Requires: GPML toolbox.
%
%   See also: fit_gp_gpml_naive, run_naive_vs_gpml_fixed_mean

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

    %% Pseudo-obs: trajectory returns to baseline around day 11
    t_pseudo = [10; 11; 12];
    y_pseudo = baseline * ones(3, 1);
    t_aug = [t_obs; t_pseudo];
    y_aug = [y_obs; y_pseudo];

    %% Both models use GPML
    try
        rng(rng_seed);
        [ymu_naive, ys2_naive, hyp_naive] = fit_gp_gpml_naive(t_obs, y_obs, t_gt, gpml_path);
        ystd_naive = sqrt(ys2_naive);

        [ymu_pseudo, ys2_pseudo, hyp_pseudo] = fit_gp_gpml_naive(t_aug, y_aug, t_gt, gpml_path);
        ystd_pseudo = sqrt(ys2_pseudo);
        gpml_ok = true;
    catch me
        warning('GPML failed: %s. Install from https://gitlab.com/hnickisch/gpml-matlab', me.message);
        ymu_naive = nan(size(t_gt)); ystd_naive = nan(size(t_gt)); hyp_naive = [];
        ymu_pseudo = nan(size(t_gt)); ystd_pseudo = nan(size(t_gt)); hyp_pseudo = [];
        gpml_ok = false;
    end

    %% Metrics (evaluate on ground truth; y_obs for probabilistic baseline)
    [rmse_n, mae_n, r2_n] = metric_helpers('point_estimate', y_gt, ymu_naive);
    [nlpd_n, msll_n, crps_n] = metric_helpers('probabilistic', y_gt, ymu_naive, ystd_naive, y_obs);
    [smse_n, cov_n, ~] = metric_helpers('calibration', y_gt, ymu_naive, ystd_naive, y_obs, []);
    nlml_n = gpml_nlz_if_ok(gpml_ok, hyp_naive, t_obs, y_obs);

    if gpml_ok
        [rmse_p, mae_p, r2_p] = metric_helpers('point_estimate', y_gt, ymu_pseudo);
        [nlpd_p, msll_p, crps_p] = metric_helpers('probabilistic', y_gt, ymu_pseudo, ystd_pseudo, y_obs);
        [smse_p, cov_p, ~] = metric_helpers('calibration', y_gt, ymu_pseudo, ystd_pseudo, y_obs, []);
        nlml_p = gpml_nlz_if_ok(true, hyp_pseudo, t_aug, y_aug);
    else
        [rmse_p, mae_p, r2_p, nlpd_p, msll_p, crps_p, smse_p, cov_p, nlml_p] = deal(NaN);
    end

    T_compare = table(...
        {'naive'; 'pseudo_obs'}, ...
        [rmse_n; rmse_p], [mae_n; mae_p], [r2_n; r2_p], ...
        [nlpd_n; nlpd_p], [msll_n; msll_p], [crps_n; crps_p], ...
        [smse_n; smse_p], [cov_n; cov_p], [nlml_n; nlml_p], ...
        'VariableNames', {'Model','RMSE','MAE','R2','NLPD','MSLL','CRPS','sMSE','Coverage','NLML'});

    fprintf('\nNaive vs Pseudo-Obs (baseline ~day 11, b=%.1f):\n', baseline);
    disp(T_compare);

    %% Save
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    fstem = 'results_naive_vs_pseudo_obs_day11';
    save(fullfile(out_dir, [fstem '.mat']), 'T_compare', 't_gt', 'y_gt', 't_obs', 'y_obs', ...
        't_pseudo', 'y_pseudo', 'ymu_naive', 'ystd_naive', 'ymu_pseudo', 'ystd_pseudo', ...
        'hyp_naive', 'hyp_pseudo', 'baseline', 'data_cfg');
    writetable(T_compare, fullfile(out_dir, [fstem '.csv']));

    %% Plot
    if make_plots
        figure('Name', 'Naive vs Pseudo-Observations');
        plot(t_gt, y_gt, 'k-', 'LineWidth', 2, 'DisplayName', 'Ground truth');
        hold on;
        plot(t_obs, y_obs, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
        plot(t_pseudo, y_pseudo, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', 'Pseudo-obs (t=9,10,11)');

        if gpml_ok
            plot(t_gt, ymu_naive, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Naive');
            ylo_n = ymu_naive - 1.96 * ystd_naive;
            yhi_n = ymu_naive + 1.96 * ystd_naive;
            fill([t_gt; flip(t_gt)], [ylo_n; flip(yhi_n)], 'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

            plot(t_gt, ymu_pseudo, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Pseudo-obs');
            ylo_p = ymu_pseudo - 1.96 * ystd_pseudo;
            yhi_p = ymu_pseudo + 1.96 * ystd_pseudo;
            fill([t_gt; flip(t_gt)], [ylo_p; flip(yhi_p)], 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end

        hold off;
        xlabel('t');
        ylabel('f(t)');
        title('Acute Transient: Naive vs Pseudo-Obs (returns to baseline ~day 11)');
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
