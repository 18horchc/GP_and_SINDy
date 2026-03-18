function [T_compare, hyp_naive, hyp_deriv, hyp_proper] = run_naive_vs_deriv_encoding(opts)
%% RUN_NAIVE_VS_DERIV_ENCODING  Compare naive GP vs virtual-point vs proper derivative obs.
%
%   [T_compare, hyp_naive, hyp_deriv, hyp_proper] = run_naive_vs_deriv_encoding()
%   [T_compare, hyp_naive, hyp_deriv, hyp_proper] = run_naive_vs_deriv_encoding(opts)
%
%   Encodes: "slope decreases from day 2 to day 10".
%
%   Models:
%     1. Naive: GP on function observations only.
%     2. Virtual-point: GPML with extra f-value points near t=2,10 to imply slope.
%     3. Proper deriv: fit_gp_deriv_obs with true derivative observations.
%
%   opts: .baseline (5), .rng_seed (42), .make_plots (true), .out_dir, .gpml_path
%
%   Requires: GPML toolbox for (1),(2). (3) is standalone.
%
%   See also: fit_gp_gpml_naive, fit_gp_deriv_obs, run_naive_vs_pseudo_obs

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

    %% Derivative encoding: slope decreases from day 2 to day 10
    % Use prior f(t)=b+A*t*exp(-c*t) for plausible values
    A = 20; c = 0.6;
    t_deriv = [1.5; 2.5; 9.5; 10.5];
    y_deriv = baseline + A * t_deriv .* exp(-c * t_deriv);
    % Above gives: t=1.5,2.5 ascending (peak ~1.67); t=9.5,10.5 descending (both near baseline)
    % For stronger "slope at 2 > slope at 10", ensure t=2 region has shallower slope.
    % Prior values: f(1.5)~17, f(2.5)~16 (slight decline); f(9.5)~5.7, f(10.5)~5.4 (flat).
    % To get slope(2) > slope(10): make t=10 region steeper. Use:
    y_deriv(3) = 8;   % (9.5, 8) - above baseline
    y_deriv(4) = 5.5; % (10.5, 5.5) - back toward baseline
    % So slope near t=10 is (5.5-8)/1 = -2.5, steeper than t=2 region.

    t_aug = [t_obs; t_deriv];
    y_aug = [y_obs; y_deriv];

    %% Proper derivative observations: x_df=[2;10], d_obs=slope values
    x_df = [2; 10];
    d_obs = [0.5; -0.5];  % positive at t=2, negative at t=10 -> slope decreases

    %% Model 1 & 2: GPML (naive, virtual-point)
    try
        rng(rng_seed);
        [ymu_naive, ys2_naive, hyp_naive] = fit_gp_gpml_naive(t_obs, y_obs, t_gt, gpml_path);
        ystd_naive = sqrt(ys2_naive);

        [ymu_deriv, ys2_deriv, hyp_deriv] = fit_gp_gpml_naive(t_aug, y_aug, t_gt, gpml_path);
        ystd_deriv = sqrt(ys2_deriv);
        gpml_ok = true;
    catch me
        warning('GPML failed: %s. Install from https://gitlab.com/hnickisch/gpml-matlab', me.message);
        ymu_naive = nan(size(t_gt)); ystd_naive = nan(size(t_gt)); hyp_naive = [];
        ymu_deriv = nan(size(t_gt)); ystd_deriv = nan(size(t_gt)); hyp_deriv = [];
        gpml_ok = false;
    end

    %% Model 3: proper derivative observations (standalone)
    rng(rng_seed);
    try
        [ymu_proper, ys2_proper, hyp_proper] = fit_gp_deriv_obs(t_obs, y_obs, x_df, d_obs, t_gt, [0;0;log(0.1)], baseline);
        ystd_proper = sqrt(ys2_proper);
    catch me
        warning('fit_gp_deriv_obs failed: %s', me.message);
        ymu_proper = nan(size(t_gt)); ystd_proper = nan(size(t_gt)); hyp_proper = [];
    end

    %% Metrics
    [rmse_n, mae_n, r2_n] = metric_helpers('point_estimate', y_gt, ymu_naive);
    [nlpd_n, msll_n, crps_n] = metric_helpers('probabilistic', y_gt, ymu_naive, ystd_naive, y_obs);
    [smse_n, cov_n, ~] = metric_helpers('calibration', y_gt, ymu_naive, ystd_naive, y_obs, []);
    nlml_n = gpml_nlz_if_ok(gpml_ok, hyp_naive, t_obs, y_obs);

    if gpml_ok
        [rmse_d, mae_d, r2_d] = metric_helpers('point_estimate', y_gt, ymu_deriv);
        [nlpd_d, msll_d, crps_d] = metric_helpers('probabilistic', y_gt, ymu_deriv, ystd_deriv, y_obs);
        [smse_d, cov_d, ~] = metric_helpers('calibration', y_gt, ymu_deriv, ystd_deriv, y_obs, []);
        nlml_d = gpml_nlz_if_ok(true, hyp_deriv, t_aug, y_aug);
    else
        [rmse_d, mae_d, r2_d, nlpd_d, msll_d, crps_d, smse_d, cov_d, nlml_d] = deal(NaN);
    end

    [rmse_p, mae_p, r2_p] = metric_helpers('point_estimate', y_gt, ymu_proper);
    [nlpd_p, msll_p, crps_p] = metric_helpers('probabilistic', y_gt, ymu_proper, ystd_proper, y_obs);
    [smse_p, cov_p, ~] = metric_helpers('calibration', y_gt, ymu_proper, ystd_proper, y_obs, []);

    T_compare = table(...
        {'naive'; 'deriv_encoding'; 'proper_deriv_obs'}, ...
        [rmse_n; rmse_d; rmse_p], [mae_n; mae_d; mae_p], [r2_n; r2_d; r2_p], ...
        [nlpd_n; nlpd_d; nlpd_p], [msll_n; msll_d; msll_p], [crps_n; crps_d; crps_p], ...
        [smse_n; smse_d; smse_p], [cov_n; cov_d; cov_p], [nlml_n; nlml_d; NaN], ...
        'VariableNames', {'Model','RMSE','MAE','R2','NLPD','MSLL','CRPS','sMSE','Coverage','NLML'});

    fprintf('\nNaive vs Derivative Encoding vs Proper Deriv Obs (slope decreases day 2->10):\n');
    disp(T_compare);

    %% Save
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    fstem = 'results_naive_vs_deriv_encoding';
    save(fullfile(out_dir, [fstem '.mat']), 'T_compare', 't_gt', 'y_gt', 't_obs', 'y_obs', ...
        't_deriv', 'y_deriv', 'x_df', 'd_obs', ...
        'ymu_naive', 'ystd_naive', 'ymu_deriv', 'ystd_deriv', 'ymu_proper', 'ystd_proper', ...
        'hyp_naive', 'hyp_deriv', 'hyp_proper', 'baseline', 'data_cfg');
    writetable(T_compare, fullfile(out_dir, [fstem '.csv']));

    %% Plot
    if make_plots
        figure('Name', 'Naive vs Derivative Encoding');
        plot(t_gt, y_gt, 'k-', 'LineWidth', 2, 'DisplayName', 'Ground truth');
        hold on;
        plot(t_obs, y_obs, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
        plot(t_deriv, y_deriv, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g', ...
            'DisplayName', 'Deriv guide (t=1.5,2.5,9.5,10.5)');
        plot(x_df, baseline + zeros(size(x_df)), 'c^', 'MarkerSize', 8, 'MarkerFaceColor', 'c', ...
            'DisplayName', 'Proper deriv obs (t=2,10)');

        if gpml_ok
            plot(t_gt, ymu_naive, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Naive');
            ylo_n = ymu_naive - 1.96 * ystd_naive;
            yhi_n = ymu_naive + 1.96 * ystd_naive;
            fill([t_gt; flip(t_gt)], [ylo_n; flip(yhi_n)], 'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

            plot(t_gt, ymu_deriv, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Deriv encoding');
            ylo_d = ymu_deriv - 1.96 * ystd_deriv;
            yhi_d = ymu_deriv + 1.96 * ystd_deriv;
            fill([t_gt; flip(t_gt)], [ylo_d; flip(yhi_d)], 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        plot(t_gt, ymu_proper, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Proper deriv obs');
        ylo_p = ymu_proper - 1.96 * ystd_proper;
        yhi_p = ymu_proper + 1.96 * ystd_proper;
        fill([t_gt; flip(t_gt)], [ylo_p; flip(yhi_p)], 'm', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

        hold off;
        xlabel('t');
        ylabel('f(t)');
        title('Acute Transient: Naive vs Derivative Encoding (slope \downarrow from day 2 to 10)');
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
