function [T_compare, hyp_naive, hyp_mono] = run_monotonic_ep_acute_transient(opts)
%% RUN_MONOTONIC_EP_ACUTE_TRANSIENT  Riihimäki–Vehtari monotonicity on acute transient.
%
%   Enforces NEGATIVE derivative on [2, 10] using EP with probit likelihood.
%
%   [T_compare, hyp_naive, hyp_mono] = run_monotonic_ep_acute_transient()
%   [T_compare, hyp_naive, hyp_mono] = run_monotonic_ep_acute_transient(opts)
%
%   INPUT (opts):
%     .Tm        - virtual derivative locations (default: [2, 3.5, 5, 6.5, 8, 9.5])
%     .nu        - monotonicity softness (default 1e-6)
%     .rng_seed  - RNG seed (42)
%     .make_plots - plot results (true)
%     .out_dir   - output directory
%     .gpml_path - for naive GP (optional)
%     .exclude_t - times to exclude from observations, e.g. 5 or [5] (default [])
%
%   See: Riihimäki & Vehtari (2010) "Gaussian processes with monotonicity
%        information", fit_gp_monotonic_ep.

    if nargin < 1, opts = struct(); end

    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    addpath(fullfile(script_dir, '..'));

    Tm        = get_opt(opts, 'Tm', [2; 3.5; 5; 6.5; 8; 9.5]);
    nu        = get_opt(opts, 'nu', 1e-6);
    rng_seed  = get_opt(opts, 'rng_seed', 42);
    make_plots = get_opt(opts, 'make_plots', true);
    out_dir   = get_opt(opts, 'out_dir', fullfile(script_dir, 'results'));
    gpml_path = get_opt(opts, 'gpml_path', []);
    exclude_t = get_opt(opts, 'exclude_t', []);

    rng(rng_seed);

    %% 1. Data (irregular T in [0,15])
    t_points_default = [0; 1; 2; 3; 5; 7; 14];
    if ~isempty(exclude_t)
        t_points_default = t_points_default(~ismember(t_points_default, exclude_t(:)));
    end
    data_cfg = struct('t_points', t_points_default, 'n_replicates', 2, ...
        'noise_pct', 0.05, 't_domain', [0, 15]);
    params = struct('t_domain', [0, 15]);
    [t_gt, y_gt, ~] = ground_truth_bio('acute_transient', 500, params);
    [t_obs, y_obs] = generate_bio_data(t_gt, y_gt, data_cfg, rng_seed);

    %% Naive GP (baseline)
    try
        [ymu_naive, ys2_naive, hyp_naive] = fit_gp_gpml_naive(t_obs, y_obs, t_gt, gpml_path);
        ystd_naive = sqrt(ys2_naive);
        naive_ok = true;
    catch me
        warning('GPML (naive) failed: %s', me.message);
        ymu_naive = nan(size(t_gt)); ystd_naive = nan(size(t_gt)); hyp_naive = [];
        naive_ok = false;
    end

    %% Monotonic EP (negative derivative on [2,10])
    ep_opts = struct('Tm', Tm, 'nu', nu);
    [ymu_mono, ys2_mono, hyp_mono, post_d_mean, post_d_var, mu_site, sigma2_site] = ...
        fit_gp_monotonic_ep(t_obs, y_obs, [2, 10], t_gt, ep_opts);
    ystd_mono = sqrt(ys2_mono);

    %% Diagnostic: P(f'(t) > 0) on dense grid in [2,10]
    t_check = linspace(2, 10, 100)';
    [d_mean_check, d_var_check] = predict_deriv_monotonic_ep(...
        t_obs, y_obs, Tm, hyp_mono, mu_site, sigma2_site, t_check);
    prob_positive = 1 - normcdf(0, d_mean_check, sqrt(d_var_check));
    max_prob_pos = max(prob_positive);

    %% Metrics
    [rmse_n, mae_n, r2_n] = metric_helpers('point_estimate', y_gt, ymu_naive);
    [nlpd_n, msll_n, crps_n] = metric_helpers('probabilistic', y_gt, ymu_naive, ystd_naive, y_obs);
    [smse_n, cov_n, ~] = metric_helpers('calibration', y_gt, ymu_naive, ystd_naive, y_obs, []);

    [rmse_m, mae_m, r2_m] = metric_helpers('point_estimate', y_gt, ymu_mono);
    [nlpd_m, msll_m, crps_m] = metric_helpers('probabilistic', y_gt, ymu_mono, ystd_mono, y_obs);
    [smse_m, cov_m, ~] = metric_helpers('calibration', y_gt, ymu_mono, ystd_mono, y_obs, []);

    if ~naive_ok
        [rmse_n, mae_n, r2_n, nlpd_n, msll_n, crps_n, smse_n, cov_n] = deal(NaN);
    end

    T_compare = table(...
        {'naive'; 'monotonic_ep'}, ...
        [rmse_n; rmse_m], [mae_n; mae_m], [r2_n; r2_m], ...
        [nlpd_n; nlpd_m], [msll_n; msll_m], [crps_n; crps_m], ...
        [smse_n; smse_m], [cov_n; cov_m], ...
        'VariableNames', {'Model','RMSE','MAE','R2','NLPD','MSLL','CRPS','sMSE','Coverage'});
    T_compare.Properties.RowNames = {'naive','monotonic_ep'};

    if isempty(exclude_t)
        fprintf('\nAcute transient: Naive vs Monotonic EP (negative deriv on [2,10]):\n');
    else
        fprintf('\nAcute transient (exclude t=%s): Naive vs Monotonic EP (negative deriv on [2,10]):\n', mat2str(exclude_t(:)'));
    end
    disp(T_compare);
    fprintf('Diagnostic: max P(f''>0) in [2,10] = %.4f\n', max_prob_pos);

    %% Save
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    if isempty(exclude_t)
        fstem = 'results_monotonic_ep_acute_transient';
    else
        fstem = sprintf('results_monotonic_ep_acute_transient_exclude_t%s', strrep(mat2str(exclude_t(:)'), ' ', '_'));
    end
    save(fullfile(out_dir, [fstem '.mat']), 'T_compare', 't_gt', 'y_gt', 't_obs', 'y_obs', ...
        'Tm', 'nu', 'exclude_t', 'ymu_naive', 'ystd_naive', 'ymu_mono', 'ystd_mono', ...
        'hyp_naive', 'hyp_mono', 'post_d_mean', 'post_d_var', ...
        'prob_positive', 't_check', 'max_prob_pos');
    writetable(T_compare, fullfile(out_dir, [fstem '.csv']));

    %% Plot
    if make_plots
        figure('Name', 'Monotonic EP: Acute Transient');
        plot(t_gt, y_gt, 'k-', 'LineWidth', 2, 'DisplayName', 'Ground truth');
        hold on;
        plot(t_obs, y_obs, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
        if naive_ok
            plot(t_gt, ymu_naive, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Naive');
            fill([t_gt; flip(t_gt)], [ymu_naive-1.96*ystd_naive; flip(ymu_naive+1.96*ystd_naive)], ...
                'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        plot(t_gt, ymu_mono, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Monotonic EP');
        fill([t_gt; flip(t_gt)], [ymu_mono-1.96*ystd_mono; flip(ymu_mono+1.96*ystd_mono)], ...
            'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        ylims = ylim;
        plot(Tm, repmat(ylims(1), size(Tm)), 'c|', 'MarkerSize', 12, 'DisplayName', 'Virtual deriv pts');
        hold off;
        xlabel('t'); ylabel('f(t)');
        excl_str = iif_empty(exclude_t, '');
        title(['Acute Transient' excl_str ': Monotonic EP (f''<0 on [2,10])']);
        legend('Location', 'best');
        grid on; xlim([0 15]);

        % Derivative diagnostic
        figure('Name', 'Derivative diagnostic');
        plot(t_check, d_mean_check, 'r-', 'LineWidth', 1.5, 'DisplayName', 'E[f''(t)]');
        hold on;
        d_lo = d_mean_check - 1.96 * sqrt(d_var_check);
        d_hi = d_mean_check + 1.96 * sqrt(d_var_check);
        fill([t_check; flip(t_check)], [d_lo; flip(d_hi)], 'r', 'FaceAlpha', 0.2, ...
            'EdgeColor', 'none', 'DisplayName', '95% CI');
        plot([2 10], [0 0], 'k--');
        xlabel('t'); ylabel('f''(t)');
        title('Posterior derivative (should be negative in [2,10])');
        legend('Location', 'best'); grid on; xlim([2 10]);
    end
end

function s = iif_empty(exclude_t, default)
    if isempty(exclude_t)
        s = default;
    else
        s = [' (exclude t=' mat2str(exclude_t(:)') ')'];
    end
end

function v = get_opt(opts, name, default)
    if isempty(opts) || ~isfield(opts, name)
        v = default;
    else
        v = opts.(name);
    end
end
