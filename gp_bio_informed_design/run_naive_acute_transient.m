function [T_results, gpr] = run_naive_acute_transient(opts)
%% RUN_NAIVE_ACUTE_TRANSIENT  Naive GP on acute transient.
%
%   [T_results, gpr] = run_naive_acute_transient()
%   [T_results, gpr] = run_naive_acute_transient(opts)
%
%   Data: t = 0,1,2,3,5,7,14 (7 points), 2 replicates each, 5% noise. Domain [0,15].
%   GP: Squared Exponential kernel, constant mean (naive—no biological encoding).
%
%   opts: .rng_seed (42), .make_plots (true), .out_dir (results folder)
%
%   Returns: T_results (metrics table), gpr (fitted RegressionGP model).
%
%   See also: ground_truth_bio, generate_bio_data, metric_helpers

    if nargin < 1, opts = struct(); end

    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    addpath(fullfile(script_dir, '..'));

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

    %% Naive GP: SE kernel, constant mean
    rng(rng_seed);  % Reproducibility: reset before fitrgp (its optimization uses randomness)
    gpr = fitrgp(t_obs, y_obs, 'KernelFunction', 'squaredexponential', ...
        'BasisFunction', 'constant', 'FitMethod', 'exact', 'PredictMethod', 'exact', ...
        'Standardize', true, 'OptimizeHyperparameters', 'auto', ...
        'HyperparameterOptimizationOptions', struct('ShowPlots', false, 'Verbose', 0, ...
        'UseParallel', false));  % UseParallel=false for deterministic hyperparameter search

    [ymu, ystd, ~] = predict(gpr, t_gt);

    %% Metrics
    [rmse, mae, r2] = metric_helpers('point_estimate', y_gt, ymu);
    [nlpd, msll, crps] = metric_helpers('probabilistic', y_gt, ymu, ystd, y_obs);
    [smse, coverage, nlml] = metric_helpers('calibration', y_gt, ymu, ystd, y_obs, gpr);

    T_results = table(rmse, mae, r2, nlpd, msll, crps, smse, coverage, nlml, ...
        'VariableNames', {'RMSE','MAE','R2','NLPD','MSLL','CRPS','sMSE','Coverage','NLML'});

    fprintf('\nNaive GP (acute transient, SE kernel):\n');
    disp(T_results);

    %% Save
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    fstem = 'results_naive_acute_transient_t01235714_2rep_5noise';
    save(fullfile(out_dir, [fstem '.mat']), 'T_results', 't_gt', 'y_gt', 't_obs', 'y_obs', ...
        'ymu', 'ystd', 'gpr', 'data_cfg');
    writetable(T_results, fullfile(out_dir, [fstem '.csv']));

    %% Plot with uncertainty bounds
    if make_plots
        figure('Name', 'Naive GP: Acute Transient');
        plot(t_gt, y_gt, 'k-', 'LineWidth', 2, 'DisplayName', 'Ground truth');
        hold on;
        plot(t_obs, y_obs, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
        plot(t_gt, ymu, 'b-', 'LineWidth', 1.5, 'DisplayName', 'GP mean (naive)');

        ylo = ymu - 1.96 * ystd;
        yhi = ymu + 1.96 * ystd;
        fill([t_gt; flip(t_gt)], [ylo; flip(yhi)], 'b', 'FaceAlpha', 0.2, ...
            'EdgeColor', 'none', 'DisplayName', '95% CI');

        hold off;
        xlabel('t');
        ylabel('f(t)');
        title('Acute Transient: Naive GP (SE kernel), t=0,1,2,3,5,7,14, 2 rep, 5% noise');
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
