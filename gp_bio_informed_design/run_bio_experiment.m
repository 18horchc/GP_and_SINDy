function [T_results, out_path] = run_bio_experiment(toy_problem, data_cfg, encoding_ids, opts)
%% RUN_BIO_EXPERIMENT  Run bio-informed GP experiment for acute transient (and future toy problems).
%
%   [T_results, out_path] = run_bio_experiment('acute_transient')
%   [T_results, out_path] = run_bio_experiment('acute_transient', data_cfg)
%   [T_results, out_path] = run_bio_experiment('acute_transient', data_cfg, encoding_ids)
%   [T_results, out_path] = run_bio_experiment('acute_transient', data_cfg, encoding_ids, opts)
%
%   toy_problem:   'acute_transient' (domain t in [0, 20])
%   data_cfg:      struct with n_points (10), n_replicates (2), noise_pct (0.05),
%                  irregular (true), t_domain ([0,20]). Default: 10 irregular, 2 rep, 5% noise.
%   encoding_ids:  cell array of encoding IDs, or 'all' for all. Default: 6 key configs.
%   opts:          .out_dir, .make_plots, .rng_seed (default 42)
%
%   Output: T_results (table), out_path (where .mat/.csv saved).
%
%   See also: encoding_registry, fit_gp_bio, metric_helpers

    if nargin < 2, data_cfg = struct(); end
    if nargin < 3 || isempty(encoding_ids)
        encoding_ids = {'vanilla','mean_only','matern_only','pseudo_only','deriv_only','full'};
    end
    if nargin < 4, opts = struct(); end

    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    addpath(fullfile(script_dir, '..'));

    out_dir   = get_opt(opts, 'out_dir', fullfile(script_dir, 'results'));
    make_plots = get_opt(opts, 'make_plots', true);
    rng_seed  = get_opt(opts, 'rng_seed', 42);

    rng(rng_seed);

    %% Default data config: 10 irregular, 2 rep, 5% noise
    data_cfg = merge_data_cfg(data_cfg);

    %% Ground truth
    [t_gt, y_gt, ~] = ground_truth_bio(toy_problem, 500);

    %% Generate observations
    [t_obs, y_obs] = generate_bio_data(t_gt, y_gt, data_cfg, rng_seed);

    %% Encoding configs
    if ischar(encoding_ids) || isstring(encoding_ids)
        if strcmpi(encoding_ids, 'all')
            configs = encoding_registry();
        else
            configs = {encoding_registry(encoding_ids)};
        end
    else
        configs = cell(1, length(encoding_ids));
        for i = 1:length(encoding_ids)
            configs{i} = encoding_registry(encoding_ids{i});
        end
    end

    %% Fit each encoding, compute metrics
    results = [];
    pred_ymu = cell(length(configs), 1);
    pred_ystd = cell(length(configs), 1);

    for ic = 1:length(configs)
        enc = configs{ic};
        [ymu, ystd, gpr] = fit_gp_bio(t_obs, y_obs, t_gt, enc);

        pred_ymu{ic} = ymu;
        pred_ystd{ic} = ystd;

        [rmse, mae, r2] = metric_helpers('point_estimate', y_gt, ymu);
        [nlpd, msll, crps] = metric_helpers('probabilistic', y_gt, ymu, ystd, y_obs);
        [smse, coverage, nlml] = metric_helpers('calibration', y_gt, ymu, ystd, y_obs, gpr);

        results = [results; {enc.id, enc.label, rmse, mae, r2, nlpd, msll, crps, smse, coverage, nlml}];
    end

    T_results = cell2table(results, 'VariableNames', ...
        {'Encoding_ID','Label','RMSE','MAE','R2','NLPD','MSLL','CRPS','sMSE','Coverage','NLML'});

    %% Save
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    fstem = get_output_fstem(toy_problem, data_cfg);
    out_path = fullfile(out_dir, [fstem '.mat']);
    save(out_path, 'T_results', 't_gt', 'y_gt', 't_obs', 'y_obs', 'configs', 'pred_ymu', 'pred_ystd', 'data_cfg');
    writetable(T_results, fullfile(out_dir, [fstem '.csv']));

    %% Plots
    if make_plots
        figure('Name', 'Bio-Informed GP: Acute Transient');
        plot(t_gt, y_gt, 'k-', 'LineWidth', 2, 'DisplayName', 'Ground truth');
        hold on;
        plot(t_obs, y_obs, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Observed');

        colors = lines(length(configs));
        for ic = 1:length(configs)
            ymu = pred_ymu{ic};
            ystd = pred_ystd{ic};
            if ~all(isnan(ymu))
                plot(t_gt, ymu, '-', 'LineWidth', 1.2, 'Color', colors(ic,:), ...
                    'DisplayName', configs{ic}.label);
            end
        end
        hold off; xlabel('t'); ylabel('f(t)');
        title(sprintf('Acute Transient: N=%d, %d rep, %.0f%% noise', ...
            data_cfg.n_points, data_cfg.n_replicates, data_cfg.noise_pct*100));
        legend('Location', 'best', 'FontSize', 8);
        grid on; xlim([0 20]);
    end

    fprintf('\nBio-informed GP results (%s):\n', fstem);
    disp(T_results);
end

function cfg = merge_data_cfg(cfg)
    if isempty(cfg), cfg = struct(); end
    d = struct('n_points', 10, 'n_replicates', 2, 'noise_pct', 0.05, ...
        'irregular', true, 't_domain', [0, 20]);
    f = fieldnames(d);
    for i = 1:length(f)
        if ~isfield(cfg, f{i}) || isempty(cfg.(f{i}))
            cfg.(f{i}) = d.(f{i});
        end
    end
end

function fstem = get_output_fstem(toy_problem, data_cfg)
    n = data_cfg.n_points;
    rep = data_cfg.n_replicates;
    noise = round(data_cfg.noise_pct * 100);
    fstem = sprintf('results_bio_%s_N%d_%drep_%dnoise', toy_problem, n, rep, noise);
end

function v = get_opt(opts, name, default)
    if isempty(opts) || ~isfield(opts, name)
        v = default;
    else
        v = opts.(name);
    end
end
