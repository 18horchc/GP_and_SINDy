%% RUN_ALL_LOGISTIC_EXPERIMENTS  Run all GP logistic experiments and aggregate metrics.
%
%   run_all_logistic_experiments
%   run_all_logistic_experiments(opts)
%
%   Runs all 20 experiment configurations (2 sampling x 5 noise x 2 optimization)
%   via run_logistic_experiment, then aggregate_logistic_metrics.
%
%   opts (optional): struct with
%     .make_plots   - true/false (default: true). Set false for batch/headless.
%     .run_subset   - indices to run, e.g. [1 5 10] or 1:20 (default: 1:20).
%     .do_aggregate - run aggregate_logistic_metrics at end (default: true).
%     .out_dir      - where to save (default: gp_logistic_design folder).
%
%   Experiment configs (same as original log_exp_01..20):
%     Exp 01-05:  Regular, 0/1/5/10/20% noise, default optimization
%     Exp 06-10:  Irregular, 0/1/5/10/20% noise, default optimization
%     Exp 11-15:  Regular, 0/1/5/10/20% noise, auto hyperparameter optimization
%     Exp 16-20:  Irregular, 0/1/5/10/20% noise, auto optimization

function run_all_logistic_experiments(opts)
    if nargin < 1, opts = struct(); end
    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    addpath(fullfile(script_dir, '..'));

    make_plots = get_opt(opts, 'make_plots', true);
    run_subset = get_opt(opts, 'run_subset', 1:20);
    do_aggregate = get_opt(opts, 'do_aggregate', true);
    out_dir = get_opt(opts, 'out_dir', script_dir);

    %% Build experiment configs: {exp_id, noise_pct, is_regular, is_auto}
    configs = build_logistic_configs();

    %% Run selected experiments
    n_total = length(configs);
    for ix = run_subset
        if ix < 1 || ix > n_total
            warning('Skipping invalid index %d (valid: 1..%d)', ix, n_total);
            continue
        end
        cfg = configs{ix};
        fprintf('Running experiment %d/%d: exp_%02d %s %.0f%% noise %s\n', ...
            ix, n_total, cfg.exp_id, ...
            iif(cfg.is_regular, 'regular', 'irregular'), ...
            cfg.noise_pct * 100, ...
            iif(cfg.is_auto, 'auto', 'default'));
        run_logistic_experiment(cfg, struct('out_dir', out_dir, 'make_plots', make_plots));
    end

    %% Aggregate
    if do_aggregate
        fprintf('Running aggregate_logistic_metrics.m\n');
        run(fullfile(script_dir, 'aggregate_logistic_metrics.m'));
    end
    fprintf('Done.\n');
end

function configs = build_logistic_configs()
    % Same 20 experiments as original log_exp_01..20
    noise_levels = [0, 0.01, 0.05, 0.10, 0.20];
    configs = {};
    exp_id = 0;
    % 01-05: regular, default
    for n = noise_levels
        exp_id = exp_id + 1;
        configs{end+1} = struct('exp_id', exp_id, 'noise_pct', n, 'is_regular', true, 'is_auto', false);
    end
    % 06-10: irregular, default
    for n = noise_levels
        exp_id = exp_id + 1;
        configs{end+1} = struct('exp_id', exp_id, 'noise_pct', n, 'is_regular', false, 'is_auto', false);
    end
    % 11-15: regular, auto
    for n = noise_levels
        exp_id = exp_id + 1;
        configs{end+1} = struct('exp_id', exp_id, 'noise_pct', n, 'is_regular', true, 'is_auto', true);
    end
    % 16-20: irregular, auto
    for n = noise_levels
        exp_id = exp_id + 1;
        configs{end+1} = struct('exp_id', exp_id, 'noise_pct', n, 'is_regular', false, 'is_auto', true);
    end
end

function v = get_opt(opts, name, default)
    if isfield(opts, name), v = opts.(name); else, v = default; end
end

function v = iif(cond, a, b)
    if cond, v = a; else, v = b; end
end
