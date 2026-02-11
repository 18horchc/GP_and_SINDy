%% RUN_ALL_LV_EXPERIMENTS  Run GP Lotka-Volterra experiments (aggregation is separate).
%
%   run_all_LV_experiments
%   run_all_LV_experiments(opts)
%
%   Runs experiment configurations via run_LV_experiment. Does NOT run
%   aggregate_LV_metrics — call that explicitly when you want to build
%   LV_metrics_all.csv/.mat.
%
%   opts (optional): struct with
%     .make_plots    - true/false (default: true). Set false for batch/headless.
%     .run_subset    - indices into the (possibly filtered) config list (default: all).
%     .config_filter - filter configs before running. See examples below.
%     .do_aggregate  - DEPRECATED/ignored. Call aggregate_LV_metrics separately.
%     .out_dir       - where to save (default: gp_LV_design folder).
%
%   config_filter: function handle @(cfg) true/false, or struct for field matching.
%     Examples:
%       opts.config_filter = @(c) c.noise_pct == 0;           % zero noise only
%       opts.config_filter = @(c) c.is_regular;              % regular sampling only
%       opts.config_filter = @(c) ~c.is_auto;                % default optimization only
%       opts.config_filter = @(c) c.noise_pct <= 0.05;       % low noise
%       opts.config_filter = struct('noise_pct', 0);         % match noise_pct==0
%
%   Experiment configs (same as original LV_exp_01..20):
%     Exp 01-05:  Regular, 0/1/5/10/20% noise, default optimization
%     Exp 06-10:  Irregular, 0/1/5/10/20% noise, default optimization
%     Exp 11-15:  Regular, 0/1/5/10/20% noise, auto hyperparameter optimization
%     Exp 16-20:  Irregular, 0/1/5/10/20% noise, auto optimization
%
%   Note: Periodic kernel uses KernelParameters (no auto opt). Six kernels total:
%   5 built-in + periodicKernel. Two states: Prey, Predator.

function run_all_LV_experiments(opts)
    if nargin < 1, opts = struct(); end
    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    addpath(fullfile(script_dir, '..'));

    make_plots = get_opt(opts, 'make_plots', true);
    run_subset = get_opt(opts, 'run_subset', []);
    config_filter = get_opt(opts, 'config_filter', []);
    out_dir = get_opt(opts, 'out_dir', script_dir);

    %% Build experiment configs: {exp_id, noise_pct, is_regular, is_auto}
    configs = build_LV_configs();

    %% Apply config_filter if provided
    if ~isempty(config_filter)
        configs = apply_config_filter(configs, config_filter);
    end

    if isempty(run_subset)
        run_subset = 1:length(configs);
    end

    %% Run selected experiments (indices into the config list, which may be filtered)
    n_total = length(configs);
    for ii = 1:length(run_subset)
        ix = run_subset(ii);
        if ix < 1 || ix > n_total
            warning('Skipping invalid index %d (valid: 1..%d)', ix, n_total);
            continue
        end
        cfg = configs{ix};
        fprintf('Running LV experiment %d/%d: exp_%02d %s %.0f%% noise %s\n', ...
            ii, length(run_subset), cfg.exp_id, ...
            iif(cfg.is_regular, 'regular', 'irregular'), ...
            cfg.noise_pct * 100, ...
            iif(cfg.is_auto, 'auto', 'default'));
        run_LV_experiment(cfg, struct('out_dir', out_dir, 'make_plots', make_plots));
    end

    fprintf('Done. Call aggregate_LV_metrics() to build LV_metrics_all.\n');
end

function filtered = apply_config_filter(configs, filter_spec)
    if isa(filter_spec, 'function_handle')
        keep = cellfun(filter_spec, configs);
        filtered = configs(keep);
    elseif isstruct(filter_spec)
        fnames = fieldnames(filter_spec);
        keep = true(1, length(configs));
        for i = 1:length(configs)
            for j = 1:length(fnames)
                fn = fnames{j};
                if ~isfield(configs{i}, fn) || configs{i}.(fn) ~= filter_spec.(fn)
                    keep(i) = false;
                    break
                end
            end
        end
        filtered = configs(keep);
    else
        filtered = configs;
    end
end

function configs = build_LV_configs()
    noise_levels = [0, 0.01, 0.05, 0.10, 0.20];
    configs = {};
    exp_id = 0;
    for n = noise_levels
        exp_id = exp_id + 1;
        configs{end+1} = struct('exp_id', exp_id, 'noise_pct', n, 'is_regular', true, 'is_auto', false);
    end
    for n = noise_levels
        exp_id = exp_id + 1;
        configs{end+1} = struct('exp_id', exp_id, 'noise_pct', n, 'is_regular', false, 'is_auto', false);
    end
    for n = noise_levels
        exp_id = exp_id + 1;
        configs{end+1} = struct('exp_id', exp_id, 'noise_pct', n, 'is_regular', true, 'is_auto', true);
    end
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
