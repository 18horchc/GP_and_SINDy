%% RUN_ALL_LOGISTIC_EXPERIMENTS  Run GP logistic experiments (aggregation is separate).
%
%   run_all_logistic_experiments
%   run_all_logistic_experiments(opts)
%
%   Runs experiment configurations via run_logistic_experiment. Does NOT run
%   aggregate_logistic_metrics — call that explicitly when you want to build
%   logistic_metrics_all.csv/.mat.
%
%   opts (optional): struct with
%     .make_plots    - true/false (default: true). Set false for batch/headless.
%     .run_subset    - indices into the (possibly filtered) config list (default: all).
%     .config_filter - filter configs before running. See examples below.
%     .do_aggregate  - DEPRECATED/ignored. Call aggregate_logistic_metrics separately.
%     .out_dir       - where to save (default: gp_logistic_design folder).
%
%   config_filter: function handle @(cfg) true/false, or struct for field matching.
%     Examples:
%       opts.config_filter = @(c) c.noise_pct == 0;           % zero noise only
%       opts.config_filter = @(c) c.is_regular;               % regular sampling only
%       opts.config_filter = @(c) ~c.is_auto;                 % default optimization only
%       opts.config_filter = @(c) c.noise_pct <= 0.05;        % low noise
%       opts.config_filter = struct('noise_pct', 0);           % match noise_pct==0
%       opts.config_filter = struct('is_regular', true);      % match is_regular==true
%
%   Experiment configs (60 total):
%     Exp 01-20:  Base (1 obs per time); regular/irregular, 5 noise, default/auto
%     Exp 21-40:  Same conditions, 3 replicates per time point
%     Exp 41-60:  Same conditions, 8 replicates per time point

function run_all_logistic_experiments(opts)
    if nargin < 1, opts = struct(); end
    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    addpath(fullfile(script_dir, '..'));

    make_plots = get_opt(opts, 'make_plots', true);
    run_subset = get_opt(opts, 'run_subset', []);
    config_filter = get_opt(opts, 'config_filter', []);
    out_dir = get_opt(opts, 'out_dir', script_dir);

    %% Build experiment configs: {exp_id, noise_pct, is_regular, is_auto}
    configs = build_logistic_configs();

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
        n_rep = get_opt(cfg, 'n_replicates', 1);
        rep_str = iif(n_rep > 1, sprintf(' %drep', n_rep), '');
        fprintf('Running experiment %d/%d: exp_%02d %s %.0f%% noise %s%s\n', ...
            ii, length(run_subset), cfg.exp_id, ...
            iif(cfg.is_regular, 'regular', 'irregular'), ...
            cfg.noise_pct * 100, ...
            iif(cfg.is_auto, 'auto', 'default'), rep_str);
        run_logistic_experiment(cfg, struct('out_dir', out_dir, 'make_plots', make_plots));
    end

    fprintf('Done. Call aggregate_logistic_metrics() to build logistic_metrics_all.\n');
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

function configs = build_logistic_configs()
    % 60 experiments: 20 base (1 rep) + 20 with 3 reps + 20 with 8 reps
    noise_levels = [0, 0.01, 0.05, 0.10, 0.20];
    configs = {};
    exp_id = 0;
    for n_rep = [1, 3, 8]
        % 01-05 / 21-25 / 41-45: regular, default
        for n = noise_levels
            exp_id = exp_id + 1;
            configs{end+1} = struct('exp_id', exp_id, 'noise_pct', n, 'is_regular', true, 'is_auto', false, 'n_replicates', n_rep);
        end
        % 06-10 / 26-30 / 46-50: irregular, default
        for n = noise_levels
            exp_id = exp_id + 1;
            configs{end+1} = struct('exp_id', exp_id, 'noise_pct', n, 'is_regular', false, 'is_auto', false, 'n_replicates', n_rep);
        end
        % 11-15 / 31-35 / 51-55: regular, auto
        for n = noise_levels
            exp_id = exp_id + 1;
            configs{end+1} = struct('exp_id', exp_id, 'noise_pct', n, 'is_regular', true, 'is_auto', true, 'n_replicates', n_rep);
        end
        % 16-20 / 36-40 / 56-60: irregular, auto
        for n = noise_levels
            exp_id = exp_id + 1;
            configs{end+1} = struct('exp_id', exp_id, 'noise_pct', n, 'is_regular', false, 'is_auto', true, 'n_replicates', n_rep);
        end
    end
end

function v = get_opt(s, name, default)
    if isfield(s, name), v = s.(name); else, v = default; end
end

function v = iif(cond, a, b)
    if cond, v = a; else, v = b; end
end
