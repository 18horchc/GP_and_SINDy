%% MASTER_RUN_SCRIPTS  Run GP SINDy experiments (aggregation is separate).
%
%   run('gp_SINDy_design/master_run_scripts.m')
%
%   Runs all 20 experiments via run_all_SINDy_experiments.
%   Does NOT run aggregate_SINDy_metrics — call that separately when ready.
%
%   Examples:
%     % Run all experiments
%     run_all_SINDy_experiments();
%
%     % Run subset, no plots, then aggregate
%     run_all_SINDy_experiments(struct('make_plots', false, 'run_subset', 1:5));
%     run('gp_SINDy_design/aggregate_SINDy_metrics.m');
%
%     % Run only zero-noise experiments (config filter)
%     run_all_SINDy_experiments(struct('config_filter', @(c) c.noise_pct == 0));
%
%     % Run only regular sampling
%     run_all_SINDy_experiments(struct('config_filter', @(c) c.is_regular));

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..'));

run_all_SINDy_experiments();
