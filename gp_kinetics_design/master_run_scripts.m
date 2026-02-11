%% MASTER_RUN_SCRIPTS  Run GP kinetics experiments (aggregation is separate).
%
%   run('gp_kinetics_design/master_run_scripts.m')
%
%   Runs all 20 experiments via run_all_kinetics_experiments.
%   Does NOT run aggregate_kinetics_metrics — call that separately when ready.
%
%   Examples:
%     % Run all experiments
%     run_all_kinetics_experiments();
%
%     % Run subset, no plots, then aggregate
%     run_all_kinetics_experiments(struct('make_plots', false, 'run_subset', 1:5));
%     run('gp_kinetics_design/aggregate_kinetics_metrics.m');
%
%     % Run only zero-noise experiments (config filter)
%     run_all_kinetics_experiments(struct('config_filter', @(c) c.noise_pct == 0));
%
%     % Run only regular sampling
%     run_all_kinetics_experiments(struct('config_filter', @(c) c.is_regular));

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..'));

run_all_kinetics_experiments();
