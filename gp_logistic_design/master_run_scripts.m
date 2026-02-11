%% MASTER_RUN_SCRIPTS  Run all GP logistic experiments and aggregate metrics.
%
%   run('gp_logistic_design/master_run_scripts.m')
%
%   Runs all 20 experiments via run_all_logistic_experiments, then aggregate_logistic_metrics.
%   No state file needed—runs in one session.
%
%   For options (e.g. no plots, subset of experiments):
%     run_all_logistic_experiments(struct('make_plots', false, 'run_subset', 1:5));

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..'));

run_all_logistic_experiments();
