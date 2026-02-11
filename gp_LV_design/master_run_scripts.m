%% MASTER_RUN_SCRIPTS  Run all GP Lotka-Volterra experiments and aggregate metrics.
%
%   run('gp_LV_design/master_run_scripts.m')
%
%   Runs all 20 experiments via run_all_LV_experiments, then aggregate_LV_metrics.
%   No state file needed—runs in one session.
%
%   For options (e.g. no plots, subset of experiments):
%     run_all_LV_experiments(struct('make_plots', false, 'run_subset', 1:10));

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..'));

run_all_LV_experiments();
