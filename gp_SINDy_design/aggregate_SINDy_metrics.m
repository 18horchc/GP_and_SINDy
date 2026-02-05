%% Aggregate all SINDy experiment metrics into one table/file
%
% Loads saved results from SINDy_exp_01 through SINDy_exp_05 (0% to 20% noise),
% merges point-estimate, probabilistic, and calibration tables, adds a
% Noise column, and saves one combined CSV (and .mat) for comparison.
%
% Prerequisite: Run SINDy_exp_01_regular_0noise.m through SINDy_exp_05_regular_20noise.m
%               with the Save (optional) section uncommented so .mat files exist.
%
% Output: SINDy_metrics_all.csv and SINDy_metrics_all.mat in this folder.
%         Columns: Noise, N, Kernel, State, Regular, Optimization, N_per_Time, RMSE, MAE, R2, NLPD, MSLL, CRPS, sMSE, Coverage, NLML.

clear; clc;

script_dir = fileparts(mfilename('fullpath'));

% Experiment file stems and noise level (fraction)
exp_files = {
    'results_SINDy_exp_01_regular_0noise.mat',  0;
    'results_SINDy_exp_02_regular_1noise.mat', 0.01;
    'results_SINDy_exp_03_regular_5noise.mat', 0.05;
    'results_SINDy_exp_04_regular_10noise.mat', 0.10;
    'results_SINDy_exp_05_regular_20noise.mat', 0.20
};

T_all = [];
for ix = 1:size(exp_files, 1)
    fname = exp_files{ix, 1};
    noise = exp_files{ix, 2};
    fpath = fullfile(script_dir, fname);
    if ~isfile(fpath)
        warning('Missing %s — run SINDy_exp_%02d_regular_*noise.m first with Save uncommented.', fname, ix);
        continue
    end
    S = load(fpath, 'T_point', 'T_prob', 'T_calib');
    T_point = S.T_point;
    T_prob  = S.T_prob;
    T_calib = S.T_calib;
    T_point.Noise = repmat(noise, height(T_point), 1);
    T_merged = join(T_point, T_prob, 'Keys', {'Kernel', 'N', 'State'});
    T_merged = join(T_merged, T_calib, 'Keys', {'Kernel', 'N', 'State'});
    T_merged.Regular = repmat(categorical("Yes"), height(T_merged), 1);
    T_merged.Optimization = repmat(categorical("Default"), height(T_merged), 1);
    T_merged.N_per_Time = ones(height(T_merged), 1);
    T_merged = T_merged(:, {'Noise', 'N', 'Kernel', 'State', 'Regular', 'Optimization', 'N_per_Time', 'RMSE', 'MAE', 'R2', 'NLPD', 'MSLL', 'CRPS', 'sMSE', 'Coverage', 'NLML'});
    T_all = [T_all; T_merged];
end

if isempty(T_all)
    error('No experiment files found. Run SINDy_exp_01 through SINDy_exp_05 with Save uncommented.');
end

out_csv = fullfile(script_dir, 'SINDy_metrics_all.csv');
out_mat = fullfile(script_dir, 'SINDy_metrics_all.mat');
writetable(T_all, out_csv);
save(out_mat, 'T_all');
fprintf('Saved %d rows to %s and %s\n', height(T_all), out_csv, out_mat);
