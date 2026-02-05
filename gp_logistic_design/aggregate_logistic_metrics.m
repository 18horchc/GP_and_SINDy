%% Aggregate all logistic growth experiment metrics into one table/file
%
% Loads saved results from log_exp_01 through log_exp_05 (0% to 20% noise),
% merges point-estimate, probabilistic, and calibration tables, adds a
% Noise column, and saves one combined CSV (and .mat) for comparison.
%
% Prerequisite: Run log_exp_01_regular_0noise.m through log_exp_05_regular_20noise.m
%               with the Save (optional) section uncommented so .mat files exist.
%
% Output: logistic_metrics_all.csv and logistic_metrics_all.mat in this folder.
%         Columns: Noise, N, Kernel, RMSE, MAE, R2, NLPD, MSLL, CRPS, sMSE, Coverage, NLML.

clear; clc;

script_dir = fileparts(mfilename('fullpath'));

% Experiment file stems and noise level (fraction)
exp_files = {
    'results_exp_01_regular_0noise.mat',  0;
    'results_exp_02_regular_1noise.mat', 0.01;
    'results_exp_03_regular_5noise.mat', 0.05;
    'results_exp_04_regular_10noise.mat', 0.10;
    'results_exp_05_regular_20noise.mat', 0.20
};

T_all = [];
for ix = 1:size(exp_files, 1)
    fname = exp_files{ix, 1};
    noise = exp_files{ix, 2};
    fpath = fullfile(script_dir, fname);
    if ~isfile(fpath)
        warning('Missing %s — run %s first with Save uncommented.', fname, strrep(fname, 'results_', 'log_'));
        continue
    end
    S = load(fpath, 'T_point', 'T_prob', 'T_calib');
    T_point = S.T_point;
    T_prob  = S.T_prob;
    T_calib = S.T_calib;
    T_point.Noise = repmat(noise, height(T_point), 1);
    T_merged = join(T_point, T_prob, 'Keys', {'Kernel', 'N'});
    T_merged = join(T_merged, T_calib, 'Keys', {'Kernel', 'N'});
    % Order columns: Noise, N, Kernel, then point, prob, calib metrics
    T_merged = T_merged(:, {'Noise', 'N', 'Kernel', 'RMSE', 'MAE', 'R2', 'NLPD', 'MSLL', 'CRPS', 'sMSE', 'Coverage', 'NLML'});
    T_all = [T_all; T_merged];
end

if isempty(T_all)
    error('No experiment files found. Run log_exp_01 through log_exp_05 with Save uncommented.');
end

out_csv = fullfile(script_dir, 'logistic_metrics_all.csv');
out_mat = fullfile(script_dir, 'logistic_metrics_all.mat');
writetable(T_all, out_csv);
save(out_mat, 'T_all');
fprintf('Saved %d rows to %s and %s\n', height(T_all), out_csv, out_mat);
