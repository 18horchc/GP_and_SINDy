%% Aggregate all SINDy experiment metrics into one table/file
%
% Loads saved results from SINDy_exp_01 through SINDy_exp_10 (regular 01-05,
% irregular 06-10), merges point-estimate, probabilistic, and calibration
% tables, adds Noise and Regular columns, and saves one combined CSV (and .mat).
%
% Prerequisite: Run SINDy_exp_01 through SINDy_exp_10 so .mat result files exist.
%
% Output: SINDy_metrics_all.csv and SINDy_metrics_all.mat in this folder.
%         Columns: Noise, N, Kernel, State, Regular, Optimization, N_per_Time, RMSE, MAE, R2, NLPD, MSLL, CRPS, sMSE, Coverage, NLML.

clear; clc;

script_dir = fileparts(mfilename('fullpath'));

% Experiment file stems: {filename, noise_fraction, is_regular}
exp_files = {
    'results_SINDy_exp_01_regular_0noise.mat',  0,    true;
    'results_SINDy_exp_02_regular_1noise.mat', 0.01, true;
    'results_SINDy_exp_03_regular_5noise.mat', 0.05, true;
    'results_SINDy_exp_04_regular_10noise.mat', 0.10, true;
    'results_SINDy_exp_05_regular_20noise.mat', 0.20, true;
    'results_SINDy_exp_06_irregular_0noise.mat',  0,    false;
    'results_SINDy_exp_07_irregular_1noise.mat', 0.01, false;
    'results_SINDy_exp_08_irregular_5noise.mat', 0.05, false;
    'results_SINDy_exp_09_irregular_10noise.mat', 0.10, false;
    'results_SINDy_exp_10_irregular_20noise.mat', 0.20, false
};

T_all = [];
for ix = 1:size(exp_files, 1)
    fname = exp_files{ix, 1};
    noise = exp_files{ix, 2};
    is_regular = exp_files{ix, 3};
    fpath = fullfile(script_dir, fname);
    if ~isfile(fpath)
        warning('Missing %s — run the corresponding SINDy_exp_*.m first.', fname);
        continue
    end
    S = load(fpath, 'T_point', 'T_prob', 'T_calib');
    T_point = S.T_point;
    T_prob  = S.T_prob;
    T_calib = S.T_calib;
    T_point.Noise = repmat(noise, height(T_point), 1);
    T_merged = join(T_point, T_prob, 'Keys', {'Kernel', 'N', 'State'});
    T_merged = join(T_merged, T_calib, 'Keys', {'Kernel', 'N', 'State'});
    if is_regular
        T_merged.Regular = repmat(categorical("Yes"), height(T_merged), 1);
    else
        T_merged.Regular = repmat(categorical("No"), height(T_merged), 1);
    end
    T_merged.Optimization = repmat(categorical("Default"), height(T_merged), 1);
    T_merged.N_per_Time = ones(height(T_merged), 1);
    T_merged = T_merged(:, {'Noise', 'N', 'Kernel', 'State', 'Regular', 'Optimization', 'N_per_Time', 'RMSE', 'MAE', 'R2', 'NLPD', 'MSLL', 'CRPS', 'sMSE', 'Coverage', 'NLML'});
    T_all = [T_all; T_merged];
end

if isempty(T_all)
    error('No experiment files found. Run SINDy_exp_01 through SINDy_exp_10 so .mat files exist.');
end

out_csv = fullfile(script_dir, 'SINDy_metrics_all.csv');
out_mat = fullfile(script_dir, 'SINDy_metrics_all.mat');
writetable(T_all, out_csv);
save(out_mat, 'T_all');
fprintf('Saved %d rows to %s and %s\n', height(T_all), out_csv, out_mat);
