%% Aggregate all SINDy experiment metrics into one table/file
%
% Loads saved results from SINDy_exp 01-60 (20 base + 20 with 3 reps + 20 with 8 reps).
% Merges point-estimate, probabilistic, and calibration tables; adds Noise,
% Regular, Optimization, N_replicates.
%
% Prerequisite: Run experiments so .mat result files exist.
%
% Output: SINDy_metrics_all.csv and SINDy_metrics_all.mat in this folder.
%         Columns: Noise, N, N_replicates, Kernel, State, Regular, Optimization,
%         N_per_Time, RMSE, MAE, R2, NLPD, MSLL, CRPS, sMSE, Coverage, NLML.

clear; clc;

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..'));

exp_files = build_aggregate_file_list();

T_all = [];
for ix = 1:size(exp_files, 1)
    fname = exp_files{ix, 1};
    noise = exp_files{ix, 2};
    is_regular = exp_files{ix, 3};
    is_auto = exp_files{ix, 4};
    fpath = fullfile(script_dir, fname);
    if ~isfile(fpath)
        warning('Missing %s — run the corresponding experiment first.', fname);
        continue
    end
    S = load(fpath, 'T_point', 'T_prob', 'T_calib');
    T_point = S.T_point;
    T_prob  = S.T_prob;
    T_calib = S.T_calib;

    if ~ismember('N_replicates', T_point.Properties.VariableNames)
        T_point.N_replicates = ones(height(T_point), 1);
        T_prob.N_replicates = ones(height(T_prob), 1);
        T_calib.N_replicates = ones(height(T_calib), 1);
    end

    T_point.Noise = repmat(noise, height(T_point), 1);
    join_keys = {'Kernel', 'N', 'N_replicates', 'State'};
    T_merged = join(T_point, T_prob, 'Keys', join_keys);
    T_merged = join(T_merged, T_calib, 'Keys', join_keys);
    if is_regular
        T_merged.Regular = repmat(categorical("Yes"), height(T_merged), 1);
    else
        T_merged.Regular = repmat(categorical("No"), height(T_merged), 1);
    end
    if is_auto
        T_merged.Optimization = repmat(categorical("Auto"), height(T_merged), 1);
    else
        T_merged.Optimization = repmat(categorical("Default"), height(T_merged), 1);
    end
    t_max = 50;
    T_merged.N_per_Time = T_merged.N / t_max;
    T_merged = T_merged(:, {'Noise', 'N', 'N_replicates', 'Kernel', 'State', 'Regular', 'Optimization', 'N_per_Time', 'RMSE', 'MAE', 'R2', 'NLPD', 'MSLL', 'CRPS', 'sMSE', 'Coverage', 'NLML'});
    T_all = [T_all; T_merged];
end

if isempty(T_all)
    error('No experiment files found. Run experiments so .mat files exist.');
end

out_csv = fullfile(script_dir, 'SINDy_metrics_all.csv');
out_mat = fullfile(script_dir, 'SINDy_metrics_all.mat');
writetable(T_all, out_csv);
save(out_mat, 'T_all');
fprintf('Saved %d rows to %s and %s\n', height(T_all), out_csv, out_mat);

function exp_files = build_aggregate_file_list()
    noise_levels = [0, 0.01, 0.05, 0.10, 0.20];
    exp_files = {};
    exp_id = 0;
    for n_rep = [1, 3, 8]
        for n = noise_levels
            exp_id = exp_id + 1;
            exp_files = [exp_files; {get_fname(exp_id, n, true, false, n_rep), n, true, false}];
        end
        for n = noise_levels
            exp_id = exp_id + 1;
            exp_files = [exp_files; {get_fname(exp_id, n, false, false, n_rep), n, false, false}];
        end
        for n = noise_levels
            exp_id = exp_id + 1;
            exp_files = [exp_files; {get_fname(exp_id, n, true, true, n_rep), n, true, true}];
        end
        for n = noise_levels
            exp_id = exp_id + 1;
            exp_files = [exp_files; {get_fname(exp_id, n, false, true, n_rep), n, false, true}];
        end
    end
end

function fname = get_fname(exp_id, noise, is_regular, is_auto, n_rep)
    sam = iif(is_regular, 'regular', 'irregular');
    if noise == 0, ns = '0noise'; elseif noise == 0.01, ns = '1noise';
    elseif noise == 0.05, ns = '5noise'; elseif noise == 0.10, ns = '10noise';
    elseif noise == 0.20, ns = '20noise'; else, ns = sprintf('%.0fnoise', noise*100); end
    auto = iif(is_auto, '_auto', '');
    rep_suf = iif(n_rep > 1, sprintf('_%drep', n_rep), '');
    fname = sprintf('results_SINDy_exp_%02d_%s_%s%s%s.mat', exp_id, sam, ns, auto, rep_suf);
end

function v = iif(cond, a, b)
    if cond, v = a; else, v = b; end
end
