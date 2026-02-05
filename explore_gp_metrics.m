%% EXPLORE_GP_METRICS  Load, filter, and plot aggregated GP experiment metrics
%
% Prerequisite: Run the aggregate_*_metrics.m script in each design folder so that
%   logistic_metrics_all.mat, LV_metrics_all.mat, kinetics_metrics_all.mat,
%   SINDy_metrics_all.mat (and/or .csv) exist.
%
% Usage: Set the CONFIG block below, then run the script. One figure is produced.
%   To get a different plot, change CONFIG and run again (or use the EXAMPLES section).
%
% Columns in aggregated tables: Noise, N, Kernel, [State], Regular, Optimization,
%   N_per_Time, RMSE, MAE, R2, NLPD, MSLL, CRPS, sMSE, Coverage, NLML.

clear; clc; close all;

root_dir = fileparts(mfilename('fullpath'));

%% ========== CONFIG — change these, then run the script ==========
MODEL            = 'LV';     % 'logistic' | 'LV' | 'kinetics' | 'SINDy'
METRIC           = 'R2';         % 'RMSE' | 'MAE' | 'R2' | 'NLPD' | 'MSLL' | 'CRPS' | 'sMSE' | 'Coverage' | 'NLML'
X_VAR            = 'Noise';        % e.g. 'Noise' or 'N'
LINE_VAR         = 'Kernel';       % 'Kernel' or 'State' (for multi-state models)
FILTER_N         = 25;             % [] = all N; or 50, or [5 10 25 50]
FILTER_REGULAR   = 'Yes';          % 'Yes' | 'No' | [] (all)
FILTER_OPT       = 'Default';      % 'Default' | 'Auto' | [] (all)
FILTER_NPT       = 1;             % 1 | 3 | 5 | [] (all)
FILTER_STATE     = [];             % [] = all states; or {'M1','M2'}, {'Prey'}, {'X','Z'}, etc.
FILTER_NOISE     = [];             % [] = all noise; or 0, 0.01, 0.05, 0.10, 0.20, or vector
SUBPLOT_BY_STATE = true;          % true = one subplot per state (LV/kinetics/SINDy); false = one plot
%% ================================================================

% Build filter spec from CONFIG
filter_spec = struct();
filter_spec.N            = FILTER_N;
filter_spec.Regular      = FILTER_REGULAR;
filter_spec.Optimization = FILTER_OPT;
filter_spec.N_per_Time   = FILTER_NPT;
filter_spec.State        = FILTER_STATE;
filter_spec.Noise        = FILTER_NOISE;

% Load, filter, plot
T = load_metrics(root_dir, MODEL);
if isempty(T)
    error('No table loaded. Run aggregate_%s_metrics.m in gp_%s_design first.', lower(MODEL), lower(MODEL));
end
T = filter_metrics(T, filter_spec, MODEL);
if isempty(T)
    warning('No rows after filtering. Relax filters and try again.');
else
    plot_metric(T, MODEL, X_VAR, METRIC, LINE_VAR, SUBPLOT_BY_STATE);
end

%% ========== HELPER FUNCTIONS ==========

function T = load_metrics(root_dir, model)
    % Load aggregated metrics table for the given model. Prefer .mat, fall back to .csv.
    [folder, fname] = model_paths(root_dir, model);
    mat_path = fullfile(folder, [fname '.mat']);
    csv_path = fullfile(folder, [fname '.csv']);
    if isfile(mat_path)
        S = load(mat_path, 'T_all');
        T = S.T_all;
    elseif isfile(csv_path)
        T = readtable(csv_path);
    else
        T = [];
    end
end

function [folder, fname] = model_paths(root_dir, model)
    switch lower(model)
        case 'logistic'
            folder = fullfile(root_dir, 'gp_logistic_design');
            fname  = 'logistic_metrics_all';
        case 'lv'
            folder = fullfile(root_dir, 'gp_LV_design');
            fname  = 'LV_metrics_all';
        case 'kinetics'
            folder = fullfile(root_dir, 'gp_kinetics_design');
            fname  = 'kinetics_metrics_all';
        case 'sindy'
            folder = fullfile(root_dir, 'gp_SINDy_design');
            fname  = 'SINDy_metrics_all';
        otherwise
            error('Unknown model: %s. Use logistic, LV, kinetics, or SINDy.', model);
    end
end

function T = filter_metrics(T, spec, model)
    % Apply filters. Empty spec fields mean "no filter".
    if isfield(spec, 'N') && ~isempty(spec.N)
        if isscalar(spec.N)
            T = T(T.N == spec.N, :);
        else
            T = T(ismember(T.N, spec.N), :);
        end
    end
    if isfield(spec, 'Noise') && ~isempty(spec.Noise)
        if isscalar(spec.Noise)
            T = T(T.Noise == spec.Noise, :);
        else
            T = T(ismember(T.Noise, spec.Noise), :);
        end
    end
    if isfield(spec, 'Regular') && ~isempty(spec.Regular)
        T = T(T.Regular == spec.Regular | T.Regular == string(spec.Regular), :);
    end
    if isfield(spec, 'Optimization') && ~isempty(spec.Optimization)
        T = T(T.Optimization == spec.Optimization | T.Optimization == string(spec.Optimization), :);
    end
    if isfield(spec, 'N_per_Time') && ~isempty(spec.N_per_Time)
        T = T(T.N_per_Time == spec.N_per_Time, :);
    end
    if isfield(spec, 'State') && ~isempty(spec.State) && ismember('State', T.Properties.VariableNames)
        T = T(ismember(string(T.State), string(spec.State)), :);
    end
end

function plot_metric(T, model, x_var, y_var, line_var, subplot_by_state)
    % Plot y_var vs x_var with one line per line_var group. If model has State and subplot_by_state, one subplot per state.
    has_state = ismember('State', T.Properties.VariableNames);
    if has_state && subplot_by_state
        states = unique(T.State, 'stable');
        nst = numel(states);
        for s = 1:nst
            Ts = T(T.State == states(s) | string(T.State) == string(states(s)), :);
            subplot(1, nst, s);
            plot_one(Ts, x_var, y_var, line_var);
            title(sprintf('%s — %s', model, string(states(s))));
        end
        sgtitle(sprintf('%s: %s vs %s (by %s)', model, y_var, x_var, line_var));
    else
        plot_one(T, x_var, y_var, line_var);
        title(sprintf('%s: %s vs %s (by %s)', model, y_var, x_var, line_var));
    end
    xlabel(x_var);
    ylabel(y_var);
    grid on;
end

function plot_one(T, x_var, y_var, line_var)
    % Single axes: y_var vs x_var, lines grouped by line_var.
    groups = unique(T.(line_var), 'stable');
    if isnumeric(groups), groups = num2cell(groups); end
    hold on;
    for g = 1:numel(groups)
        gr = groups(g);
        if iscell(gr), gr = gr{1}; end
        if isnumeric(T.(line_var))
            mask = T.(line_var) == gr;
        else
            % Convert to string so comparison works for categorical, string, or cell
            mask = strcmp(string(T.(line_var)), string(gr));
        end
        Tg = T(mask, :);
        Tg = sortrows(Tg, x_var);
        plot(Tg.(x_var), Tg.(y_var), '-o', 'LineWidth', 1.5, 'DisplayName', string(gr));
    end
    legend('Location', 'best');
    hold off;
end

%% ========== EXAMPLES: different plots ==========
%
% Copy a block below into the CONFIG section (replace the existing CONFIG values), then run.
%
% --- Example 1: RMSE vs Noise for logistic, regular sampling, N=50, one line per kernel ---
%   MODEL            = 'logistic';
%   METRIC           = 'RMSE';
%   X_VAR            = 'Noise';
%   LINE_VAR         = 'Kernel';
%   FILTER_N         = 50;
%   FILTER_REGULAR   = 'Yes';
%   FILTER_OPT       = 'Default';
%   FILTER_NPT       = 1;
%   FILTER_STATE     = [];
%   FILTER_NOISE     = [];
%   SUBPLOT_BY_STATE = true;
%   → One figure: RMSE on y, Noise on x, five lines (one per kernel).
%
% --- Example 2: RMSE vs N for logistic, 0% noise only, lines = kernel ---
%   MODEL            = 'logistic';
%   METRIC           = 'RMSE';
%   X_VAR            = 'N';
%   LINE_VAR         = 'Kernel';
%   FILTER_N         = [];
%   FILTER_REGULAR   = 'Yes';
%   FILTER_OPT       = 'Default';
%   FILTER_NPT       = 1;
%   FILTER_STATE     = [];
%   FILTER_NOISE     = 0;
%   SUBPLOT_BY_STATE = true;
%   → One figure: RMSE vs N at 0% noise, one line per kernel.
%
% --- Example 3: SINDy, Coverage vs Noise, N=50, one subplot per state (M1, M2), lines = kernel ---
%   MODEL            = 'SINDy';
%   METRIC           = 'Coverage';
%   X_VAR            = 'Noise';
%   LINE_VAR         = 'Kernel';
%   FILTER_N         = 50;
%   FILTER_REGULAR   = 'Yes';
%   FILTER_OPT       = 'Default';
%   FILTER_NPT       = 1;
%   FILTER_STATE     = [];
%   FILTER_NOISE     = [];
%   SUBPLOT_BY_STATE = true;
%   → One figure with two subplots (M1, M2); each subplot has Coverage vs Noise, one line per kernel.
%
% --- Example 4: Kinetics, NLPD vs N for 5% noise, regular, lines = kernel, one subplot per state (X,Y,Z) ---
%   MODEL            = 'kinetics';
%   METRIC           = 'NLPD';
%   X_VAR            = 'N';
%   LINE_VAR         = 'Kernel';
%   FILTER_N         = [];
%   FILTER_REGULAR   = 'Yes';
%   FILTER_OPT       = 'Default';
%   FILTER_NPT       = 1;
%   FILTER_STATE     = [];
%   FILTER_NOISE     = 0.05;
%   SUBPLOT_BY_STATE = true;
%   → One figure with three subplots (X, Y, Z); each has NLPD vs N at 5% noise, one line per kernel.
%
% --- Example 5: LV, R² vs Noise, N=25, Prey only, lines = kernel ---
%   MODEL            = 'LV';
%   METRIC           = 'R2';
%   X_VAR            = 'Noise';
%   LINE_VAR         = 'Kernel';
%   FILTER_N         = 25;
%   FILTER_REGULAR   = 'Yes';
%   FILTER_OPT       = 'Default';
%   FILTER_NPT       = 1;
%   FILTER_STATE     = {'Prey'};
%   FILTER_NOISE     = [];
%   SUBPLOT_BY_STATE = false;
%   → Single plot: R² vs Noise for Prey only, one line per kernel.
%
% --- Example 6: Compare kernels by MAE for logistic, all N, 10% noise ---
%   MODEL            = 'logistic';
%   METRIC           = 'MAE';
%   X_VAR            = 'N';
%   LINE_VAR         = 'Kernel';
%   FILTER_N         = [];
%   FILTER_REGULAR   = 'Yes';
%   FILTER_OPT       = 'Default';
%   FILTER_NPT       = 1;
%   FILTER_STATE     = [];
%   FILTER_NOISE     = 0.10;
%   SUBPLOT_BY_STATE = true;
%   → One figure: MAE vs N at 10% noise, one line per kernel.
