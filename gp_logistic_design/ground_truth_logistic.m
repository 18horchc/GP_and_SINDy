function [t_gt, y_gt] = ground_truth_logistic(n_points, params)
%GROUND_TRUTH_LOGISTIC High-resolution logistic growth curve for GP experiments.
%
%   [t_gt, y_gt] = ground_truth_logistic()
%       Returns 500 points over [0, 50] with default logistic parameters.
%
%   [t_gt, y_gt] = ground_truth_logistic(n_points)
%       Returns n_points (default 500) over [0, 50].
%
%   [t_gt, y_gt] = ground_truth_logistic(n_points, params)
%       Uses custom params (same as logistic_growth: r, K, P0, tspan).
%
%   Used as the "true" curve for RMSE and other metrics. Domain t = [0, 50].
%
%   See also: logistic_growth

    if nargin < 1 || isempty(n_points)
        n_points = 500;
    end
    if nargin < 2
        params = struct();
    end

    % Get curve from existing logistic_growth (no plot)
    [t_raw, P_raw] = logistic_growth(params, false);

    % Interpolate to exactly n_points on [0, 50]
    t_gt = linspace(0, 50, n_points)';
    y_gt = interp1(t_raw, P_raw, t_gt, 'linear', 'extrap');
end
