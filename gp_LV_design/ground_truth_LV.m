function [t_gt, prey_gt, pred_gt] = ground_truth_LV(n_points, params)
%GROUND_TRUTH_LV High-resolution Lotka-Volterra prey and predator curves for GP experiments.
%
%   [t_gt, prey_gt, pred_gt] = ground_truth_LV()
%       Returns 500 points over [0, 50] with default LV parameters (same as lotka_volterra_model).
%
%   [t_gt, prey_gt, pred_gt] = ground_truth_LV(n_points)
%       Returns n_points (default 500) over [0, 50].
%
%   [t_gt, prey_gt, pred_gt] = ground_truth_LV(n_points, params)
%       Uses custom params (same as lotka_volterra_model: alpha, beta, delta, gamma, x0, y0, tspan).
%
%   Used as the "true" curves for RMSE and other metrics. Domain t = [0, 50].
%
%   See also: lotka_volterra_model

    if nargin < 1 || isempty(n_points)
        n_points = 500;
    end
    if nargin < 2
        params = struct();
    end

    % Get curves from existing lotka_volterra_model (no plot)
    [t_raw, prey_raw, pred_raw] = lotka_volterra_model(params, false);

    % Interpolate to exactly n_points on [0, 50]
    t_gt = linspace(0, 50, n_points)';
    prey_gt = interp1(t_raw, prey_raw, t_gt, 'linear', 'extrap');
    pred_gt = interp1(t_raw, pred_raw, t_gt, 'linear', 'extrap');
end
