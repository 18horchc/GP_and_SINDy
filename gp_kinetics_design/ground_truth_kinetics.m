function [t_gt, x_gt, y_gt, z_gt] = ground_truth_kinetics(n_points, params)
%GROUND_TRUTH_KINETICS High-resolution Robertson reaction kinetics curves for GP experiments.
%
%   [t_gt, x_gt, y_gt, z_gt] = ground_truth_kinetics()
%       Returns 500 points over [0, 50] with default Robertson parameters.
%
%   [t_gt, x_gt, y_gt, z_gt] = ground_truth_kinetics(n_points)
%       Returns n_points (default 500) over [0, 50].
%
%   [t_gt, x_gt, y_gt, z_gt] = ground_truth_kinetics(n_points, params)
%       Uses custom params (same as reaction_kinetics: x0, y0, z0, tspan).
%
%   Used as the "true" curves for RMSE and other metrics. Domain t = [0, 50].
%
%   See also: reaction_kinetics

    if nargin < 1 || isempty(n_points)
        n_points = 500;
    end
    if nargin < 2
        params = struct();
    end

    % Get curves from existing reaction_kinetics (no plot)
    [t_raw, C] = reaction_kinetics(params, false);
    x_raw = C(:, 1);
    y_raw = C(:, 2);
    z_raw = C(:, 3);

    % Interpolate to exactly n_points on [0, 50]
    t_gt = linspace(0, 50, n_points)';
    x_gt = interp1(t_raw, x_raw, t_gt, 'linear', 'extrap');
    y_gt = interp1(t_raw, y_raw, t_gt, 'linear', 'extrap');
    z_gt = interp1(t_raw, z_raw, t_gt, 'linear', 'extrap');
end
