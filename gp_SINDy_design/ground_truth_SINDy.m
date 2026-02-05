function [t_gt, M1_gt, M2_gt] = ground_truth_SINDy(n_points, params)
%GROUND_TRUTH_SINDY High-resolution SINDy+MCMC microglia curves for GP experiments.
%
%   [t_gt, M1_gt, M2_gt] = ground_truth_SINDy()
%       Returns 500 points over [0, 50] with default SINDy parameters.
%
%   [t_gt, M1_gt, M2_gt] = ground_truth_SINDy(n_points)
%       Returns n_points (default 500) over [0, 50].
%
%   [t_gt, M1_gt, M2_gt] = ground_truth_SINDy(n_points, params)
%       Uses custom params (same as sindy_mcmc_microglia: theta, M0, tspan).
%
%   Used as the "true" curves for RMSE and other metrics. Domain t = [0, 50].
%
%   See also: sindy_mcmc_microglia

    if nargin < 1 || isempty(n_points)
        n_points = 500;
    end
    if nargin < 2
        params = struct();
    end

    % Get curves from existing sindy_mcmc_microglia (no plot)
    [t_raw, M1_raw, M2_raw, ~] = sindy_mcmc_microglia(params, false);

    % Interpolate to exactly n_points on [0, 50]
    t_gt = linspace(0, 50, n_points)';
    M1_gt = interp1(t_raw, M1_raw, t_gt, 'linear', 'extrap');
    M2_gt = interp1(t_raw, M2_raw, t_gt, 'linear', 'extrap');
end
