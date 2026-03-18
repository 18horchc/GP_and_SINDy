function [t_gt, y_gt, params] = ground_truth_bio(toy_problem, n_points, params)
%GROUND_TRUTH_BIO  High-resolution ground truth for bio-informed GP experiments.
%
%   [t_gt, y_gt, params] = ground_truth_bio('acute_transient')
%   [t_gt, y_gt, params] = ground_truth_bio('acute_transient', n_points)
%   [t_gt, y_gt, params] = ground_truth_bio('acute_transient', n_points, params)
%
%   Returns (t_gt, y_gt) over t in [0, 20] for acute_transient. For other toy
%   problems, the domain may differ. params returns the toy equation parameters.
%
%   See also: bio_informed_gp_priors

    if nargin < 2 || isempty(n_points)
        n_points = 500;
    end
    if nargin < 3
        params = struct();
    end

    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root, '..'));

    switch lower(toy_problem)
        case 'acute_transient'
            [t_raw, y_raw, params] = bio_informed_gp_priors('acute_transient', params);
            t_dom = get_opt(params, 't_domain', [0, 20]);
            t_gt = linspace(t_dom(1), t_dom(2), n_points)';
            y_gt = interp1(t_raw, y_raw, t_gt, 'linear', 'extrap');
        otherwise
            error('ground_truth_bio:unknown', 'Unknown toy problem: %s', toy_problem);
    end
end

function v = get_opt(s, fld, default)
    if isempty(s) || ~isfield(s, fld)
        v = default;
    else
        v = s.(fld);
    end
end
