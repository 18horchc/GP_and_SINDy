function [t_obs, y_obs] = generate_bio_data(t_gt, y_gt, data_cfg, rng_seed)
%GENERATE_BIO_DATA  Generate (t, y) observations from ground truth for bio-informed GP experiments.
%
%   [t_obs, y_obs] = generate_bio_data(t_gt, y_gt, data_cfg)
%   [t_obs, y_obs] = generate_bio_data(t_gt, y_gt, data_cfg, rng_seed)
%
%   data_cfg: struct with fields
%     .n_points     - number of time points (default 10)
%     .n_replicates - observations per time point (default 2)
%     .noise_pct    - proportional noise, sigma = noise_pct * |y_true| (default 0.05)
%     .irregular    - true for irregular sampling, false for regular (default true)
%     .t_domain     - [t_min, t_max] for sampling (default [0, 20])
%     .t_points     - if provided, use these exact times (overrides n_points, irregular)
%
%   Returns column vectors t_obs, y_obs.
%
%   See also: ground_truth_bio

    if nargin < 4 || isempty(rng_seed)
        rng_seed = 42;
    end
    rng(rng_seed);

    n_pts   = get_opt(data_cfg, 'n_points', 10);
    n_rep   = get_opt(data_cfg, 'n_replicates', 2);
    noise   = get_opt(data_cfg, 'noise_pct', 0.05);
    irregular = get_opt(data_cfg, 'irregular', true);
    t_dom   = get_opt(data_cfg, 't_domain', [0, 20]);
    t_pts   = get_opt(data_cfg, 't_points', []);

    t_min = t_dom(1);
    t_max = t_dom(2);

    if ~isempty(t_pts)
        t_points = t_pts(:);
    elseif irregular
        t_points = sort(t_min + (t_max - t_min) * rand(n_pts, 1));
    else
        t_points = linspace(t_min, t_max, n_pts)';
    end

    y_true_at_t = interp1(t_gt, y_gt, t_points, 'linear', 'extrap');
    t_obs = repelem(t_points, n_rep);
    y_obs = repelem(y_true_at_t, n_rep);

    if noise > 0
        sigma_floor = 1e-10 * (std(y_gt) + 1e-10);
        for i = 1:length(t_points)
            y_true_i = y_true_at_t(i);
            sigma_i = max(noise * abs(y_true_i), sigma_floor);
            idx = (i - 1) * n_rep + (1:n_rep);
            y_obs(idx) = y_true_i + sigma_i * randn(n_rep, 1);
        end
    end

    t_obs = t_obs(:);
    y_obs = y_obs(:);
end

function v = get_opt(cfg, name, default)
    if isempty(cfg) || ~isfield(cfg, name)
        v = default;
    else
        v = cfg.(name);
    end
end
