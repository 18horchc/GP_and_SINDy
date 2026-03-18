function [ymu, ys2, hyp] = fit_gp_deriv_obs(x_f, y_f, x_df, d_obs, xs, hyp_init, mean_val)
%FIT_GP_DERIV_OBS  Fit GP with function and derivative observations (standalone, no GPML).
%
%   [ymu, ys2, hyp] = fit_gp_deriv_obs(x_f, y_f, x_df, d_obs, xs)
%   [ymu, ys2, hyp] = fit_gp_deriv_obs(x_f, y_f, x_df, d_obs, xs, hyp_init)
%   [ymu, ys2, hyp] = fit_gp_deriv_obs(x_f, y_f, x_df, d_obs, xs, hyp_init, mean_val)
%
%   x_f, y_f:  function observations (n x 1)
%   x_df:      derivative observation locations (m x 1)
%   d_obs:     observed derivative values (m x 1)
%   xs:        test points for prediction
%   hyp_init:  [log(ell); log(sigma_f); log(sigma_n)] (optional)
%   mean_val:  constant prior mean m(x)=mean_val (default 0); dm/dx=0 for deriv obs
%
%   Encodes "slope decreases from day 2 to 10": x_df=[2;10], d_obs=[0.5;-0.5]
%
%   See also: run_naive_vs_deriv_encoding, docs/Derivative_Observations_Guide.md

    if nargin < 6, hyp_init = [0; 0; log(0.1)]; end
    if nargin < 7, mean_val = 0; end

    x_f = x_f(:); y_f = y_f(:);
    x_df = x_df(:); d_obs = d_obs(:);
    xs = xs(:);

    n_f = length(x_f);
    n_df = length(x_df);

    %% Optimize hyperparameters via NLML
    opt = optimoptions('fminunc', 'Display', 'off', 'MaxFunctionEvaluations', 500);
    hyp = fminunc(@(theta) nlz(theta, x_f, y_f, x_df, d_obs, mean_val), hyp_init, opt);

    ell = exp(hyp(1));
    sf = exp(hyp(2));
    sn = exp(hyp(3));

    %% Build full covariance K and solve
    [K, m_vec] = build_K(x_f, y_f, x_df, d_obs, ell, sf, sn, mean_val);
    z = [y_f; d_obs];
    alpha = K \ (z - m_vec);

    %% Predict at xs
    ymu = zeros(size(xs));
    ys2 = zeros(size(xs));
    for i = 1:length(xs)
        k_star = build_k_star(xs(i), x_f, x_df, ell, sf);
        ymu(i) = mean_val + k_star' * alpha;
        v_star = sf^2 - k_star' * (K \ k_star);
        ys2(i) = max(v_star, 1e-10);
    end
end

function n = nlz(theta, x_f, y_f, x_df, d_obs, mean_val)
    if nargin < 6, mean_val = 0; end
    ell = exp(theta(1));
    sf = exp(theta(2));
    sn = exp(theta(3));
    [K, m_vec] = build_K(x_f, y_f, x_df, d_obs, ell, sf, sn, mean_val);
    z = [y_f; d_obs];
    nn = size(K, 1);
    L = chol(K + 1e-8*eye(nn), 'lower');
    v = L \ (z - m_vec);
    n = 0.5*(v'*v) + sum(log(diag(L))) + 0.5*nn*log(2*pi);
end

function [K, m_vec] = build_K(x_f, y_f, x_df, d_obs, ell, sf, sn, mean_val)
    if nargin < 8, mean_val = 0; end
    n_f = length(x_f);
    n_df = length(x_df);
    n = n_f + n_df;

    K = zeros(n, n);
    m_vec = [mean_val * ones(n_f, 1); zeros(n_df, 1)];  % dm/dx = 0 for constant mean

    % K_ff
    for i = 1:n_f
        for j = 1:n_f
            r = x_f(i) - x_f(j);
            K(i,j) = sf^2 * exp(-0.5*r^2/ell^2);
        end
    end
    % K_fdf: Cov(f(x_i), df(x_j')/dx')
    for i = 1:n_f
        for j = 1:n_df
            r = x_f(i) - x_df(j);
            k_val = sf^2 * exp(-0.5*r^2/ell^2);
            K(i, n_f+j) = k_val * r / ell^2;  % d/dx' k(x,x') = k * (x-x')/ell^2
            K(n_f+j, i) = K(i, n_f+j);
        end
    end
    % K_dfdf: Cov(df(x_i)/dx, df(x_j)/dx) = ∂²k/∂x∂x'
    for i = 1:n_df
        for j = 1:n_df
            r = x_df(i) - x_df(j);
            k_val = sf^2 * exp(-0.5*r^2/ell^2);
            K(n_f+i, n_f+j) = k_val * (1/ell^2 - r^2/ell^4);
        end
    end

    K = K + sn^2 * eye(n);
end

function k_star = build_k_star(xs, x_f, x_df, ell, sf)
    n_f = length(x_f);
    n_df = length(x_df);
    k_star = zeros(n_f + n_df, 1);
    for i = 1:n_f
        r = xs - x_f(i);
        k_star(i) = sf^2 * exp(-0.5*r^2/ell^2);
    end
    for i = 1:n_df
        r = xs - x_df(i);
        k_val = sf^2 * exp(-0.5*r^2/ell^2);
        k_star(n_f+i) = k_val * r / ell^2;  % d/dx' k(xs, x_df)
    end
end
