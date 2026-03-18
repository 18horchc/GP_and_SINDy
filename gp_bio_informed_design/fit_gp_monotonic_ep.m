function [ymu, ys2, hyp, post_d_mean, post_d_var, mu_site, sigma2_site] = fit_gp_monotonic_ep(T, Y, I_mono, xs, opts)
%FIT_GP_MONOTONIC_EP  GP with monotonicity via EP (Riihimäki & Vehtari 2010).
%
%   Enforces negative derivative on interval I_mono = [a,b] using virtual
%   derivative points and probit likelihood p(m|d)=Phi(-d/nu) approximated by EP.
%
%   [ymu, ys2, hyp] = fit_gp_monotonic_ep(T, Y, I_mono, xs)
%   [ymu, ys2, hyp] = fit_gp_monotonic_ep(T, Y, I_mono, xs, opts)
%
%   INPUT:
%     T       - observed times (n x 1)
%     Y       - observed values (n x 1)
%     I_mono  - monotone interval [t_min, t_max], derivative negative inside
%     xs      - test points for prediction
%     opts    - struct with:
%       .Tm        - virtual derivative locations (default: evenly spaced in I_mono)
%       .nu        - monotonicity softness (default 1e-6)
%       .M         - number of virtual points if Tm not given (default 6)
%       .hyp_init  - [log(eta); log(ell); log(sigma)] (default auto)
%       .ep_max    - max EP iterations (default 100)
%       .ep_tol    - convergence tolerance on site params (default 1e-4)
%
%   OUTPUT:
%     ymu       - predictive mean at xs
%     ys2       - predictive variance at xs
%     hyp       - fitted [log(eta); log(ell); log(sigma)]
%     post_d_mean, post_d_var - posterior mean/variance of derivatives at Tm
%     mu_site, sigma2_site - EP site params (for predict_deriv_monotonic_ep)
%
%   See: Riihimäki & Vehtari (2010) "Gaussian processes with monotonicity
%        information", AISTATS.

    if nargin < 5, opts = struct(); end

    T = T(:);
    Y = Y(:);
    xs = xs(:);
    a = I_mono(1);
    b = I_mono(2);

    Tm      = get_opt(opts, 'Tm', []);
    nu      = get_opt(opts, 'nu', 1e-6);
    M       = get_opt(opts, 'M', 6);
    hyp_init = get_opt(opts, 'hyp_init', []);
    ep_max  = get_opt(opts, 'ep_max', 100);
    ep_tol  = get_opt(opts, 'ep_tol', 1e-4);

    if isempty(Tm)
        Tm = linspace(a, b, M)';
    end
    Tm = Tm(:);
    M = length(Tm);

    N = length(T);

    % Initial hyperparameters
    if isempty(hyp_init)
        sig_y = std(Y);
        ell_init = (max(T)-min(T)) / 4;
        hyp_init = [log(sig_y+1e-2); log(ell_init); log(0.05*sig_y)];
    end

    %% Optimize hyperparameters
    opt = optimoptions('fminunc', 'Display', 'off', 'MaxFunctionEvaluations', 300);
    hyp = fminunc(@(theta) neg_log_Z(theta, T, Y, Tm, nu, ep_max, ep_tol), hyp_init, opt);

    eta = exp(hyp(1));
    ell = exp(hyp(2));
    sigma = exp(hyp(3));

    %% Run EP to get site parameters
    [~, mu_site, sigma2_site] = run_ep(T, Y, Tm, eta, ell, sigma, nu, ep_max, ep_tol);

    %% Build joint matrix for prediction
    [K_ff, K_fd, K_dd] = build_cov_blocks(T, Tm, eta, ell);
    K_joint = [K_ff + sigma^2*eye(N), K_fd; K_fd', K_dd];
    tilde_Sigma = diag(sigma2_site);
    tilde_Sigma_joint = blkdiag(sigma^2*eye(N), tilde_Sigma);
    tilde_mu_joint = [Y; mu_site];
    W = K_joint + tilde_Sigma_joint;
    iW = inv(W);

    %% Predict at xs
    nxs = length(xs);
    ymu = zeros(nxs, 1);
    ys2 = zeros(nxs, 1);
    for i = 1:nxs
        k_star = build_k_star(xs(i), T, Tm, eta, ell);
        ymu(i) = k_star' * iW * tilde_mu_joint;
        kss = eta^2;  % k(x*,x*) for SE
        ys2(i) = max(kss - k_star' * iW * k_star, 1e-10);
    end

    if nargout >= 4
        post_d_mean = mu_site;
        post_d_var = sigma2_site;
    end
    % mu_site, sigma2_site always computed (outputs 6-7 for predict_deriv_monotonic_ep)
end

%% ========== EP and covariance helpers ==========

function n = neg_log_Z(theta, T, Y, Tm, nu, ep_max, ep_tol)
    eta = exp(theta(1));
    ell = exp(theta(2));
    sigma = exp(theta(3));
    [logZ, ~, ~] = run_ep(T, Y, Tm, eta, ell, sigma, nu, ep_max, ep_tol);
    n = -logZ;
end

function [logZ, mu_site, sigma2_site] = run_ep(T, Y, Tm, eta, ell, sigma, nu, ep_max, ep_tol)
    N = length(T);
    M = length(Tm);

    [K_ff, K_fd, K_dd] = build_cov_blocks(T, Tm, eta, ell);
    K_joint = [K_ff + sigma^2*eye(N), K_fd; K_fd', K_dd];

    % Initialize sites: uninformative (large variance) so posterior = prior initially
    sigma2_site = 1e6 * ones(M, 1);  % large = weak prior
    mu_site = zeros(M, 1);

    for iter = 1:ep_max
        sigma2_old = sigma2_site;
        mu_old = mu_site;

        tilde_Sigma = diag(sigma2_site);
        tilde_Sigma_joint = blkdiag(sigma^2*eye(N), tilde_Sigma);
        tilde_mu_joint = [Y; mu_site];

        W = K_joint + tilde_Sigma_joint;
        alpha = W \ tilde_mu_joint;
        mu_post = K_joint * alpha;  % posterior mean (GP regression identity)
        Sigma_post = K_joint - K_joint * (W \ K_joint);  % posterior covariance

        for j = 1:M
            idx = N + j;
            mu_j = mu_post(idx);
            sigma2_j = Sigma_post(idx, idx);

            % Cavity: q_{-j}(d_j) by removing site j
            prec_site = 1 / sigma2_site(j);
            prec_post = 1 / sigma2_j;
            prec_cavity = prec_post - prec_site;
            if prec_cavity <= 1e-10
                continue;
            end
            sigma2_cavity = 1 / prec_cavity;
            mu_cavity = sigma2_cavity * (mu_j / sigma2_j - mu_site(j) / sigma2_site(j));

            % Moment match tilted = cavity * Phi(-d/nu) for NEGATIVE derivative
            % For Phi(d/nu) paper adds; for Phi(-d/nu) we subtract to push d downward
            zi = -mu_cavity / (nu * sqrt(1 + sigma2_cavity / nu^2));
            zi = max(min(zi, 35), -35);  % avoid overflow in Phi

            nzi = normpdf(zi);
            Phi_zi = normcdf(zi);

            if Phi_zi < 1e-100
                Phi_zi = 1e-100;
            end

            % Site mean: subtract (not add) for Phi(-d/nu) so tilted mean shifts toward negative
            mu_site(j) = mu_cavity - sigma2_cavity * nzi / (Phi_zi * nu * sqrt(1 + sigma2_cavity/nu^2));
            sigma2_site(j) = sigma2_cavity - sigma2_cavity^2 * nzi / (Phi_zi * (nu^2 + sigma2_cavity)) * (zi + nzi/Phi_zi);

            if sigma2_site(j) <= 0
                sigma2_site(j) = 1e-10;
            end
        end

        if max(abs(sigma2_site - sigma2_old)) < ep_tol && max(abs(mu_site - mu_old)) < ep_tol
            break;
        end
    end

    % Log marginal likelihood (EP approximation)
    tilde_Sigma = diag(sigma2_site);
    tilde_Sigma_joint = blkdiag(sigma^2*eye(N), tilde_Sigma);
    tilde_mu_joint = [Y; mu_site];
    W = K_joint + tilde_Sigma_joint;

    log_det_W = sum(log(abs(eig(W))));
    if ~isfinite(log_det_W)
        log_det_W = 2 * sum(log(abs(diag(chol(W)))));
    end

    term1 = -0.5 * log_det_W - 0.5 * tilde_mu_joint' * (W \ tilde_mu_joint);

    % Cavity terms (need posterior marginals for each derivative)
    alpha_Z = W \ tilde_mu_joint;
    mu_post = K_joint * alpha_Z;
    Sigma_post = K_joint - K_joint * (W \ K_joint);
    term2 = 0;
    term3 = 0;
    for j = 1:M
        idx = N + j;
        mu_j = mu_post(idx);
        sigma2_j = Sigma_post(idx, idx);
        prec_site = 1 / sigma2_site(j);
        prec_post = 1 / sigma2_j;
        prec_cavity = prec_post - prec_site;
        if prec_cavity > 1e-10
            sigma2_cavity = 1 / prec_cavity;
            mu_cavity = sigma2_cavity * (mu_j / sigma2_j - mu_site(j) / sigma2_site(j));
            zi = -mu_cavity / (nu * sqrt(1 + sigma2_cavity/nu^2));
            zi = max(min(zi, 35), -35);
            term2 = term2 + (mu_cavity - mu_site(j))^2 / (2 * (sigma2_cavity + sigma2_site(j)));
            term3 = term3 + log(normcdf(zi)) + 0.5*log(sigma2_cavity + sigma2_site(j));
        end
    end
    logZ = term1 + term2 + term3;
end

function [K_ff, K_fd, K_dd] = build_cov_blocks(T, Tm, eta, ell)
    N = length(T);
    M = length(Tm);

    K_ff = zeros(N, N);
    for i = 1:N
        for j = 1:N
            r = T(i) - T(j);
            K_ff(i,j) = eta^2 * exp(-0.5*r^2/ell^2);
        end
    end

    K_fd = zeros(N, M);
    for i = 1:N
        for j = 1:M
            r = T(i) - Tm(j);
            k_val = eta^2 * exp(-0.5*r^2/ell^2);
            K_fd(i,j) = k_val * r / ell^2;  % Cov(f(T), df(Tm)/dx') = d/dx' k = k*(x-x')/ell^2
        end
    end

    K_dd = zeros(M, M);
    for i = 1:M
        for j = 1:M
            r = Tm(i) - Tm(j);
            k_val = eta^2 * exp(-0.5*r^2/ell^2);
            K_dd(i,j) = k_val * (1/ell^2 - r^2/ell^4);  % d2/dxdx' k
        end
    end
end

function k_star = build_k_star(xs, T, Tm, eta, ell)
    N = length(T);
    M = length(Tm);
    k_star = zeros(N + M, 1);
    for i = 1:N
        r = xs - T(i);
        k_star(i) = eta^2 * exp(-0.5*r^2/ell^2);
    end
    for i = 1:M
        r = xs - Tm(i);
        k_val = eta^2 * exp(-0.5*r^2/ell^2);
        k_star(N+i) = k_val * r / ell^2;  % Cov(f(xs), df(Tm)/dx')
    end
end

function v = get_opt(opts, name, default)
    if isempty(opts) || ~isfield(opts, name)
        v = default;
    else
        v = opts.(name);
    end
end
