function [d_mean, d_var] = predict_deriv_monotonic_ep(T, Y, Tm, hyp, mu_site, sigma2_site, ts)
%PREDICT_DERIV_MONOTONIC_EP  Posterior mean/variance of df/dx at points ts.
%
%   [d_mean, d_var] = predict_deriv_monotonic_ep(T, Y, Tm, hyp, mu_site, sigma2_site, ts)
%
%   Use after fit_gp_monotonic_ep; pass mu_site, sigma2_site from its output.
%
%   See also: fit_gp_monotonic_ep

    T = T(:); Y = Y(:); Tm = Tm(:); ts = ts(:);
    N = length(T); M = length(Tm);
    eta = exp(hyp(1)); ell = exp(hyp(2)); sigma = exp(hyp(3));
    [K_ff, K_fd, K_dd] = build_cov_blocks(T, Tm, eta, ell);
    K_joint = [K_ff + sigma^2*eye(N), K_fd; K_fd', K_dd];
    tilde_Sigma = diag(sigma2_site);
    tilde_Sigma_joint = blkdiag(sigma^2*eye(N), tilde_Sigma);
    tilde_mu_joint = [Y; mu_site];
    W = K_joint + tilde_Sigma_joint;
    alpha = W \ tilde_mu_joint;
    nts = length(ts);
    d_mean = zeros(nts, 1);
    d_var = zeros(nts, 1);
    for i = 1:nts
        kd = build_k_star_deriv(ts(i), T, Tm, eta, ell);
        d_mean(i) = kd' * alpha;
        kdd_ss = eta^2 / ell^2;
        d_var(i) = max(kdd_ss - kd' * (W \ kd), 1e-10);
    end
end

function [K_ff, K_fd, K_dd] = build_cov_blocks(T, Tm, eta, ell)
    N = length(T); M = length(Tm);
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
            K_fd(i,j) = k_val * r / ell^2;
        end
    end
    K_dd = zeros(M, M);
    for i = 1:M
        for j = 1:M
            r = Tm(i) - Tm(j);
            k_val = eta^2 * exp(-0.5*r^2/ell^2);
            K_dd(i,j) = k_val * (1/ell^2 - r^2/ell^4);
        end
    end
end

function kd = build_k_star_deriv(xs, T, Tm, eta, ell)
    N = length(T); M = length(Tm);
    kd = zeros(N + M, 1);
    for i = 1:N
        r = xs - T(i);
        k_val = eta^2 * exp(-0.5*r^2/ell^2);
        kd(i) = k_val * (-r) / ell^2;
    end
    for i = 1:M
        r = xs - Tm(i);
        k_val = eta^2 * exp(-0.5*r^2/ell^2);
        kd(N+i) = k_val * (1/ell^2 - r^2/ell^4);
    end
end
