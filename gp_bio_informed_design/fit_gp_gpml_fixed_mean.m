function [ymu, ys2, hyp] = fit_gp_gpml_fixed_mean(x, y, xs, mean_val, gpml_path)
%FIT_GP_GPML_FIXED_MEAN  Fit GP with fixed prior mean using GPML.
%
%   [ymu, ys2, hyp] = fit_gp_gpml_fixed_mean(x, y, xs, mean_val)
%   [ymu, ys2, hyp] = fit_gp_gpml_fixed_mean(x, y, xs, mean_val, gpml_path)
%
%   x, y:      training inputs (n x d) and outputs (n x 1)
%   xs:        test inputs for prediction (m x d)
%   mean_val:  fixed prior mean m(x) = mean_val (e.g. 5)
%   gpml_path: path to GPML toolbox (optional; try common locations if empty)
%
%   Returns: ymu (predictive mean at xs), ys2 (predictive variance), hyp (full struct).
%
%   Requires: GPML toolbox (http://www.gaussianprocess.org/gpml/code/matlab/)
%   Download from: https://gitlab.com/hnickisch/gpml-matlab

    if nargin < 5, gpml_path = []; end
    ensure_gpml_path(gpml_path);

    %% Ensure column vectors
    x = x(:);
    y = y(:);
    xs = xs(:);

    %% GPML setup: meanConst (fixed), covSEiso, likGauss
    meanfunc = @meanConst;
    covfunc = @covSEiso;
    likfunc = @likGauss;

    % hyp.mean is FIXED at mean_val; we optimize only cov and lik
    hyp_fixed.mean = mean_val;
    hyp_fixed.cov = [0; 0];   % log(ell), log(sf) - initial
    hyp_fixed.lik = log(0.1);  % log(sn)

    %% Optimize only cov and lik (keep mean fixed)
    theta0 = [hyp_fixed.cov; hyp_fixed.lik];
    opt = optimoptions('fminunc', 'Display', 'off', 'MaxFunctionEvaluations', 500);

    theta_opt = fminunc(@(theta) gp_nlz(theta, mean_val, x, y, meanfunc, covfunc, likfunc), ...
        theta0, opt);

    hyp.mean = mean_val;
    hyp.cov = theta_opt(1:2);
    hyp.lik = theta_opt(3);

    %% Predict
    [ymu, ys2] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, x, y, xs);
    ys2 = max(ys2, 1e-10);
end

function nlZ = gp_nlz(theta, mean_val, x, y, meanfunc, covfunc, likfunc)
    hyp.mean = mean_val;
    hyp.cov = theta(1:2);
    hyp.lik = theta(3);
    [nlZ, ~] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, x, y);
end
