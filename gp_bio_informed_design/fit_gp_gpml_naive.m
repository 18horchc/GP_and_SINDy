function [ymu, ys2, hyp] = fit_gp_gpml_naive(x, y, xs, gpml_path)
%FIT_GP_GPML_NAIVE  Fit naive GP with GPML: meanConst + covSEiso, all hyperparameters optimized.
%
%   [ymu, ys2, hyp] = fit_gp_gpml_naive(x, y, xs)
%   [ymu, ys2, hyp] = fit_gp_gpml_naive(x, y, xs, gpml_path)
%
%   Naive: constant mean (fitted), SE kernel. Same model class as fitrgp naive,
%   but implemented in GPML for direct comparison with GPML fixed-mean.
%
%   Returns: ymu (predictive mean), ys2 (predictive variance), hyp (fitted struct).
%
%   Requires: GPML toolbox.

    if nargin < 4, gpml_path = []; end
    ensure_gpml_path(gpml_path);

    x = x(:);
    y = y(:);
    xs = xs(:);

    meanfunc = @meanConst;
    covfunc = @covSEiso;
    likfunc = @likGauss;

    hyp = struct('mean', mean(y), 'cov', [0; 0], 'lik', log(0.1));
    hyp = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, x, y);

    [ymu, ys2] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, x, y, xs);
    ys2 = max(ys2, 1e-10);
end
