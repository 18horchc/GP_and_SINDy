function [ymu, ystd, gpr] = fit_gp_bio(t_obs, y_obs, t_gt, enc_config)
%FIT_GP_BIO  Fit GP with bio-informed encoding for acute transient.
%
%   [ymu, ystd, gpr] = fit_gp_bio(t_obs, y_obs, t_gt, enc_config)
%
%   enc_config: from encoding_registry (struct with .informative_mean,
%   .kernel_matern32, .late_pseudo_obs, .derivative_constraints).
%
%   Returns predictive mean and std at t_gt, plus the fitted RegressionGP model.
%
%   Pseudo-observations and derivative guides are appended to (t_obs, y_obs)
%   before fitting. All use the same observation noise (fitrgp has one Sigma).
%
%   See also: encoding_registry, run_bio_experiment

    baseline = 0.1;  % b for acute transient
    t_max = 20;

    t_fit = t_obs(:);
    y_fit = y_obs(:);

    % --- Informative mean: pseudo-obs at boundaries to encode baseline b=0.1 ---
    % At t=0, f(0)=b exactly. At t=20, f(20)≈0.2 for default params; we use b as
    % "asymptotic baseline" (curve returns toward b as t→∞).
    if enc_config.informative_mean
        t_fit = [t_fit; 0; t_max];
        y_fit = [y_fit; baseline; baseline];
    end

    % --- Late-time pseudo-obs: (16, b), (18, b), (20, b) ---
    if enc_config.late_pseudo_obs
        t_late = [16; 18; t_max];
        y_late = baseline * ones(3, 1);
        % Avoid duplicate at t_max if informative_mean already added it
        if enc_config.informative_mean
            t_late = [16; 18];  % 20 already in from informative_mean
            y_late = baseline * ones(2, 1);
        end
        t_fit = [t_fit; t_late];
        y_fit = [y_fit; y_late];
    end

    % --- Derivative shape guide: virtual points to encourage rise-then-fall ---
    % Uses prior shape f(t)=b+A*t*exp(-c*t) with b=0.1,A=2,c=0.3 for plausible values.
    % Rising (t~1.5-4) and falling (t~12-18) phases.
    if enc_config.derivative_constraints
        A = 2; c = 0.3;
        t_rise = [1.5; 3.0];
        y_rise = baseline + A * t_rise .* exp(-c * t_rise);
        t_fall = [12; 17];
        y_fall = baseline + A * t_fall .* exp(-c * t_fall);
        t_fit = [t_fit; t_rise; t_fall];
        y_fit = [y_fit; y_rise; y_fall];
    end

    % --- Kernel choice ---
    if enc_config.kernel_matern32
        kern = 'matern32';
    else
        kern = 'squaredexponential';
    end

    % --- Fit with fitrgp ---
    base_opts = {'BasisFunction', 'constant', 'FitMethod', 'exact', ...
        'PredictMethod', 'exact', 'Standardize', true, ...
        'OptimizeHyperparameters', 'auto', ...
        'HyperparameterOptimizationOptions', struct('ShowPlots', false, 'Verbose', 0)};

    try
        gpr = fitrgp(t_fit, y_fit, 'KernelFunction', kern, base_opts{:});
        [ymu, ystd, ~] = predict(gpr, t_gt);
    catch me
        warning('fit_gp_bio:fit_failed', 'GP fit failed: %s', me.message);
        ymu = nan(size(t_gt));
        ystd = nan(size(t_gt));
        gpr = [];
    end
end
