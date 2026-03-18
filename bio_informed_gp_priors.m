%% BIO_INFORMED_GP_PRIORS - Biologically Informed Gaussian Process Prior Categories
%
% Organizes biologically informed GP priors into three categories:
%   1. SHAPE PRIORS     - Constrain the allowable form of the latent trajectory
%   2. STRUCTURAL PRIORS - Encode broader temporal organization (event timing, periodicity)
%   3. OBSERVATION-MODEL PRIORS - Reflect statistical properties of the measurement process
%
% Each prior class includes:
%   - A clean biological example
%   - Qualitative knowledge to encode
%   - A toy equation for the latent truth
%   - GP encoding method
%
% Author: Cordelia Horch
% Date: 2026

%% ========================================================================
%  PRIOR CLASS 1: SHAPE PRIORS
%  ========================================================================
%
% (A) ACUTE TRANSIENT INFLAMMATORY RESPONSE
%     Example: Cytokine concentration or activated microglia after injury
%     Qualitative: Response starts near baseline, rises after insult, peaks, declines
%     Toy eq:  f(t) = b + A*t*exp(-c*t),  with b, A, c > 0
%     GP encoding: Informative mean (baseline b) + smooth kernel (Matérn/SqExp).
%                  Optional: pseudo-observations at late times near baseline,
%                  derivative-sign constraints for rise-then-fall.
%
% (B) BOUNDED ACTIVATION FRACTION
%     Example: Proportion of activated cells, lesion area fraction
%     Qualitative: Output in [0,1]; often starts near 0, may approach 1 or peak below
%     Toy eq:  f(t) = 1 / (1 + exp(-k*(t - t_0)))  (logistic)
%     GP encoding: Latent GP g(t) with logistic link f(t)=sigma(g(t)).
%                  Alternative: Beta likelihood for noisy proportions.
%                  Encode baseline/saturation timing via latent mean.
%
% ========================================================================
%  PRIOR CLASS 2: STRUCTURAL PRIORS
%  ========================================================================
%
% (C) DELAYED ACTIVATION
%     Example: Gene induction after drug dosing, immune activation after stroke
%     Qualitative: Known event time t_last; near baseline before, response turns on after
%     Toy eq:  f(t) = b                     for t < t_last
%              f(t) = b + A*(1 - exp(-k*(t - t_last)))  for t >= t_last
%     GP encoding: Changepoint mean or changepoint kernel (separate pre/post structures).
%                  Pseudo-observations before t_last to enforce flat baseline;
%                  different length scales before and after.
%
% (D) CIRCADIAN / RHYTHMIC BIOMARKER
%     Example: Clock gene expression, hormone cycling
%     Qualitative: Approximate period known; trajectories repeat smoothly, mild modulation
%     Toy eq:  f(t) = b + A*sin((2*pi/P)*t + phi)
%     GP encoding: Periodic kernel if regular; quasi-periodic (periodic × SqExp/Matérn)
%                  if amplitude or baseline drifts. Encode baseline in mean; fix or
%                  weakly regularize period P.
%
% ========================================================================
%  PRIOR CLASS 3: OBSERVATION-MODEL PRIORS
%  ========================================================================
%
% (E) TIME-VARYING ASSAY NOISE
%     Example: Greater variability at late inflammatory times, lower precision at high signal
%     Qualitative: Noise variance not constant; some regions more reliable than others
%     Toy eq:  y_j = f(t) + epsilon,  epsilon ~ N(0, sigma^2(t))
%     GP encoding: Heteroscedastic GP likelihood; model log(sigma^2(t)) with second GP.
%                  Or supply input-dependent observation variances if assay precision known.
%
% (F) COUNT-VALUED CELL MEASUREMENTS
%     Example: Activated-cell counts per FOV, transcript counts
%     Qualitative: Non-negative integers; variance often scales with mean (not Gaussian)
%     Toy eq:  y_j ~ Poisson(Lambda(t)),  Lambda(t) = exp(f(t))
%     GP encoding: GP prior on latent log-rate g(t). Use negative binomial for overdispersion.
%                  Enforces positivity and discrete observation process.
%
% ========================================================================

%% Demo: generate toy ground truth for each prior type
%   [t, f, params] = bio_informed_gp_priors('acute_transient');
%   [t, f, params] = bio_informed_gp_priors('bounded_activation');
%   [t, f, params] = bio_informed_gp_priors('delayed_activation');
%   [t, f, params] = bio_informed_gp_priors('circadian');
%   y = bio_informed_gp_priors('heteroscedastic_obs', t, f);
%   y = bio_informed_gp_priors('count_obs', t, f);

%% Demo: plotting (A-D)
% bio_informed_gp_priors('plot')
% bio_informed_gp_priors('plot', 0, 30)   % t from 0 to 30

%% ========================================================================
%  PRIMARY ENTRY POINT
%  ========================================================================

function varargout = bio_informed_gp_priors(mode, varargin)
%BIO_INFORMED_GP_PRIORS  Generate ground truth for biologically informed prior categories.
%
%   [t, f, params] = bio_informed_gp_priors('acute_transient', [params])
%   [t, f, params] = bio_informed_gp_priors('bounded_activation', [params])
%   [t, f, params] = bio_informed_gp_priors('delayed_activation', [params])
%   [t, f, params] = bio_informed_gp_priors('circadian', [params])
%   y = bio_informed_gp_priors('heteroscedastic_obs', t, f, [sigma_fun])
%   y = bio_informed_gp_priors('count_obs', t, f, [overdispersed])
%
%   bio_informed_gp_priors('plot', [t_min], [t_max])  - Plot all four trajectories on t in [t_min, t_max] (default [0, 20])
%
%   For shape/structural priors, returns (t, f, params). For observation priors, pass (t, f).
    switch lower(mode)
        case 'acute_transient'
            p = [];
            if nargin >= 2 && isstruct(varargin{1}), p = varargin{1}; end
            [t, f, params] = local_acute_transient([], p);
            varargout = {t, f, params};
        case 'bounded_activation'
            p = [];
            if nargin >= 2 && isstruct(varargin{1}), p = varargin{1}; end
            [t, f, params] = local_bounded_activation([], p);
            varargout = {t, f, params};
        case 'delayed_activation'
            p = [];
            if nargin >= 2 && isstruct(varargin{1}), p = varargin{1}; end
            [t, f, params] = local_delayed_activation([], p);
            varargout = {t, f, params};
        case 'circadian'
            p = [];
            if nargin >= 2 && isstruct(varargin{1}), p = varargin{1}; end
            [t, f, params] = local_circadian([], p);
            varargout = {t, f, params};
        case 'heteroscedastic_obs'
            t = varargin{1};
            f = varargin{2};
            sigma_fun = [];
            if nargin >= 4, sigma_fun = varargin{3}; end
            if isempty(sigma_fun)
                y = local_heteroscedastic_obs(t, f);
            else
                y = local_heteroscedastic_obs(t, f, sigma_fun);
            end
            varargout = {y};
        case 'count_obs'
            t = varargin{1};
            f = varargin{2};
            od = false;
            if nargin >= 4, od = varargin{3}; end
            y = local_count_obs(t, f, od);
            varargout = {y};
        case 'plot'
            t_min = 0; t_max = 20;
            if nargin >= 2 && isnumeric(varargin{1}), t_min = varargin{1}; end
            if nargin >= 3 && isnumeric(varargin{2}), t_max = varargin{2}; end
            local_plot_all(t_min, t_max);
            varargout = {};
        otherwise
            error('bio_informed_gp_priors:unknown', ...
                'Unknown mode. Use: acute_transient, bounded_activation, delayed_activation, circadian, heteroscedastic_obs, count_obs, plot.');
    end
end

%% ========================================================================
%  LOCAL FUNCTIONS (Toy equations for latent truth)
%  ========================================================================

function [t, f, params] = local_acute_transient(t, params)
% ACUTE_TRANSIENT  f(t) = b + A*t*exp(-c*t), rise-then-fall shape.
%   params.b  - baseline (default 5)
%   params.A  - amplitude scale (default 20)
%   params.c  - decay rate (default 0.6)
    if nargin < 1 || isempty(t)
        t = linspace(0, 20, 500)';
    end
    if nargin < 2, params = struct(); end
    b = get_field(params, 'b', 5);
    A = get_field(params, 'A', 20);
    c = get_field(params, 'c', 0.6);
    f = b + A * t(:) .* exp(-c * t(:));
end

function [t, f, params] = local_bounded_activation(t, params)
% BOUNDED_ACTIVATION  f(t) = 1/(1+exp(-k*(t-t_0))), logistic S-curve in [0,1].
%   params.k   - steepness (default 0.5)
%   params.t_0 - inflection time (default 10)
    if nargin < 1 || isempty(t)
        t = linspace(0, 25, 500)';
    end
    if nargin < 2, params = struct(); end
    k = get_field(params, 'k', 0.5);
    t0 = get_field(params, 't_0', 10);
    f = 1 ./ (1 + exp(-k * (t(:) - t0)));
end

function [t, f, params] = local_delayed_activation(t, params)
% DELAYED_ACTIVATION  Piecewise: flat baseline before t_last, exponential rise after.
%   params.b      - baseline (default 0.2)
%   params.A      - amplitude (default 1.5)
%   params.k      - rise rate (default 0.4)
%   params.t_last - event time (default 5)
    if nargin < 1 || isempty(t)
        t = linspace(0, 20, 500)';
    end
    if nargin < 2, params = struct(); end
    b = get_field(params, 'b', 0.2);
    A = get_field(params, 'A', 1.5);
    k = get_field(params, 'k', 0.4);
    t_last = get_field(params, 't_last', 5);
    tt = t(:);
    f = b * ones(size(tt));
    post = tt >= t_last;
    f(post) = b + A * (1 - exp(-k * (tt(post) - t_last)));
end

function [t, f, params] = local_circadian(t, params)
% CIRCADIAN  f(t) = b + A*sin((2*pi/P)*t + phi), rhythmic biomarker.
%   params.b   - baseline (default 1)
%   params.A   - amplitude (default 0.5)
%   params.P   - period (default 12)
%   params.phi - phase (default 0)
    if nargin < 1 || isempty(t)
        t = linspace(0, 72, 500)';  % 3 cycles for ~24h period
    end
    if nargin < 2, params = struct(); end
    b = get_field(params, 'b', 1);
    A = get_field(params, 'A', 0.5);
    P = get_field(params, 'P', 12);
    phi = get_field(params, 'phi', 0);
    f = b + A * sin((2*pi/P) * t(:) + phi);
end

function y = local_heteroscedastic_obs(t, f, sigma_fun)
% HETEROSCEDASTIC_OBS  y_j = f(t) + epsilon, epsilon ~ N(0, sigma^2(t)).
%   sigma_fun(t) returns std at each t. Default: sigma increases at late times.
%   Example: sigma(t) = 0.05 + 0.02*t  (more noise late)
    if nargin < 3
        sigma_fun = @(tt) 0.05 + 0.02*tt;  % time-varying noise
    end
    sigma = sigma_fun(t(:));
    y = f(:) + sigma .* randn(size(f));
end

function y = local_count_obs(t, f, overdispersed)
% COUNT_OBS  y_j ~ Poisson(exp(f(t))) or Negative Binomial if overdispersed.
%   f(t) is log-rate; ensures positive integer counts.
    if nargin < 3, overdispersed = false; end
    lambda = exp(f(:));
    lambda = max(lambda, 1e-10);  % numerical stability
    if overdispersed
        % Negative binomial: var = mu + mu^2/r; use r = 5 for moderate overdispersion
        r = 5;
        p = r ./ (r + lambda);
        y = nbinrnd(r, p);
    else
        y = poissrnd(lambda);
    end
end

function local_plot_all(t_min, t_max)
% Plot all four trajectory types on t in [t_min, t_max]. Default [0, 20].
    if nargin < 1, t_min = 0; end
    if nargin < 2, t_max = 20; end
    t = linspace(t_min, t_max, 500)';

    [~, f1] = local_acute_transient(t, struct());
    [~, f2] = local_bounded_activation(t, struct());
    [~, f3] = local_delayed_activation(t, struct());
    [~, f4] = local_circadian(t, struct());

    figure('Name', 'Biologically Informed GP Priors - Toy Trajectories');
    subplot(2, 2, 1);
    plot(t, f1, 'b-', 'LineWidth', 1.5);
    xlabel('t'); ylabel('f(t)'); title('(A) Acute Transient');
    xlim([t_min t_max]); grid on;

    subplot(2, 2, 2);
    plot(t, f2, 'r-', 'LineWidth', 1.5);
    xlabel('t'); ylabel('f(t)'); title('(B) Bounded Activation');
    xlim([t_min t_max]); ylim([0 1]); grid on;

    subplot(2, 2, 3);
    plot(t, f3, 'g-', 'LineWidth', 1.5);
    xlabel('t'); ylabel('f(t)'); title('(C) Delayed Activation');
    xlim([t_min t_max]); grid on;

    subplot(2, 2, 4);
    plot(t, f4, 'm-', 'LineWidth', 1.5);
    xlabel('t'); ylabel('f(t)'); title('(D) Circadian / Rhythmic');
    xlim([t_min t_max]); grid on;

    sgtitle(sprintf('Biologically Informed GP Priors - Latent Trajectories (t \\in [%g, %g])', t_min, t_max));
end

function v = get_field(s, fld, default)
    if isempty(s) || ~isfield(s, fld)
        v = default;
    else
        v = s.(fld);
    end
end
