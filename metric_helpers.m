function varargout = metric_helpers(metric_name, varargin)
% METRIC_HELPERS  Compute GP metrics per PDF "GP Toy Problem Project - Metrics".
%
% Shared by gp_logistic_design and gp_LV_design. Keep in project root (GP_and_SINDy).
%
% Point-estimate:
%   [rmse, mae, r2] = metric_helpers('point_estimate', y_true, y_pred)
%
% Probabilistic (predictive mean and std at test points, ground truth, training y for baseline):
%   [nlpd, msll, crps] = metric_helpers('probabilistic', y_true, ymu, ystd, y_train)
%
% Calibration (same y_true, ymu, ystd; optional gpr model for NLML):
%   [smse, coverage, nlml] = metric_helpers('calibration', y_true, ymu, ystd, y_train, gpr)
%
% NLML only (from fitted gpr model):
%   nlml = metric_helpers('nlml', gpr)

switch lower(metric_name)
    case 'point_estimate'
        % [RMSE, MAE, R2]; args: y_true, y_pred
        y_true = varargin{1};
        y_pred = varargin{2};
        n = numel(y_true);
        res = y_true(:) - y_pred(:);
        rmse = sqrt(mean(res.^2));
        mae = mean(abs(res));
        ss_res = sum(res.^2); %SSE
        ss_tot = sum((y_true(:) - mean(y_true(:))).^2); %SST
        if ss_tot > 0
            r2 = 1 - ss_res / ss_tot;
        else
            r2 = NaN;
        end
        varargout = {rmse, mae, r2};

    case 'probabilistic'
        % [NLPD, MSLL, CRPS]; args: y_true, ymu, ystd, y_train
        y_true = varargin{1}(:);
        ymu = varargin{2}(:);
        ystd = varargin{3}(:);
        y_train = varargin{4}(:);
        n = numel(y_true);
        % Avoid zero/negative std
        ystd = max(ystd, 1e-10);

        % NLPD = -sum log p(y_i* | X,y,x_i*); log p = log N(y; mu, sigma^2)
        log_p = -0.5*log(2*pi) - log(ystd) - 0.5*((y_true - ymu)./ystd).^2;
        nlpd = -sum(log_p);
        %Might be overcomplicating above calculation? Could use:
        % -sum(normlogpdf(y_true, ymu, sqrt(ystd)));


        % MSLL = (1/n)*sum [ -log p(y_*^i|GP) + log p(y_*^i|train_data) ]
        % Baseline: N(mean(y_train), var(y_train))
        m_train = mean(y_train);
        v_train = var(y_train);
        if v_train <= 0
            v_train = 1e-10;
        end
        log_p_baseline = -0.5*log(2*pi) - 0.5*log(v_train) - 0.5*((y_true - m_train).^2 / v_train);
        msll = mean(-log_p + log_p_baseline);

        % CRPS(μ,σ,y) = σ * [ (y-μ)/σ * (2*Φ((y-μ)/σ) - 1) + 2*φ((y-μ)/σ) - 1/√π ]
        z = (y_true - ymu) ./ ystd;
        crps_per = ystd .* (z .* (2*normcdf(z) - 1) + 2*normpdf(z) - 1/sqrt(pi));
        crps = mean(crps_per);
        varargout = {nlpd, msll, crps};

    case 'calibration'
        % [sMSE, Coverage, NLML]; args: y_true, ymu, ystd, y_train, gpr (optional)
        y_true = varargin{1}(:);
        ymu = varargin{2}(:);
        ystd = varargin{3}(:);
        y_train = varargin{4}(:);
        gpr = [];
        if nargin >= 6
            gpr = varargin{5};
        end

        % Standardized MSE = MSE / Var(y_true)
        mse = mean((y_true - ymu).^2);
        vy = var(y_true);
        if vy > 0
            smse = mse / vy;
        else
            smse = NaN;
        end

        % Coverage: fraction of y_true in 95% interval
        ylo = ymu - 1.96 * ystd;
        yhi = ymu + 1.96 * ystd;
        coverage = mean((y_true >= ylo) & (y_true <= yhi));

        % NLML from model if provided
        if ~isempty(gpr)
            nlml = metric_helpers('nlml', gpr);
        else
            nlml = NaN;
        end
        varargout = {smse, coverage, nlml};

    case 'nlml'
        % Negative log marginal likelihood from fitted RegressionGP
        gpr = varargin{1};
        X = gpr.X;
        y = gpr.Y(:);
        n = numel(y);
        sigma = gpr.Sigma;
        params = gpr.KernelInformation.KernelParameters;
        kname = lower(gpr.KernelInformation.Name);

        % If Standardize was true, fitrgp fits on standardized data; use stored internals
        % RegressionGP stores original X,Y; kernel is evaluated on standardized X
        if isprop(gpr, 'Standardize') && gpr.Standardize
            X_std = (X - mean(X)) ./ (std(X) + 1e-10);
            y_std = (y - mean(y)) / (std(y) + 1e-10);
        else
            X_std = X;
            y_std = y;
        end

        K = build_gram(X_std, kname, params);
        Ky = K + (sigma^2) * eye(n);
        try
            L = chol(Ky, 'lower');
            z = L \ y_std;
            nlml = 0.5*(z'*z) + sum(log(diag(L))) + (n/2)*log(2*pi);
        catch
            nlml = NaN;
        end
        varargout = {nlml};

    otherwise
        error('metric_helpers:unknown', 'Unknown metric name: %s', metric_name);
end
end

function K = build_gram(X, kname, params)
% Build Gram matrix K(X,X) for kernel kname with params. X is n x d.
n = size(X, 1);
K = zeros(n, n);
if contains(kname, 'squared') || contains(kname, 'squaredexponential')
    % [SigmaL, SigmaF]
    ell = params(1); sf = params(2);
    for i = 1:n
        for j = i:n
            r2 = sum((X(i,:) - X(j,:)).^2);
            K(i,j) = sf^2 * exp(-0.5*r2/ell^2);
            K(j,i) = K(i,j);
        end
    end
elseif contains(kname, 'exponential') && ~contains(kname, 'matern')
    % Matern 1/2: [SigmaL, SigmaF]
    ell = params(1); sf = params(2);
    for i = 1:n
        for j = i:n
            r = sqrt(sum((X(i,:) - X(j,:)).^2));
            K(i,j) = sf^2 * exp(-r/ell);
            K(j,i) = K(i,j);
        end
    end
elseif contains(kname, 'matern32')
    ell = params(1); sf = params(2);
    for i = 1:n
        for j = i:n
            r = sqrt(sum((X(i,:) - X(j,:)).^2));
            r_ell = sqrt(3)*r/ell;
            K(i,j) = sf^2 * (1 + r_ell) * exp(-r_ell);
            K(j,i) = K(i,j);
        end
    end
elseif contains(kname, 'matern52')
    ell = params(1); sf = params(2);
    for i = 1:n
        for j = i:n
            r = sqrt(sum((X(i,:) - X(j,:)).^2));
            r_ell = sqrt(5)*r/ell;
            K(i,j) = sf^2 * (1 + r_ell + r_ell^2/3) * exp(-r_ell);
            K(j,i) = K(i,j);
        end
    end
elseif contains(kname, 'rational')
    % [SigmaL, Alpha, SigmaF]
    ell = params(1); alpha = params(2); sf = params(3);
    for i = 1:n
        for j = i:n
            r2 = sum((X(i,:) - X(j,:)).^2);
            K(i,j) = sf^2 * (1 + r2/(2*alpha*ell^2))^(-alpha);
            K(j,i) = K(i,j);
        end
    end
elseif contains(kname, 'periodic') || (contains(kname, 'custom') && numel(params) == 3)
    % Periodic: theta(1)=sigmaF, theta(2)=period, theta(3)=lengthScale
    sigmaF = params(1); p = params(2); l = params(3);
    for i = 1:n
        for j = i:n
            dist = abs(X(i,:) - X(j,:));
            if numel(dist) > 1, dist = sqrt(sum(dist.^2)); end
            K(i,j) = sigmaF^2 * exp(-2 * sin(pi * dist / p).^2 / l^2);
            K(j,i) = K(i,j);
        end
    end
else
    error('metric_helpers:kernel', 'Unsupported kernel: %s', kname);
end
end
