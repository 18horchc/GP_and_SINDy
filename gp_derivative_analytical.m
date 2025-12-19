function [deriv_mean, deriv_var] = gp_derivative_analytical(gprMdl, X_new)
% Compute analytical derivative of a Gaussian Process
%
% For a GP with squared exponential kernel k(x, x') = σ² exp(-(x-x')²/(2l²))
% The derivative GP has covariance:
%   k'(x, x') = σ²/l² * (1 - (x-x')²/l²) * exp(-(x-x')²/(2l²))
%
% Inputs:
%   gprMdl: Fitted GP model from fitrgp
%   X_new:  New input points where to evaluate derivative (column vector)
%
% Outputs:
%   deriv_mean: Mean of derivative GP at X_new
%   deriv_var:  Variance of derivative GP at X_new (optional)

    % Extract training data
    X_train = gprMdl.X;
    y_train = gprMdl.Y;
    
    % Get kernel parameters
    % For squared exponential kernel, we need sigma and length scale
    kernelParams = gprMdl.KernelInformation.KernelParameters;
    
    % Extract parameters
    % MATLAB stores them as [SigmaL, SigmaF] where:
    %   SigmaL = length scale (l)
    %   SigmaF = signal standard deviation (sigma_f)
    if length(kernelParams) >= 2
        l = kernelParams(1);         % Length scale (SigmaL)
        sigma_f = kernelParams(2);  % Signal standard deviation (SigmaF)
    else
        error('Expected at least 2 kernel parameters for squared exponential');
    end
    
    % Get noise variance
    sigma_n = gprMdl.Sigma;  % Observation noise
    
    % Number of training and test points
    n_train = size(X_train, 1);
    n_test = size(X_new, 1);
    
    % Compute covariance matrices using vectorized operations for efficiency
    % K(X, X): covariance between training points
    % k(x, x') = σ² exp(-(x-x')²/(2l²))
    X_train_mat = repmat(X_train, 1, n_train);
    X_train_mat_T = X_train_mat';
    d_XX = X_train_mat - X_train_mat_T;
    K_XX = sigma_f^2 * exp(-d_XX.^2 / (2*l^2));
    K_XX = K_XX + sigma_n^2 * eye(n_train);  % Add noise
    
    % K'(X_new, X): cross-covariance between derivative at test points 
    % and function at training points
    % ∂k(x*, x)/∂x* = -σ²/l² * (x* - x) * exp(-(x*-x)²/(2l²))
    X_new_mat = repmat(X_new, 1, n_train);
    X_train_mat_for_cross = repmat(X_train', n_test, 1);
    d_deriv_X = X_new_mat - X_train_mat_for_cross;
    K_deriv_X = -sigma_f^2 / l^2 * d_deriv_X .* exp(-d_deriv_X.^2 / (2*l^2));
    
    % K'(X_new, X_new): covariance of derivative at test points
    % ∂²k(x*, x*')/(∂x* ∂x*') = σ²/l² * (1 - (x*-x*')²/l²) * exp(-(x*-x*')²/(2l²))
    X_new_mat_2 = repmat(X_new, 1, n_test);
    X_new_mat_2_T = X_new_mat_2';
    d_deriv_XX = X_new_mat_2 - X_new_mat_2_T;
    K_deriv_XX = sigma_f^2 / l^2 * (1 - d_deriv_XX.^2 / l^2) .* exp(-d_deriv_XX.^2 / (2*l^2));
    
    % Compute mean of derivative GP
    % μ' = K'(X_new, X) * K(X,X)^(-1) * y
    alpha = K_XX \ y_train;
    deriv_mean = K_deriv_X * alpha;
    
    % Compute variance if requested
    if nargout > 1
        % Var[f'] = K'(X_new, X_new) - K'(X_new, X) * K(X,X)^(-1) * K'(X, X_new)
        K_deriv_X_transpose = K_deriv_X';
        deriv_var = diag(K_deriv_XX - K_deriv_X * (K_XX \ K_deriv_X_transpose));
    end
end

