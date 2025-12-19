function [deriv_mean, deriv_var] = gp_derivative_analytical(gprMdl, X_new)
% Compute analytical derivative of a Gaussian Process
%
% Supports multiple kernel types:
% - Squared exponential: k(x, x') = σ² exp(-(x-x')²/(2l²))
% - Matern32: k(r) = σ²(1 + √3r/l)exp(-√3r/l) where r = |x-x'|
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
    
    % Get kernel type and parameters
    kernelName = gprMdl.KernelInformation.Name;
    kernelParams = gprMdl.KernelInformation.KernelParameters;
    
    % Extract parameters (MATLAB stores as [SigmaL, SigmaF])
    if length(kernelParams) >= 2
        l = kernelParams(1);         % Length scale (SigmaL)
        sigma_f = kernelParams(2);  % Signal standard deviation (SigmaF)
    else
        error('Expected at least 2 kernel parameters');
    end
    
    % Get noise variance
    sigma_n = gprMdl.Sigma;  % Observation noise
    
    % Number of training and test points
    n_train = size(X_train, 1);
    n_test = size(X_new, 1);
    
    % Compute covariance matrices based on kernel type
    X_train_mat = repmat(X_train, 1, n_train);
    X_train_mat_T = X_train_mat';
    d_XX = X_train_mat - X_train_mat_T;
    
    X_new_mat = repmat(X_new, 1, n_train);
    X_train_mat_for_cross = repmat(X_train', n_test, 1);
    d_deriv_X = X_new_mat - X_train_mat_for_cross;
    
    X_new_mat_2 = repmat(X_new, 1, n_test);
    X_new_mat_2_T = X_new_mat_2';
    d_deriv_XX = X_new_mat_2 - X_new_mat_2_T;
    
    if strcmpi(kernelName, 'SquaredExponential')
        % Squared exponential kernel
        K_XX = sigma_f^2 * exp(-d_XX.^2 / (2*l^2));
        K_deriv_X = -sigma_f^2 / l^2 * d_deriv_X .* exp(-d_deriv_X.^2 / (2*l^2));
        K_deriv_XX = sigma_f^2 / l^2 * (1 - d_deriv_XX.^2 / l^2) .* exp(-d_deriv_XX.^2 / (2*l^2));
        
    elseif strcmpi(kernelName, 'Matern32')
        % Matern 3/2 kernel: k(r) = σ²(1 + √3r/l)exp(-√3r/l)
        sqrt3 = sqrt(3);
        r_XX = abs(d_XX);
        r_deriv_X = abs(d_deriv_X);
        r_deriv_XX = abs(d_deriv_XX);
        
        K_XX = sigma_f^2 * (1 + sqrt3 * r_XX / l) .* exp(-sqrt3 * r_XX / l);
        % Derivative of Matern32: more complex, using numerical approximation for cross-covariance
        % For exact derivative, we'd need: ∂k/∂x* = -σ²√3/l * sign(x*-x) * (1 + √3r/l) * exp(-√3r/l)
        K_deriv_X = -sigma_f^2 * sqrt3 / l * sign(d_deriv_X) .* (1 + sqrt3 * r_deriv_X / l) .* exp(-sqrt3 * r_deriv_X / l);
        % Second derivative for K_deriv_XX (simplified)
        K_deriv_XX = sigma_f^2 * 3 / l^2 * (1 - sqrt3 * r_deriv_XX / l) .* exp(-sqrt3 * r_deriv_XX / l);
        
    elseif strcmpi(kernelName, 'RationalQuadratic')
        % Rational Quadratic kernel: k(x,x') = σ²(1 + (x-x')²/(2αl²))^(-α)
        % MATLAB stores as [SigmaL, AlphaRQ, SigmaF]
        if length(kernelParams) >= 3
            l = kernelParams(1);         % Length scale (SigmaL)
            alpha = kernelParams(2);     % Shape parameter (AlphaRQ)
            sigma_f = kernelParams(3);   % Signal std (SigmaF)
        else
            error('RationalQuadratic kernel requires 3 parameters: [l, alpha, sigma_f]');
        end
        
        d2_XX = d_XX.^2;
        d2_deriv_X = d_deriv_X.^2;
        d2_deriv_XX = d_deriv_XX.^2;
        
        K_XX = sigma_f^2 * (1 + d2_XX / (2*alpha*l^2)).^(-alpha);
        % Derivative of RQ kernel: ∂k/∂x* = -σ² * (x*-x) / l² * (1 + (x*-x)²/(2αl²))^(-α-1)
        % Note: The α cancels in the numerator, leaving just (x*-x)/l²
        K_deriv_X = -sigma_f^2 * d_deriv_X / l^2 .* (1 + d2_deriv_X / (2*alpha*l^2)).^(-alpha-1);
        % Second derivative: ∂²k/(∂x*∂x*') = σ²/l² * [1 - (α+1)(x*-x*')²/(2αl²)] * (1 + (x*-x*')²/(2αl²))^(-α-1)
        K_deriv_XX = sigma_f^2 / l^2 * (1 - (alpha+1)*d2_deriv_XX / (2*alpha*l^2)) .* (1 + d2_deriv_XX / (2*alpha*l^2)).^(-alpha-1);
        
    else
        warning('Kernel %s not fully supported for analytical derivatives. Using finite differences approximation.', kernelName);
        % Fallback: use finite differences on GP mean
        [mu, ~] = predict(gprMdl, X_new);
        dt = mean(diff(X_new));
        if length(X_new) > 1
            deriv_mean = gradient(mu, dt);
        else
            deriv_mean = 0;
        end
        if nargout > 1
            deriv_var = zeros(size(deriv_mean));
        end
        return;
    end
    
    K_XX = K_XX + sigma_n^2 * eye(n_train);  % Add noise
    
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

