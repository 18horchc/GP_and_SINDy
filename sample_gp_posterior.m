function [samples] = sample_gp_posterior(gprMdl, X_new, num_samples)
% Sample curves from GP posterior distribution
%
% This properly samples from the multivariate Gaussian posterior of the GP,
% giving us actual GP curves (not just noisy mean estimates).
%
% Inputs:
%   gprMdl: Fitted GP model from fitrgp
%   X_new:  New input points where to sample (column vector)
%   num_samples: Number of sample curves to generate (default: 1)
%
% Output:
%   samples: Matrix of size [length(X_new), num_samples] where each column
%            is one sample curve from the GP posterior

    if nargin < 3
        num_samples = 1;
    end

    % Get GP predictions (mean and standard deviation)
    [mu, sigma] = predict(gprMdl, X_new);
    
    % For proper GP sampling, we need the full covariance matrix, not just
    % the diagonal variance. However, predict() only gives us the diagonal.
    % We need to compute the full covariance matrix.
    
    % Extract training data and parameters
    X_train = gprMdl.X;
    y_train = gprMdl.Y;
    kernelName = gprMdl.KernelInformation.Name;
    kernelParams = gprMdl.KernelInformation.KernelParameters;
    sigma_n = gprMdl.Sigma;      % Observation noise
    
    n_train = size(X_train, 1);
    n_test = size(X_new, 1);
    
    % Compute covariance matrices based on kernel type
    X_new_mat = repmat(X_new, 1, n_test);
    X_new_mat_T = X_new_mat';
    d_XX = X_new_mat - X_new_mat_T;
    
    X_new_mat_cross = repmat(X_new, 1, n_train);
    X_train_mat_cross = repmat(X_train', n_test, 1);
    d_cross = X_new_mat_cross - X_train_mat_cross;
    
    X_train_mat = repmat(X_train, 1, n_train);
    X_train_mat_T = X_train_mat';
    d_train = X_train_mat - X_train_mat_T;
    
    if strcmpi(kernelName, 'SquaredExponential')
        l = kernelParams(1);
        sigma_f = kernelParams(2);
        K_XX_new = sigma_f^2 * exp(-d_XX.^2 / (2*l^2));
        K_cross = sigma_f^2 * exp(-d_cross.^2 / (2*l^2));
        K_train = sigma_f^2 * exp(-d_train.^2 / (2*l^2)) + sigma_n^2 * eye(n_train);
        
    elseif strcmpi(kernelName, 'Matern32')
        l = kernelParams(1);
        sigma_f = kernelParams(2);
        sqrt3 = sqrt(3);
        r_XX = abs(d_XX);
        r_cross = abs(d_cross);
        r_train = abs(d_train);
        K_XX_new = sigma_f^2 * (1 + sqrt3 * r_XX / l) .* exp(-sqrt3 * r_XX / l);
        K_cross = sigma_f^2 * (1 + sqrt3 * r_cross / l) .* exp(-sqrt3 * r_cross / l);
        K_train = sigma_f^2 * (1 + sqrt3 * r_train / l) .* exp(-sqrt3 * r_train / l) + sigma_n^2 * eye(n_train);
        
    elseif strcmpi(kernelName, 'RationalQuadratic')
        % RQ kernel: k(x,x') = σ²(1 + (x-x')²/(2αl²))^(-α)
        % MATLAB stores as [SigmaL, AlphaRQ, SigmaF]
        l = kernelParams(1);         % Length scale (SigmaL)
        alpha = kernelParams(2);      % Shape parameter (AlphaRQ)
        sigma_f = kernelParams(3);    % Signal std (SigmaF)
        d2_XX = d_XX.^2;
        d2_cross = d_cross.^2;
        d2_train = d_train.^2;
        K_XX_new = sigma_f^2 * (1 + d2_XX / (2*alpha*l^2)).^(-alpha);
        K_cross = sigma_f^2 * (1 + d2_cross / (2*alpha*l^2)).^(-alpha);
        K_train = sigma_f^2 * (1 + d2_train / (2*alpha*l^2)).^(-alpha) + sigma_n^2 * eye(n_train);
        
    else
        warning('Kernel %s not fully supported in sample_gp_posterior. Using predict() approximation.', kernelName);
        % Fallback: use diagonal approximation from predict()
        [mu, sigma] = predict(gprMdl, X_new);
        samples = zeros(n_test, num_samples);
        for i = 1:num_samples
            samples(:, i) = mu + sigma .* randn(n_test, 1);
        end
        return;
    end
    
    % Posterior covariance at X_new
    % Cov[f*|X,y] = K(X*, X*) - K(X*, X) * K(X,X)^(-1) * K(X, X*)
    K_cross_T = K_cross';
    K_inv = K_train \ eye(n_train);
    K_posterior = K_XX_new - K_cross * (K_inv * K_cross_T);
    
    % Ensure K_posterior is symmetric and positive definite (numerical stability)
    K_posterior = (K_posterior + K_posterior') / 2;
    K_posterior = K_posterior + 1e-6 * eye(n_test);  % Add small jitter
    
    % Sample from multivariate Gaussian: N(mu, K_posterior)
    % Use Cholesky decomposition for numerical stability
    try
        L = chol(K_posterior, 'lower');
    catch
        % If Cholesky fails, use eigendecomposition
        [V, D] = eig(K_posterior);
        D = diag(max(real(diag(D)), 1e-10));  % Ensure positive
        L = V * sqrt(D);
    end
    
    % Generate samples
    samples = zeros(n_test, num_samples);
    for i = 1:num_samples
        z = randn(n_test, 1);
        samples(:, i) = mu + L * z;
    end
end

