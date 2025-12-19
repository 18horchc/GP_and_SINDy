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
    kernelParams = gprMdl.KernelInformation.KernelParameters;
    l = kernelParams(1);         % Length scale
    sigma_f = kernelParams(2);   % Signal standard deviation
    sigma_n = gprMdl.Sigma;      % Observation noise
    
    n_train = size(X_train, 1);
    n_test = size(X_new, 1);
    
    % Compute full covariance matrix K(X_new, X_new)
    % k(x, x') = σ² exp(-(x-x')²/(2l²))
    X_new_mat = repmat(X_new, 1, n_test);
    X_new_mat_T = X_new_mat';
    d_XX = X_new_mat - X_new_mat_T;
    K_XX_new = sigma_f^2 * exp(-d_XX.^2 / (2*l^2));
    
    % Compute cross-covariance K(X_new, X_train)
    X_new_mat_cross = repmat(X_new, 1, n_train);
    X_train_mat_cross = repmat(X_train', n_test, 1);
    d_cross = X_new_mat_cross - X_train_mat_cross;
    K_cross = sigma_f^2 * exp(-d_cross.^2 / (2*l^2));
    
    % Compute K(X_train, X_train)
    X_train_mat = repmat(X_train, 1, n_train);
    X_train_mat_T = X_train_mat';
    d_train = X_train_mat - X_train_mat_T;
    K_train = sigma_f^2 * exp(-d_train.^2 / (2*l^2)) + sigma_n^2 * eye(n_train);
    
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

