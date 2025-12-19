function K = composite_kernel(X1, X2, params)
% Composite kernel: Matern32 + RationalQuadratic
% This combines local smoothness (Matern) with multi-scale behavior (RQ)
%
% K = w1 * Matern32 + w2 * RationalQuadratic
%
% Inputs:
%   X1, X2: Input vectors
%   params: [l1, sigma_f1, l2, alpha, sigma_f2, w1, w2]
%           where w1 + w2 = 1 (weights for each kernel component)
%
% Note: This is a custom kernel function for fitrgp
% Usage: gprMdl = fitrgp(X, y, 'KernelFunction', @(X1, X2) composite_kernel(X1, X2, params), ...)

    % Extract parameters
    l1 = params(1);        % Matern32 length scale
    sigma_f1 = params(2);   % Matern32 signal std
    l2 = params(3);         % RQ length scale
    alpha = params(4);      % RQ shape parameter
    sigma_f2 = params(5);   % RQ signal std
    w1 = params(6);         % Weight for Matern32
    w2 = params(7);        % Weight for RQ
    
    % Ensure weights sum to 1
    w_total = w1 + w2;
    w1 = w1 / w_total;
    w2 = w2 / w_total;
    
    % Compute Matern32 kernel
    sqrt3 = sqrt(3);
    d = abs(X1 - X2);
    r = d / l1;
    K_matern = sigma_f1^2 * (1 + sqrt3 * r) .* exp(-sqrt3 * r);
    
    % Compute Rational Quadratic kernel
    d2 = (X1 - X2).^2;
    K_rq = sigma_f2^2 * (1 + d2 / (2 * alpha * l2^2)).^(-alpha);
    
    % Composite kernel
    K = w1 * K_matern + w2 * K_rq;
end

