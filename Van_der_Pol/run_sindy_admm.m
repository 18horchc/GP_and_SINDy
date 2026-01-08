function [Xi, stats] = run_sindy_admm(Theta, dXdt, lambda, options)
%RUN_SINDY_ADMM Run ADMM (Alternating Direction Method of Multipliers) for LASSO
%
% Solves the LASSO problem: min ||ΘΞ - Ẋ||² + λ||Ξ||₁
% Using ADMM algorithm as in Hsin et al. (2025)
%
% Inputs:
%   Theta: Library matrix [N x num_terms]
%   dXdt: Derivative matrix [N x num_states]
%   lambda: Sparsity regularization parameter (default: 0.01)
%   options: (optional) Structure with options:
%     - rho: ADMM penalty parameter (default: 1.0)
%     - max_iterations: Maximum ADMM iterations (default: 1000)
%     - abs_tol: Absolute tolerance for convergence (default: 1e-4)
%     - rel_tol: Relative tolerance for convergence (default: 1e-2)
%     - verbose: Print convergence info (default: false)
%
% Outputs:
%   Xi: Coefficient matrix [num_terms x num_states]
%   stats: Structure with statistics (sparsity, iterations, residuals, etc.)

    if nargin < 4
        options = struct();
    end
    if nargin < 3
        lambda = 0.01;
    end
    
    % Set default options
    if ~isfield(options, 'rho')
        options.rho = 1.0;  % ADMM penalty parameter
    end
    if ~isfield(options, 'max_iterations')
        options.max_iterations = 1000;
    end
    if ~isfield(options, 'abs_tol')
        options.abs_tol = 1e-4;
    end
    if ~isfield(options, 'rel_tol')
        options.rel_tol = 1e-2;
    end
    if ~isfield(options, 'verbose')
        options.verbose = false;
    end
    
    [N, num_terms] = size(Theta);
    [~, num_states] = size(dXdt);
    
    % Initialize variables
    % Start with least squares solution
    Xi = Theta \ dXdt;
    z = Xi;  % Auxiliary variable
    u = zeros(size(Xi));  % Dual variable
    
    % Precompute (Θ'Θ + ρI)^(-1) for efficiency (same for all state variables)
    rho = options.rho;
    ThetaT_Theta = Theta' * Theta;
    A = ThetaT_Theta + rho * eye(num_terms);
    
    % Check if A is well-conditioned, otherwise use pseudo-inverse
    if cond(A) > 1e12
        A_inv = pinv(A);
    else
        A_inv = inv(A);
    end
    
    % Store initial sparsity
    initial_sparsity = sum(Xi(:) == 0) / numel(Xi);
    
    % ADMM iterations
    converged = false;
    for iter = 1:options.max_iterations
        % Store previous values for convergence check
        Xi_prev = Xi;
        z_prev = z;
        
        % Update Xi: minimize ||ΘΞ - Ẋ||² + (ρ/2)||Ξ - z + u||²
        % Solution: Ξ = (Θ'Θ + ρI)^(-1) * (Θ'Ẋ + ρ(z - u))
        for j = 1:num_states
            b = Theta' * dXdt(:, j) + rho * (z(:, j) - u(:, j));
            Xi(:, j) = A_inv * b;
        end
        
        % Update z: soft thresholding
        % z = soft_threshold(Ξ + u, λ/ρ)
        threshold = lambda / rho;
        z = soft_threshold(Xi + u, threshold);
        
        % Update u: dual variable update
        u = u + Xi - z;
        
        % Check convergence
        % Primal residual: r = ||Ξ - z||
        % Dual residual: s = ρ||z - z_prev||
        r_norm = norm(Xi - z, 'fro');
        s_norm = rho * norm(z - z_prev, 'fro');
        
        % Tolerance checks
        eps_pri = sqrt(num_terms * num_states) * options.abs_tol + ...
                  options.rel_tol * max(norm(Xi, 'fro'), norm(z, 'fro'));
        eps_dual = sqrt(num_terms * num_states) * options.abs_tol + ...
                   options.rel_tol * norm(rho * u, 'fro');
        
        if r_norm < eps_pri && s_norm < eps_dual
            converged = true;
            if options.verbose
                fprintf('  ADMM converged at iteration %d\n', iter);
            end
            break;
        end
        
        if options.verbose && mod(iter, 100) == 0
            fprintf('  Iter %d: r_norm=%.2e, s_norm=%.2e\n', iter, r_norm, s_norm);
        end
    end
    
    if ~converged && options.verbose
        warning('ADMM did not converge within %d iterations', options.max_iterations);
    end
    
    % Compute statistics
    final_sparsity = sum(Xi(:) == 0) / numel(Xi);
    num_nonzero = sum(Xi(:) ~= 0);
    
    % Compute prediction error
    dXdt_pred = Theta * Xi;
    prediction_error = dXdt - dXdt_pred;
    rmse = sqrt(mean(prediction_error(:).^2));
    
    stats = struct();
    stats.initial_sparsity = initial_sparsity;
    stats.final_sparsity = final_sparsity;
    stats.num_nonzero = num_nonzero;
    stats.num_terms = num_terms;
    stats.num_states = num_states;
    stats.lambda = lambda;
    stats.rho = rho;
    stats.iterations = iter;
    stats.converged = converged;
    stats.rmse = rmse;
    stats.r_norm = r_norm;
    stats.s_norm = s_norm;
end

function z = soft_threshold(x, threshold)
%SOFT_THRESHOLD Soft thresholding operator
%   z = sign(x) * max(|x| - threshold, 0)
    z = sign(x) .* max(abs(x) - threshold, 0);
end

