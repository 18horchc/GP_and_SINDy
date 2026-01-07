function [Xi, stats] = run_sindy_stls(Theta, dXdt, lambda, max_iterations)
%RUN_SINDY_STLS Run Sequential Thresholded Least Squares (STLS) for SINDy
%
% Inputs:
%   Theta: Library matrix [N x num_terms]
%   dXdt: Derivative matrix [N x num_states]
%   lambda: Sparsity threshold (default: 0.01)
%   max_iterations: Maximum STLS iterations (default: 10)
%
% Outputs:
%   Xi: Coefficient matrix [num_terms x num_states]
%   stats: Structure with statistics (sparsity, iterations, etc.)

    if nargin < 4
        max_iterations = 10;
    end
    if nargin < 3
        lambda = 0.01;
    end
    
    [N, num_terms] = size(Theta);
    [~, num_states] = size(dXdt);
    
    % Initial guess: Standard least squares
    Xi = Theta \ dXdt;
    
    % Store initial sparsity
    initial_sparsity = sum(Xi(:) == 0) / numel(Xi);
    
    % Sequential thresholded least squares
    for iter = 1:max_iterations
        % Find small coefficients
        small_indices = abs(Xi) < lambda;
        
        % Threshold small coefficients to zero
        Xi(small_indices) = 0;
        
        % For each state variable, regress on non-zero terms only
        for j = 1:num_states
            big_indices = ~small_indices(:, j);
            
            if sum(big_indices) > 0
                % Regress only on significant terms
                % Use ridge regression for numerical stability
                lambda_ridge = 1e-6;
                Theta_sub = Theta(:, big_indices);
                
                if size(Theta_sub, 2) > 0 && rank(Theta_sub) == size(Theta_sub, 2)
                    % Full rank: use standard least squares
                    Xi(big_indices, j) = Theta_sub \ dXdt(:, j);
                else
                    % Rank deficient: use ridge regression
                    Xi(big_indices, j) = (Theta_sub' * Theta_sub + lambda_ridge * eye(size(Theta_sub, 2))) \ ...
                                         (Theta_sub' * dXdt(:, j));
                end
            else
                % All coefficients were thresholded - set to zero
                Xi(:, j) = 0;
            end
        end
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
    stats.iterations = max_iterations;
    stats.rmse = rmse;
end

