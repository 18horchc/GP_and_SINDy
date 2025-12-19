%% 4b. Sequential Thresholded Least Squares (STLSQ)
function Xi = pool_STLSQ(Theta, dXdt, lambda, iterations)
    % Solves Theta * Xi = dXdt with sparsity
    % Initial guess: Standard Least Squares
    Xi = Theta \ dXdt;
    
    % Sequentially threshold small coefficients
    for k = 1:iterations
        small_indices = abs(Xi) < lambda;
        Xi(small_indices) = 0;
        
        for j = 1:size(dXdt, 2) % For each state variable
            big_indices = ~small_indices(:, j);
            % Regress only on the significant indices
            % Use ridge regression (small alpha) to handle the rank deficiency
            lambda_ridge = 1e-6; 
            T_sub = Theta(:, big_indices);
            Xi(big_indices, j) = (T_sub' * T_sub + lambda_ridge * eye(size(T_sub,2))) \ (T_sub' * dXdt(:, j));        end
    end
end