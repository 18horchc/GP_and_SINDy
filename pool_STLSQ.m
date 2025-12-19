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
            Xi(big_indices, j) = Theta(:, big_indices) \ dXdt(:, j);
        end
    end
end