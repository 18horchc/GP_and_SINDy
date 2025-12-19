%% 4c. The Ensemble Wrapper
function [Xi_mean, inclusion_prob] = run_ESINDy(X, polyorder, dXdt, lambda, n_models, ensemble_size)
    % X: state data, dXdt: derivative data
    % n_models: number of ensemble realizations (e.g., 100)
    % ensemble_size: fraction of data to use in each bag (e.g., 0.8)
    
    Theta = build_library(X, polyorder); 
    [num_terms, num_states] = size(Theta \ dXdt);
    
    Xi_ensemble = zeros(num_terms, num_states, n_models);
    num_samples = size(X, 1);
    n_subset = round(num_samples * ensemble_size);

    for i = 1:n_models
        % Randomly sample data points with replacement (bootstrapping)
        idx = randi(num_samples, n_subset, 1);
        
        % Solve SINDy for this subset
        Xi_ensemble(:, :, i) = pool_STLSQ(Theta(idx, :), dXdt(idx, :), lambda, 10);
    end
    
    % Calculate Inclusion Probability (how often is a term non-zero?)
    inclusion_prob = sum(Xi_ensemble ~= 0, 3) / n_models;
    
    % Calculate Mean Coefficients (of the non-zero entries)
    Xi_mean = mean(Xi_ensemble, 3);
end