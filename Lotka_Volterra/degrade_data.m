function [t_sparse, X_noisy] = degrade_data(t_true, X_true, n_points, noise_std)
    % DEGRADE_DATA Simulates real-world data collection.
    %
    % Inputs:
    %   t_true    - High-resolution time vector
    %   X_true    - High-resolution state data
    %   n_points  - The number of samples to keep (integer)
    %   noise_std - Standard deviation of the Gaussian noise to add
    
    % 1. Determine sampling indices
    % We pick n_points equally spaced across the time range
    idx = round(linspace(1, length(t_true), n_points));
    
    t_sparse = t_true(idx);
    X_clean_sub = X_true(idx, :);
    
    % 2. Add Gaussian White Noise
    % Noise is relative to each channel (x and y)
    noise = noise_std * randn(size(X_clean_sub));
    X_noisy = X_clean_sub + noise;
    
    % Ensure no negative population values (physical constraint for LV)
    X_noisy(X_noisy < 0) = 0;
end