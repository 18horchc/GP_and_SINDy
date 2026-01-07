function [t_sparse, X_sparse, X_true_at_sparse, noise_level] = generate_sparse_noisy_data(t_true, X_true, n_points, noise_fraction, sampling_method)
%GENERATE_SPARSE_NOISY_DATA Generate sparse noisy observations from true solution
%
% Inputs:
%   t_true: Full time vector
%   X_true: Full state matrix [x(t), y(t)] of size [length(t), 2]
%   n_points: Number of sparse points to sample (default: 8-10)
%   noise_fraction: Noise level as fraction of signal std (default: 0.1 = 10%)
%   sampling_method: 'uniform', 'random', or 'strategic' (default: 'uniform')
%
% Outputs:
%   t_sparse: Time points where data is sampled
%   X_sparse: Noisy observations [x_noisy, y_noisy] of size [n_points, 2]
%   X_true_at_sparse: True values at sparse points (for reference)
%   noise_level: Actual noise standard deviation used

    if nargin < 5
        sampling_method = 'uniform';
    end
    if nargin < 4
        noise_fraction = 0.1;  % 10% noise
    end
    if nargin < 3
        n_points = 9;  % Default: 9 points (between 8-10)
    end
    
    % Sample time points based on method
    switch lower(sampling_method)
        case 'uniform'
            % Uniform spacing
            idx_sparse = round(linspace(1, length(t_true), n_points));
            idx_sparse = unique(idx_sparse);  % Ensure unique indices
            if length(idx_sparse) < n_points
                % If we lost points due to rounding, add more
                idx_sparse = round(linspace(1, length(t_true), n_points + 2));
                idx_sparse = unique(idx_sparse);
                idx_sparse = idx_sparse(1:min(n_points, length(idx_sparse)));
            end
            
        case 'random'
            % Random sampling
            idx_sparse = sort(randperm(length(t_true), n_points));
            
        case 'strategic'
            % More points where dynamics change rapidly (higher curvature)
            % Compute approximate curvature as |d^2x/dt^2|
            dt = mean(diff(t_true));
            dXdt = diff(X_true, 1) / dt;
            d2Xdt2 = diff(dXdt, 1) / dt;
            curvature = sqrt(sum(d2Xdt2.^2, 2));
            curvature = [curvature(1); curvature; curvature(end)];  % Pad
            
            % Sample with probability proportional to curvature
            prob = curvature / sum(curvature);
            cumprob = cumsum(prob);
            idx_sparse = zeros(n_points, 1);
            for i = 1:n_points
                r = rand();
                idx_sparse(i) = find(cumprob >= r, 1, 'first');
            end
            idx_sparse = unique(sort(idx_sparse));
            % If we have fewer than n_points, fill with uniform
            if length(idx_sparse) < n_points
                remaining = n_points - length(idx_sparse);
                all_idx = 1:length(t_true);
                available = setdiff(all_idx, idx_sparse);
                if ~isempty(available)
                    additional = available(round(linspace(1, length(available), remaining)));
                    idx_sparse = sort([idx_sparse(:); additional(:)]);
                end
            end
            idx_sparse = idx_sparse(1:min(n_points, length(idx_sparse)));
            
        otherwise
            error('Unknown sampling method: %s. Use ''uniform'', ''random'', or ''strategic''', sampling_method);
    end
    
    % Extract sparse time and true values
    t_sparse = t_true(idx_sparse);
    X_true_at_sparse = X_true(idx_sparse, :);
    
    % Compute noise level based on signal standard deviation
    signal_std = std(X_true, 0, 1);  % Standard deviation along time dimension
    noise_std = noise_fraction * signal_std;
    noise_level = noise_std;
    
    % Add Gaussian noise
    X_sparse = X_true_at_sparse + randn(size(X_true_at_sparse)) .* noise_std;
end

