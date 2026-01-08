function [t_aug, X_aug, dXdt_aug, augmentation_stats] = augment_data_with_gp_samples(...
    t_sparse, X_sparse, gp_models, t_dense, num_gp_samples, options)
%AUGMENT_DATA_WITH_GP_SAMPLES Augment sparse data with GP posterior samples
%
% Inputs:
%   t_sparse: Sparse time points [N x 1]
%   X_sparse: Sparse noisy observations [N x 2]
%   gp_models: Cell array of GP models {gprMdl_x, gprMdl_y}
%   t_dense: Dense time grid for GP sampling [M x 1]
%   num_gp_samples: Number of GP sample curves (default: 1)
%   options: (optional) Structure with options:
%     - combine_method: 'append' (default) or 'replace'
%     - weight_original: Weight for original data (default: 1.0)
%
% Outputs:
%   t_aug: Augmented time vector
%   X_aug: Augmented state matrix
%   dXdt_aug: Augmented derivatives (from GP analytical derivatives for ALL points)
%   augmentation_stats: Structure with statistics
%
% Note: Following Hsin et al. (2025), this function uses GP analytical derivatives
% for all augmented data points (both sparse and dense), providing consistent
% denoised derivatives throughout rather than mixing finite differences with GP derivatives.

    if nargin < 6
        options = struct();
    end
    if nargin < 5
        num_gp_samples = 1;
    end
    
    % Set default options
    if ~isfield(options, 'combine_method')
        options.combine_method = 'append';  % Append GP samples to original
    end
    if ~isfield(options, 'weight_original')
        options.weight_original = 1.0;  % No weighting by default
    end
    
    num_vars = length(gp_models);
    n_sparse = length(t_sparse);
    n_dense = length(t_dense);
    
    % Sample from GP posterior for each state variable
    X_gp_samples = zeros(n_dense, num_vars, num_gp_samples);
    
    for i = 1:num_vars
        gp_model = gp_models{i};
        
        % Sample from GP posterior
        % Note: We'll use the existing sample_gp_posterior function if available
        % Otherwise, we'll use a simpler approach
        try
            % Try to use sample_gp_posterior if it exists
            if exist('sample_gp_posterior', 'file')
                samples = sample_gp_posterior(gp_model, t_dense, num_gp_samples);
                X_gp_samples(:, i, :) = reshape(samples, [n_dense, 1, num_gp_samples]);
            else
                % Fallback: sample from predictive distribution
                [mu, sigma] = predict(gp_model, t_dense);
                for s = 1:num_gp_samples
                    X_gp_samples(:, i, s) = mu + sigma .* randn(n_dense, 1);
                end
            end
        catch
            % Fallback: use GP mean with added noise
            [mu, sigma] = predict(gp_model, t_dense);
            for s = 1:num_gp_samples
                X_gp_samples(:, i, s) = mu + sigma .* randn(n_dense, 1);
            end
        end
    end
    
    % Compute GP analytical derivatives for ALL points (sparse + dense)
    % Following Hsin et al. (2025): use GP analytical derivatives for all augmented data points
    % This provides consistent, denoised derivatives throughout
    
    % First, determine all time points where we need derivatives
    switch lower(options.combine_method)
        case 'append'
            t_all = [t_sparse; t_dense];
        case 'replace'
            t_all = t_dense;
        otherwise
            t_all = [t_sparse; t_dense];
    end
    
    n_all = length(t_all);
    
    % Compute GP analytical derivatives at all time points
    % Note: Derivatives are computed from the GP model, so they're the same for all samples
    dXdt_gp_all = zeros(n_all, num_vars);
    
    for i = 1:num_vars
        gp_model = gp_models{i};
        alpha = gp_model.Alpha;
        X_train = gp_model.X;
        kernel_params = gp_model.KernelInformation.KernelParameters;
        
        if length(kernel_params) >= 2
            L = kernel_params(1);
            sigma_f = kernel_params(2);
        else
            L = kernel_params(1);
            sigma_f = 1.0;
        end
        
        % Compute GP analytical derivative at all time points
        for j = 1:n_all
            t_star = t_all(j);
            dist = t_star - X_train;
            k_star = sigma_f^2 * exp(-0.5 * (dist.^2) / (L^2));
            dk_dt = -(dist / (L^2)) .* k_star;
            dXdt_gp_all(j, i) = dk_dt' * alpha;
        end
    end
    
    % Sample from GP posterior for dense grid only (for state values)
    X_gp_samples = zeros(n_dense, num_vars, num_gp_samples);
    
    for i = 1:num_vars
        gp_model = gp_models{i};
        
        % Sample from GP posterior
        try
            if exist('sample_gp_posterior', 'file')
                samples = sample_gp_posterior(gp_model, t_dense, num_gp_samples);
                X_gp_samples(:, i, :) = reshape(samples, [n_dense, 1, num_gp_samples]);
            else
                [mu, sigma] = predict(gp_model, t_dense);
                for s = 1:num_gp_samples
                    X_gp_samples(:, i, s) = mu + sigma .* randn(n_dense, 1);
                end
            end
        catch
            [mu, sigma] = predict(gp_model, t_dense);
            for s = 1:num_gp_samples
                X_gp_samples(:, i, s) = mu + sigma .* randn(n_dense, 1);
            end
        end
    end
    
    % Combine sparse data with GP samples
    switch lower(options.combine_method)
        case 'append'
            % Append GP samples to original sparse data
            % Use first GP sample (or average if multiple)
            if num_gp_samples == 1
                X_gp_combined = X_gp_samples(:, :, 1);
            else
                % Average multiple samples
                X_gp_combined = mean(X_gp_samples, 3);
            end
            
            % Combine: original sparse + GP samples
            t_aug = [t_sparse; t_dense];
            X_aug = [X_sparse; X_gp_combined];
            
            % Use GP analytical derivatives for ALL points (sparse + dense)
            % Split the pre-computed derivatives: first n_sparse for sparse, rest for dense
            dXdt_sparse_gp = dXdt_gp_all(1:n_sparse, :);
            dXdt_dense_gp = dXdt_gp_all(n_sparse+1:end, :);
            dXdt_aug = [dXdt_sparse_gp; dXdt_dense_gp];
            
        case 'replace'
            % Replace sparse data entirely with GP samples
            t_aug = t_dense;
            if num_gp_samples == 1
                X_aug = X_gp_samples(:, :, 1);
            else
                X_aug = mean(X_gp_samples, 3);
            end
            
            % Use GP analytical derivatives for all dense points
            % Since t_all = t_dense in this case, dXdt_gp_all contains derivatives for all points
            dXdt_aug = dXdt_gp_all;
            
        otherwise
            error('Unknown combine_method: %s. Use ''append'' or ''replace''', options.combine_method);
    end
    
    % Sort by time
    [t_aug, sort_idx] = sort(t_aug);
    X_aug = X_aug(sort_idx, :);
    dXdt_aug = dXdt_aug(sort_idx, :);
    
    % Remove duplicates (keep first occurrence)
    [t_aug, unique_idx] = unique(t_aug, 'stable');
    X_aug = X_aug(unique_idx, :);
    dXdt_aug = dXdt_aug(unique_idx, :);
    
    % Compute statistics
    augmentation_stats = struct();
    augmentation_stats.n_original = n_sparse;
    augmentation_stats.n_gp_samples = n_dense;
    augmentation_stats.n_augmented = length(t_aug);
    augmentation_stats.num_gp_samples = num_gp_samples;
    augmentation_stats.combine_method = options.combine_method;
end

