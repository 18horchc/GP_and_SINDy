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
%   dXdt_aug: Augmented derivatives (from GP analytical derivatives)
%   augmentation_stats: Structure with statistics

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
    
    % Compute GP analytical derivatives for augmented data
    % We'll compute derivatives for each GP sample
    dXdt_gp_samples = zeros(n_dense, num_vars, num_gp_samples);
    
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
        
        for s = 1:num_gp_samples
            % For derivatives, we use the GP mean derivative (same for all samples)
            % This is because derivatives are computed analytically from the GP model
            for j = 1:n_dense
                t_star = t_dense(j);
                dist = t_star - X_train;
                k_star = sigma_f^2 * exp(-0.5 * (dist.^2) / (L^2));
                dk_dt = -(dist / (L^2)) .* k_star;
                dXdt_gp_samples(j, i, s) = dk_dt' * alpha;
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
                dXdt_gp_combined = dXdt_gp_samples(:, :, 1);
            else
                % Average multiple samples
                X_gp_combined = mean(X_gp_samples, 3);
                dXdt_gp_combined = mean(dXdt_gp_samples, 3);
            end
            
            % Combine: original sparse + GP samples
            t_aug = [t_sparse; t_dense];
            X_aug = [X_sparse; X_gp_combined];
            
            % For derivatives: use finite differences on sparse, GP analytical on dense
            dXdt_sparse_fd = compute_finite_differences(t_sparse, X_sparse, 'central');
            dXdt_aug = [dXdt_sparse_fd; dXdt_gp_combined];
            
        case 'replace'
            % Replace sparse data entirely with GP samples
            t_aug = t_dense;
            if num_gp_samples == 1
                X_aug = X_gp_samples(:, :, 1);
                dXdt_aug = dXdt_gp_samples(:, :, 1);
            else
                X_aug = mean(X_gp_samples, 3);
                dXdt_aug = mean(dXdt_gp_samples, 3);
            end
            
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

