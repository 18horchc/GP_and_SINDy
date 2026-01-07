function metrics = compute_comprehensive_metrics(X_true, dXdt_true, X_pred, dXdt_pred, ...
    Xi_pred, Xi_true, library_names, method_name)
%COMPUTE_COMPREHENSIVE_METRICS Compute comprehensive metrics for a method
%
% Inputs:
%   X_true: True state matrix [N x 2]
%   dXdt_true: True derivatives [N x 2]
%   X_pred: Predicted state matrix [N x 2]
%   dXdt_pred: Predicted derivatives [N x 2] (can be empty if N/A)
%   Xi_pred: Predicted coefficients [num_terms x 2] (can be empty if N/A)
%   Xi_true: True coefficients [num_terms x 2] (can be empty if N/A)
%   library_names: Library term names
%   method_name: Name of method (for display)
%
% Outputs:
%   metrics: Structure with all computed metrics

    metrics = struct();
    metrics.method_name = method_name;
    
    % Trajectory metrics
    traj_error = X_pred - X_true;
    metrics.rmse_traj = sqrt(mean(traj_error(:).^2));
    metrics.mae_traj = mean(abs(traj_error(:)));
    
    % R² for trajectory
    ss_res = sum(traj_error(:).^2);
    ss_tot = sum((X_true(:) - mean(X_true(:))).^2);
    if ss_tot > 0
        metrics.r2_traj = 1 - ss_res / ss_tot;
    else
        metrics.r2_traj = NaN;
    end
    
    % Per-state metrics
    for i = 1:size(X_true, 2)
        state_name = ['state_' num2str(i)];
        metrics.(['rmse_' state_name]) = sqrt(mean(traj_error(:, i).^2));
        metrics.(['mae_' state_name]) = mean(abs(traj_error(:, i)));
        
        ss_res_i = sum(traj_error(:, i).^2);
        ss_tot_i = sum((X_true(:, i) - mean(X_true(:, i))).^2);
        if ss_tot_i > 0
            metrics.(['r2_' state_name]) = 1 - ss_res_i / ss_tot_i;
        else
            metrics.(['r2_' state_name]) = NaN;
        end
    end
    
    % Derivative metrics (if available)
    if ~isempty(dXdt_pred) && ~isempty(dXdt_true)
        deriv_error = dXdt_pred - dXdt_true;
        metrics.rmse_deriv = sqrt(mean(deriv_error(:).^2));
        metrics.mae_deriv = mean(abs(deriv_error(:)));
        
        % Correlation coefficient
        if numel(dXdt_true) > 1
            corr_coef = corrcoef(dXdt_true(:), dXdt_pred(:));
            metrics.corr_deriv = corr_coef(1, 2);
        else
            metrics.corr_deriv = NaN;
        end
    else
        metrics.rmse_deriv = NaN;
        metrics.mae_deriv = NaN;
        metrics.corr_deriv = NaN;
    end
    
    % Coefficient metrics (if available)
    if ~isempty(Xi_pred) && ~isempty(Xi_true)
        coeff_error = Xi_pred - Xi_true;
        metrics.rmse_coeff = sqrt(mean(coeff_error(:).^2));
        metrics.mae_coeff = mean(abs(coeff_error(:)));
        metrics.coeff_error_l2 = norm(coeff_error(:));
        
        % Count correctly identified terms (non-zero in both or zero in both)
        nz_true = (Xi_true ~= 0);
        nz_pred = (Xi_pred ~= 0);
        correct_terms = sum((nz_true == nz_pred), 'all');
        total_terms = numel(Xi_true);
        metrics.coeff_identification_rate = correct_terms / total_terms;
        
        % Count sparsity (number of non-zero terms)
        metrics.sparsity_pred = sum(Xi_pred(:) ~= 0);
        metrics.sparsity_true = sum(Xi_true(:) ~= 0);
        metrics.sparsity_match = (metrics.sparsity_pred == metrics.sparsity_true);
        
        % Per-equation coefficient metrics
        for i = 1:size(Xi_true, 2)
            eq_name = ['eq_' num2str(i)];
            coeff_error_i = coeff_error(:, i);
            metrics.(['rmse_coeff_' eq_name]) = sqrt(mean(coeff_error_i.^2));
            metrics.(['mae_coeff_' eq_name]) = mean(abs(coeff_error_i));
        end
    else
        metrics.rmse_coeff = NaN;
        metrics.mae_coeff = NaN;
        metrics.coeff_error_l2 = NaN;
        metrics.coeff_identification_rate = NaN;
        metrics.sparsity_pred = NaN;
        metrics.sparsity_true = NaN;
        metrics.sparsity_match = NaN;
    end
end

