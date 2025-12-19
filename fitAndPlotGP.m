function [Xnew_time, Y_gp, dY_dt] = fitAndPlotGP(t, X_data, n_points)
    % Multi-variable GP fitting for SINDy
    num_vars = size(X_data, 2);
    Xnew_time = linspace(min(t), max(t), n_points)';
    
    Y_gp = zeros(n_points, num_vars);
    dY_dt = zeros(n_points, num_vars);
    
    for i = 1:num_vars
        y = X_data(:,i);
        % Note: Using your sqrt transform logic from earlier
        y_sqrt_obs = sqrt(y + 1e-6); 

        gprMdl = fitrgp(t, y_sqrt_obs, 'KernelFunction', 'squaredexponential', ...
            'BasisFunction', 'constant', 'Standardize', true);

        % Predictions
        [ypred_sqrt, ~, ~] = predict(gprMdl, Xnew_time);
        
        % Analytical Derivative of the Mean (as established earlier)
        L = gprMdl.KernelInformation.KernelParameters(1); 
        sigma_f = gprMdl.KernelInformation.KernelParameters(2);
        alpha = gprMdl.Alpha;
        X_train = gprMdl.X;
        
        dmu_dx_sqrt = zeros(n_points, 1);
        for j = 1:n_points
            dist = Xnew_time(j) - X_train;
            K_val = (sigma_f^2) * exp(-0.5 * (dist.^2) / (L^2));
            dK_dx = K_val .* (-dist / L^2);
            dmu_dx_sqrt(j) = sum(alpha .* dK_dx);
        end

        % Transform back & Apply Chain Rule
        Y_gp(:,i) = max(0, ypred_sqrt).^2;
        dY_dt(:,i) = 2 * ypred_sqrt .* dmu_dx_sqrt; 
    end
end