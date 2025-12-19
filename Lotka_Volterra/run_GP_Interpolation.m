function [t_fine, X_gp, dXdt_gp] = run_GP_Interpolation(t_meas, X_meas, resample_factor)
    t_meas = t_meas(:);
    [~, D] = size(X_meas);
    t_fine = linspace(min(t_meas), max(t_meas), resample_factor)';
    
    X_gp = zeros(resample_factor, D);
    dXdt_gp = zeros(resample_factor, D);

    for i = 1:D
        y = X_meas(:, i);
        
        % Fit the GP Model
        % We disable Standardization for a moment to keep the derivative math clean 
        % unless your time scale is massive (e.g., t > 10^5).
        gpModel = fitrgp(t_meas, y, 'KernelFunction', 'squaredexponential', ...
                         'Standardize', false);
        
        % 1. Predict Smooth Trajectory
        X_gp(:, i) = predict(gpModel, t_fine);
        
        % 2. Extract Parameters for Analytical Derivative
        alpha = gpModel.Alpha;
        X_train = gpModel.X; 
        sigma_f = gpModel.KernelInformation.KernelParameters(1); 
        L = gpModel.KernelInformation.KernelParameters(2);       
        
        % 3. Analytical Derivative Calculation
        % The derivative of the Mean is: d_mu/dt = sum( dK/dt * alpha )
        for j = 1:resample_factor
            t_star = t_fine(j);
            dist = t_star - X_train; 
            
            % Kernel k(t_star, t_train)
            k_star = sigma_f^2 * exp(-0.5 * (dist.^2) / L^2);
            
            % Derivative of Kernel dk/dt_star
            dk_dt = -(dist / L^2) .* k_star;
            
            % Compute the weighted sum
            dXdt_gp(j, i) = dk_dt' * alpha;
        end
    end
end