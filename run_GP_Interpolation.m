function [t_fine, X_gp, dXdt_gp] = run_GP_Interpolation(t_meas, X_meas, resample_factor)
    % RUN_GP_INTERPOLATION Fits a GP to noisy data and returns smooth states + derivatives.
    %
    % Inputs:
    %   t_meas - Sparse, noisy time points (N x 1)
    %   X_meas - Sparse, noisy state data (N x D)
    %   resample_factor - How many points to generate (e.g., 100)
    
    [N, D] = size(X_meas);
    t_fine = linspace(min(t_meas), max(t_meas), resample_factor)';
    
    X_gp = zeros(resample_factor, D);
    dXdt_gp = zeros(resample_factor, D);

    % Process each state variable (channel) independently
    for i = 1:D
        y = X_meas(:, i);
        
        % 1. Fit the GP Model
        % We use the Squared Exponential kernel ('matern52' is also a good choice)
        gpModel = fitrgp(t_meas, y, 'KernelFunction', 'squaredexponential', ...
                         'Standardize', true);
        
        % 2. Get the Predictive Mean (Smooth Trajectory)
        X_gp(:, i) = predict(gpModel, t_fine);
        
        % 3. Calculate the Analytical Derivative
        % The GP mean is mu = K(t_fine, t_meas) * alpha
        % Where alpha = (K(t_meas, t_meas) + sigma^2*I)^-1 * y
        alpha = gpModel.Alpha;
        X_train = gpModel.X;
        sigma_f = gpModel.KernelInformation.KernelParameters(1); % Signal Std
        L = gpModel.KernelInformation.KernelParameters(2);       % Length scale
        
        % Compute the derivative of the kernel matrix dK/dt
        for j = 1:resample_factor
            t_star = t_fine(j);
            % Difference between query point and training points
            dist = t_star - X_train; 
            
            % Squared Exponential Kernel value k(t_star, t_train)
            k_star = sigma_f^2 * exp(-0.5 * (dist.^2) / L^2);
            
            % Derivative of SE kernel: dk/dt_star = - (t_star - t_train)/L^2 * k
            dk_dt = -(dist / L^2) .* k_star;
            
            % The derivative of the mean is the sum of (dk_dt * alpha)
            dXdt_gp(j, i) = sum(dk_dt .* alpha);
        end
    end
end


%Benefits:
% Noise Rejection: The fitrgp function automatically estimates the noise 
% level ($\sigma_n$) and ignores high-frequency jitter.

%The GP derivative is a weighted sum of smooth kernel derivatives, meaning
% $\dot{X}$ will be smooth even if the input was very noisy.

%Variable Upsampling: You can set resample_factor to 500 even if you only 
% started with 20 points, giving the SINDy algorithm plenty of "snapshots" 
% to perform its regression.