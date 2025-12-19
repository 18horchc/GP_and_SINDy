function [gprMdl, Xnew, ypred, dydt] = fitAndPlotGP(X, y, titleStr, ylabelStr, standardizeToggle)
    % 1. Data Transformation
    y_sqrt = sqrt(y + 1e-6); 

    % 2. Fit the Gaussian Process Model
    gprMdl = fitrgp(X, y_sqrt, ...
        'KernelFunction', 'squaredexponential', ...
        'BasisFunction', 'constant', ...
        'FitMethod', 'exact', ...
        'PredictMethod', 'exact', ...
        'Standardize', standardizeToggle);

    % 3. Generate Predictions for the Mean
    Xnew = linspace(min(X), max(X), 300)';
    [ypred_sqrt, ~, yint_sqrt] = predict(gprMdl, Xnew);

    % --- 4. ANALYTICAL DERIVATIVE CALCULATION ---
    % Extract hyperparameters: SigmaL (Length scale) and SigmaF (Signal std dev)
    L = gprMdl.KernelInformation.KernelParameters(1); 
    sigma_f = gprMdl.KernelInformation.KernelParameters(2);
    
    % Get alpha coefficients: alpha = (K + sigma_n^2*I)^-1 * y
    alpha = gprMdl.Alpha;
    X_train = gprMdl.X; % Training inputs (standardized if toggle was true)

    % Calculate the derivative of the Mean function for Squared Exponential Kernel
    % dmu/dx = sum_{i=1}^n alpha_i * dK(x, x_i)/dx
    num_new = length(Xnew);
    dmu_dx_sqrt = zeros(num_new, 1);
    
    for i = 1:num_new
        % SE Kernel Derivative: K(x,x') * (-(x - x') / L^2)
        dist = Xnew(i) - X_train;
        K_val = (sigma_f^2) * exp(-0.5 * (dist.^2) / (L^2));
        dK_dx = K_val .* (-dist / L^2);
        dmu_dx_sqrt(i) = sum(alpha .* dK_dx);
    end

    % --- 5. CHAIN RULE FOR ORIGINAL SCALE ---
    % Since y = y_sqrt^2, dy/dt = 2 * y_sqrt * d(y_sqrt)/dt
    ypred = max(0, ypred_sqrt).^2;
    dydt = 2 * ypred_sqrt .* dmu_dx_sqrt; 
    
    % Transform intervals for plotting
    ylow  = max(0, yint_sqrt(:,1)).^2;
    yhigh = max(0, yint_sqrt(:,2)).^2;

    % 6. Plotting
    hold on;
    scatter(X, y, 50, 'r', 'filled', 'DisplayName', 'Observed Data');
    plot(Xnew, ypred, 'b-', 'LineWidth', 2, 'DisplayName', 'GP Mean ($\hat{y}$)');
    plot(Xnew, ylow,  'k--', 'LineWidth', 1.4, 'DisplayName', '95% Prediction Interval');
    plot(Xnew, yhigh, 'k--', 'LineWidth', 1.4, 'HandleVisibility', 'off');
    
    xlabel('Time (days)');
    ylabel(ylabelStr);
    title(titleStr);
    legend('Location', 'best', 'Interpreter', 'latex');
    grid on;
end