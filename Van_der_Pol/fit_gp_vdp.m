function [gp_models, X_gp, dXdt_gp, gp_stats] = fit_gp_vdp(t_sparse, X_sparse, t_dense, options)
%FIT_GP_VDP Fit Gaussian Process models to Van der Pol sparse noisy data
%
% Inputs:
%   t_sparse: Time points where sparse noisy data is available [N x 1]
%   X_sparse: Noisy observations [N x 2] where columns are [x, y]
%   t_dense: Dense time grid for GP predictions [M x 1]
%   options: (optional) Structure with GP fitting options:
%     - kernel: Kernel function (default: 'squaredexponential')
%     - standardize: Whether to standardize (default: false)
%     - sigma: Initial noise std (default: auto)
%     - basis: Basis function (default: 'constant')
%
% Outputs:
%   gp_models: Cell array of GP models {gprMdl_x, gprMdl_y}
%   X_gp: GP predictions on dense grid [M x 2] where columns are [x_gp, y_gp]
%   dXdt_gp: GP analytical derivatives [M x 2] where columns are [dx/dt, dy/dt]
%   gp_stats: Structure with statistics (R², hyperparameters, etc.)

    if nargin < 4
        options = struct();
    end
    
    % Set default options
    if ~isfield(options, 'kernel')
        options.kernel = 'squaredexponential';
    end
    if ~isfield(options, 'standardize')
        options.standardize = false;  % Keep false for clean derivative math
    end
    if ~isfield(options, 'basis')
        options.basis = 'constant';
    end
    
    % Ensure column vectors
    t_sparse = t_sparse(:);
    t_dense = t_dense(:);
    
    num_vars = size(X_sparse, 2);
    num_dense = length(t_dense);
    
    % Initialize outputs
    gp_models = cell(num_vars, 1);
    X_gp = zeros(num_dense, num_vars);
    dXdt_gp = zeros(num_dense, num_vars);
    gp_stats = struct();
    gp_stats.r2 = zeros(num_vars, 1);
    gp_stats.hyperparams = cell(num_vars, 1);
    
    % Fit GP for each state variable
    for i = 1:num_vars
        y_sparse = X_sparse(:, i);
        
        % Fit GP model
        if isfield(options, 'sigma')
            gpModel = fitrgp(t_sparse, y_sparse, ...
                'KernelFunction', options.kernel, ...
                'BasisFunction', options.basis, ...
                'Standardize', options.standardize, ...
                'Sigma', options.sigma);
        else
            gpModel = fitrgp(t_sparse, y_sparse, ...
                'KernelFunction', options.kernel, ...
                'BasisFunction', options.basis, ...
                'Standardize', options.standardize);
        end
        
        gp_models{i} = gpModel;
        
        % Predict mean on dense grid
        X_gp(:, i) = predict(gpModel, t_dense);
        
        % Compute analytical derivatives
        % Extract GP parameters
        alpha = gpModel.Alpha;
        X_train = gpModel.X;
        
        % Get kernel parameters
        % Note: MATLAB's fitrgp returns parameters in order [lengthscale, signal_variance]
        kernel_params = gpModel.KernelInformation.KernelParameters;
        if length(kernel_params) >= 2
            L = kernel_params(1);  % Lengthscale
            sigma_f = kernel_params(2);  % Signal standard deviation
        else
            % Fallback if only one parameter
            L = kernel_params(1);
            sigma_f = 1.0;  % Default
        end
        
        % Compute derivative for each dense time point
        for j = 1:num_dense
            t_star = t_dense(j);
            dist = t_star - X_train;
            
            % Squared exponential kernel: k(t, t') = sigma_f^2 * exp(-0.5 * (t-t')^2 / L^2)
            k_star = sigma_f^2 * exp(-0.5 * (dist.^2) / (L^2));
            
            % Derivative of kernel with respect to t_star:
            % dk/dt = -(t_star - t_train) / L^2 * k(t_star, t_train)
            dk_dt = -(dist / (L^2)) .* k_star;
            
            % Derivative of GP mean: dμ/dt = sum(alpha_i * dk/dt)
            dXdt_gp(j, i) = dk_dt' * alpha;
        end
        
        % Compute R² for validation
        y_pred_train = predict(gpModel, t_sparse);
        ss_res = sum((y_sparse - y_pred_train).^2);
        ss_tot = sum((y_sparse - mean(y_sparse)).^2);
        if ss_tot > 0
            gp_stats.r2(i) = 1 - ss_res / ss_tot;
        else
            gp_stats.r2(i) = NaN;  % Constant data
        end
        
        % Store hyperparameters
        gp_stats.hyperparams{i} = struct();
        gp_stats.hyperparams{i}.lengthscale = L;
        gp_stats.hyperparams{i}.signal_std = sigma_f;
        gp_stats.hyperparams{i}.noise_std = gpModel.Sigma;
    end
    
    % Store variable names
    gp_stats.var_names = {'x', 'y'};
end

