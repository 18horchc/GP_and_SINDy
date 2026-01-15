% Comparison of SINDy on sparse noisy data (Path 1) vs GP-enhanced SINDy (Path 2)
% For Lotka-Volterra equations
%
% Path 1: Direct SINDy on 10 noisy sparse points (averaged from 30 measurements)
% Path 2: Fit GP to 30 points -> sample more data -> use GP derivative -> SINDy with expanded library

clear; close all; clc

%% ========== Generate Lotka-Volterra Data ==========
% Lotka-Volterra equations:
% dx/dt = alpha*x - beta*x*y
% dy/dt = delta*x*y - gamma*y

% Parameters for Lotka-18
alpha = 1.5;  % Prey growth rate
beta = 1.0;   % Predation rate
gamma = 3.0;  % Predator death rate
delta = 1.0;  % Predator growth rate

% Initial conditions
x0 = 1.0;  % Prey
y0 = 1.0;  % Predator

% Time span for generating true solution
t_span = [0, 10];
t_true = linspace(0, 10, 1000);

% Solve Lotka-Volterra ODE
[t_sol, sol] = ode45(@(t, y) lotka_volterra_ode(t, y, alpha, beta, gamma, delta), ...
                     t_span, [x0; y0], odeset('RelTol', 1e-8));

% Interpolate to get values at specific time points
x_true = interp1(t_sol, sol(:,1), t_true);
y_true = interp1(t_sol, sol(:,2), t_true);

% Select 10 random time points
num_unique_points = 10;
t_sparse_unique = sort(rand(1, num_unique_points) * 10);  % 10 random time points in [0, 10]
x_sparse_true_unique = interp1(t_sol, sol(:,1), t_sparse_unique);
y_sparse_true_unique = interp1(t_sol, sol(:,2), t_sparse_unique);

% Generate 3 independent noisy measurements at each time point
% This gives us 30 total data points (3 at each of 10 time points)
num_replicates = 3;  % Number of measurements per time point
noise_level = 0.1;  % 10% noise

% Initialize arrays for all 30 points
t_sparse = repmat(t_sparse_unique, num_replicates, 1);
t_sparse = t_sparse(:);  % Flatten to column vector
x_sparse_noisy = zeros(size(t_sparse));
y_sparse_noisy = zeros(size(t_sparse));

% Generate independent noisy measurements at each time point
for i = 1:length(t_sparse_unique)
    % True values at this time point
    x_true_i = x_sparse_true_unique(i);
    y_true_i = y_sparse_true_unique(i);
    
    % Generate 3 independent noisy measurements
    % 10% noise means noise std = 10% of signal value at each point
    for j = 1:num_replicates
        idx = (i-1)*num_replicates + j;
        noise_x = normrnd(0, noise_level * abs(x_true_i));
        noise_y = normrnd(0, noise_level * abs(y_true_i));
        
        x_sparse_noisy(idx) = x_true_i + noise_x;
        y_sparse_noisy(idx) = y_true_i + noise_y;
    end
end

% Ensure non-negative (Lotka-Volterra should be positive)
x_sparse_noisy = max(0.01, x_sparse_noisy);
y_sparse_noisy = max(0.01, y_sparse_noisy);

fprintf('Generated %d sparse noisy points from Lotka-Volterra\n', length(t_sparse));
fprintf('  - %d unique random time points\n', length(t_sparse_unique));
fprintf('  - %d replicates per time point\n', num_replicates);
fprintf('  - Total: %d data points\n', length(t_sparse));

%% ========== PATH 1: Direct SINDy on averaged data ==========
fprintf('\n========== PATH 1: Direct SINDy ==========\n');

% For Path 1, we use AVERAGED data at each unique time point
% This gives us 10 points (one average per time point) instead of 30 raw points
% Average replicates at each unique time point to get 10 sequential points
xdata_path1 = zeros(size(t_sparse_unique));
ydata_path1 = zeros(size(t_sparse_unique));

for i = 1:length(t_sparse_unique)
    idx_at_time = abs(t_sparse - t_sparse_unique(i)) < 1e-6;
    xdata_path1(i) = mean(x_sparse_noisy(idx_at_time));
    ydata_path1(i) = mean(y_sparse_noisy(idx_at_time));
end

xdata_path1 = xdata_path1';
ydata_path1 = ydata_path1';

fprintf('Path 1 uses averaged data: %d points (averages at %d unique random time points)\n', ...
        length(xdata_path1), length(t_sparse_unique));

% SINDy uses continuous-time formulation: dX/dt = f(X)
% Compute derivatives using finite differences
% Sort data by time to ensure proper derivative computation
[t_sorted, sort_idx] = sort(t_sparse_unique);
xdata_sorted = xdata_path1(sort_idx);
ydata_sorted = ydata_path1(sort_idx);

% Compute derivatives using finite differences
% Use central differences for interior points, forward/backward for endpoints
dxdt_path1 = zeros(size(xdata_sorted));
dydt_path1 = zeros(size(ydata_sorted));

for i = 1:length(t_sorted)
    if i == 1
        % Forward difference at first point
        dt = t_sorted(2) - t_sorted(1);
        dxdt_path1(i) = (xdata_sorted(2) - xdata_sorted(1)) / dt;
        dydt_path1(i) = (ydata_sorted(2) - ydata_sorted(1)) / dt;
    elseif i == length(t_sorted)
        % Backward difference at last point
        dt = t_sorted(end) - t_sorted(end-1);
        dxdt_path1(i) = (xdata_sorted(end) - xdata_sorted(end-1)) / dt;
        dydt_path1(i) = (ydata_sorted(end) - ydata_sorted(end-1)) / dt;
    else
        % Central difference for interior points (more accurate)
        dt = t_sorted(i+1) - t_sorted(i-1);
        dxdt_path1(i) = (xdata_sorted(i+1) - xdata_sorted(i-1)) / dt;
        dydt_path1(i) = (ydata_sorted(i+1) - ydata_sorted(i-1)) / dt;
    end
end

% Library functions evaluated at all points
% Theta_path1: Library functions evaluated at X(t) [current state]
% X_dot_path1: dX/dt [derivative] - this is what we're trying to predict
Theta_path1 = [ones(length(xdata_sorted), 1) ...
               xdata_sorted(:) ...
               ydata_sorted(:) ...
               xdata_sorted(:).^2 ...
               ydata_sorted(:).^2 ...
               xdata_sorted(:).*ydata_sorted(:)];

X_dot_path1 = [dxdt_path1(:) dydt_path1(:)];

fprintf('Computed derivatives using finite differences (central/forward/backward)\n');

% Standard SINDy: Sequential Thresholded Least Squares (STLS)
% 1) Initial least-squares fit: Xi = Theta \ X_dot
% 2) Threshold small coefficients and refit on remaining terms
lambda_base = 0.1;
lambda = lambda_base * std(X_dot_path1(:));  % Threshold scaled by data magnitude
max_iter = 10;
epsguess_path1 = sindy_stls(Theta_path1, X_dot_path1, lambda, max_iter);

fprintf('Path 1 Model Coefficients (SINDy-STLS):\n');
fprintf('  X: [%s]\n', num2str(epsguess_path1(:,1)'));
fprintf('  Y: [%s]\n', num2str(epsguess_path1(:,2)'));

%% ========== PATH 2: GP-enhanced SINDy ==========
fprintf('\n========== PATH 2: GP-enhanced SINDy ==========\n');

% Fit GP to sparse noisy data
X_gp = t_sparse(:);
y1_gp = x_sparse_noisy(:);
y2_gp = y_sparse_noisy(:);

% Fit GP models - try different kernels to get better fit
% 
% Options for KernelFunction:
% - 'squaredexponential': Smooth, infinitely differentiable (default)
% - 'matern32': Less smooth, allows for more variation
% - 'matern52': Intermediate smoothness
% - 'rationalquadratic': More flexible, can capture multiple scales
% - 'ardsquaredexponential': Automatic relevance determination (ARD)
%
% For Lotka-Volterra with oscillatory behavior, Matern or rational quadratic
% might work better than squared exponential.

% Try Rational Quadratic kernel for better fit
% Rational Quadratic is more flexible than Matern32 and can capture
% multiple length scales, which is good for oscillatory data
% 
% NOTE: For composite kernels, we'd need custom kernel functions, which
% makes parameter extraction for derivatives more complex. Starting with
% Rational Quadratic as it's built-in and more flexible than Matern32.

% Using Squared Exponential kernel
% Note: Smooth kernel with well-defined analytical derivatives
% The paper uses constant basis functions for simplicity
fprintf('Fitting GP models with Squared Exponential kernel...\n');
gprMdl_x = fitrgp(X_gp, y1_gp, ...
    'KernelFunction', 'squaredexponential', ...
    'BasisFunction', 'constant', ...  % Constant basis (as in Wang & Barber)
    'FitMethod', 'exact', ...
    'PredictMethod', 'exact', ...
    'Standardize', false);

gprMdl_y = fitrgp(X_gp, y2_gp, ...
    'KernelFunction', 'squaredexponential', ...
    'BasisFunction', 'constant', ...  % Constant basis (as in Wang & Barber)
    'FitMethod', 'exact', ...
    'PredictMethod', 'exact', ...
    'Standardize', false);

fprintf('GP kernel used: %s\n', gprMdl_x.KernelInformation.Name);
fprintf('GP basis function: %s\n', gprMdl_x.BasisFunction);
r2_x = compute_gp_r2(gprMdl_x, X_gp, y1_gp);
r2_y = compute_gp_r2(gprMdl_y, X_gp, y2_gp);
fprintf('GP fit quality - X: R² = %.4f, Y: R² = %.4f\n', r2_x, r2_y);

if r2_x < 0.5 || r2_y < 0.5
    warning('GP fit quality is poor (R² < 0.5). Consider adjusting hyperparameters or kernel.');
end

% Generate dense time points for posterior sampling
num_aug_points = 1000;
t_dense = linspace(0, 10, num_aug_points);
X_dense = t_dense(:);

% Predict from GP (directly, no transform needed)
[y1pred, y1pred_std, ~] = predict(gprMdl_x, X_dense);
[y2pred, y2pred_std, ~] = predict(gprMdl_y, X_dense);

% Check GP predictions
fprintf('GP prediction diagnostics:\n');
fprintf('  X mean: range=[%.4f, %.4f], std=%.4f\n', min(y1pred), max(y1pred), std(y1pred));
fprintf('  Y mean: range=[%.4f, %.4f], std=%.4f\n', min(y2pred), max(y2pred), std(y2pred));

% Ensure non-negative
y1pred = max(0.01, y1pred);
y2pred = max(0.01, y2pred);

% Sample ONE GP curve from the posterior distribution
% This yields 1000 points along a sampled trajectory
num_gp_samples = 1;
y1_gp_samples = sample_gp_posterior(gprMdl_x, X_dense, num_gp_samples);
y2_gp_samples = sample_gp_posterior(gprMdl_y, X_dense, num_gp_samples);

% Use the sampled curve
y1_sampled = y1_gp_samples(:, 1);
y2_sampled = y2_gp_samples(:, 1);

% Ensure non-negative
y1_sampled = max(0.01, y1_sampled);
y2_sampled = max(0.01, y2_sampled);

% Use the sampled curve as the augmented dataset
t_combined = t_dense;
x_combined = y1_sampled;
y_combined = y2_sampled;

% Weight points closer to the GP mean higher
eps_std = 1e-6;
z_x = (x_combined - y1pred) ./ (y1pred_std + eps_std);
z_y = (y_combined - y2pred) ./ (y2pred_std + eps_std);
weights = exp(-0.5 * (z_x.^2 + z_y.^2));
weights = weights / mean(weights);  % Normalize for numerical stability

fprintf('Augmented dataset: %d GP posterior sample points\n', length(t_combined));
fprintf('  Weight stats: min=%.4f, max=%.4f, mean=%.4f\n', min(weights), max(weights), mean(weights));

% Compute GP derivative analytically for ALL points
% Following Hsin et al. (2025): use GP analytical derivatives for all augmented data points
% The derivative of a GP is itself a GP! We can compute it by differentiating
% the covariance kernel directly, which is more principled than finite differences.
%
% Mathematical foundation (example for squared exponential kernel):
% - The derivative GP has covariance: k'(x, x') = σ²/l² * (1 - (x-x')²/l²) * exp(-(x-x')²/(2l²))
% - The mean of the derivative GP: μ' = K'(X_new, X) * K(X,X)^(-1) * y
% Note: Squared Exponential kernel is smooth; analytical derivatives are well-defined.
%
% How this fits into SINDy workflow:
% 1. We fit GP to sparse noisy data → get smooth function estimate
% 2. We compute analytical derivative of GP at ALL points (sparse + dense) → get smooth derivative estimate
% 3. We use these derivatives in continuous-time SINDy: dX/dt = f(X)
% 4. The GP derivative provides denoised, uncertainty-quantified derivatives
%    that are more reliable than finite differences, especially with sparse data
%
% Note: t_combined = t_dense, so we compute derivatives for all points in the combined dataset
% This ensures consistent GP-derived derivatives throughout, matching Hsin et al. approach
[y1_deriv_gp, y1_deriv_var] = gp_derivative_analytical(gprMdl_x, X_dense);
[y2_deriv_gp, y2_deriv_var] = gp_derivative_analytical(gprMdl_y, X_dense);

% Check GP derivative quality
fprintf('\nGP Derivative Quality Check:\n');
if any(isnan(y1_deriv_gp)) || any(isnan(y2_deriv_gp))
    warning('GP derivatives contain NaN values. This may cause poor model discovery.');
end
if any(isinf(y1_deriv_gp)) || any(isinf(y2_deriv_gp))
    warning('GP derivatives contain Inf values. This may cause poor model discovery.');
end
if std(y1_deriv_gp) < 0.01 || std(y2_deriv_gp) < 0.01
    warning('GP derivatives have very low variance. GP fit may be too smooth.');
end

% Check for stability issues and print diagnostics
fprintf('GP derivative diagnostics:\n');
fprintf('  X derivative: mean=%.4f, std=%.4f, min=%.4f, max=%.4f\n', ...
        mean(y1_deriv_gp), std(y1_deriv_gp), min(y1_deriv_gp), max(y1_deriv_gp));
fprintf('  Y derivative: mean=%.4f, std=%.4f, min=%.4f, max=%.4f\n', ...
        mean(y2_deriv_gp), std(y2_deriv_gp), min(y2_deriv_gp), max(y2_deriv_gp));

% Print kernel hyperparameters for debugging
params_x = gprMdl_x.KernelInformation.KernelParameters;
params_y = gprMdl_y.KernelInformation.KernelParameters;
kernel_name = gprMdl_x.KernelInformation.Name;
if strcmpi(kernel_name, 'SquaredExponential')
    % Squared Exponential: [SigmaL, SigmaF] = [l, sigma_f]
    fprintf('  X GP kernel params: l=%.4f, sigma_f=%.4f\n', params_x(1), params_x(2));
    fprintf('  Y GP kernel params: l=%.4f, sigma_f=%.4f\n', params_y(1), params_y(2));
elseif strcmpi(kernel_name, 'RationalQuadratic')
    % Rational Quadratic: [SigmaL, AlphaRQ, SigmaF] = [l, alpha, sigma_f]
    fprintf('  X GP kernel params: l=%.4f, alpha=%.4f, sigma_f=%.4f\n', ...
            params_x(1), params_x(2), params_x(3));
    fprintf('  Y GP kernel params: l=%.4f, alpha=%.4f, sigma_f=%.4f\n', ...
            params_y(1), params_y(2), params_y(3));
else
    % Other kernels - print all parameters
    fprintf('  X GP kernel params: %s\n', mat2str(params_x));
    fprintf('  Y GP kernel params: %s\n', mat2str(params_y));
end

% Check for problematic hyperparameters (for any kernel)
min_length_scale = 0.01;  % Minimum reasonable length scale
if params_x(1) < min_length_scale
    warning('Path 2: X GP length scale (l=%.4f) is too small. This may cause numerical instability.', params_x(1));
    fprintf('  This is likely due to the GP optimization finding a poor local minimum with sparse/noisy data.\n');
end
if params_y(1) < min_length_scale
    warning('Path 2: Y GP length scale (l=%.4f) is too small. This may cause numerical instability.', params_y(1));
end

if any(isnan(y1_deriv_gp)) || any(isinf(y1_deriv_gp))
    warning('Path 2: GP derivative for X contains NaN or Inf values. This may cause issues.');
end
if any(isnan(y2_deriv_gp)) || any(isinf(y2_deriv_gp))
    warning('Path 2: GP derivative for Y contains NaN or Inf values. This may cause issues.');
end
% Check if analytical derivatives are reliable
% If derivative variance is too low, the analytical computation may be incorrect
% (e.g., due to kernel-specific issues) - use finite differences as fallback
y1_deriv_std = std(y1_deriv_gp);
if isnan(y1_deriv_std) || y1_deriv_std < 0.01
    warning('Path 2: GP derivative for X has very low variance (std=%.4f). Using finite differences fallback.', y1_deriv_std);
    % Fallback to finite differences if analytical derivative is too small
    dt = mean(diff(X_dense));
    y1_deriv_gp = gradient(y1pred, dt);
    fprintf('  Using finite differences for X derivative (new std=%.4f)\n', std(y1_deriv_gp));
end
y2_deriv_std = std(y2_deriv_gp);
if isnan(y2_deriv_std) || y2_deriv_std < 0.01 || all(abs(y2_deriv_gp) < 1e-10)
    warning('Path 2: GP derivative for Y has very low variance (std=%.4f). Using finite differences fallback.', y2_deriv_std);
    fprintf('  Y GP fit quality (R²=%.4f) is poor - GP mean may be nearly constant.\n', r2_y);
    dt = mean(diff(X_dense));
    y2_deriv_gp = gradient(y2pred, dt);
    fprintf('  Using finite differences for Y derivative (new std=%.4f)\n', std(y2_deriv_gp));
end

% Use GP analytical derivatives for all points in t_combined
% Since t_combined = t_dense, the derivatives are already computed for all points
% This matches Hsin et al. (2025) approach: GP analytical derivatives for all augmented data
y1_deriv_combined = y1_deriv_gp;
y2_deriv_combined = y2_deriv_gp;

% For SINDy with derivatives, we'll use a continuous-time approach
% similar to SINDy but with the expanded library

% Expanded library (since we have more data, we can include more terms)
% Basic library: [1, x, y, x^2, y^2, x*y]
% Expanded: add [x^3, y^3, x^2*y, x*y^2, sin(x), sin(y), etc.]
Theta_path2 = [ones(length(t_combined), 1) ...
               x_combined(:) ...
               y_combined(:) ...
               x_combined(:).^2 ...
               y_combined(:).^2 ...
               x_combined(:).*y_combined(:) ...
               x_combined(:).^3 ...
               y_combined(:).^3 ...
               x_combined(:).^2.*y_combined(:) ...
               x_combined(:).*y_combined(:).^2];

X_dot_path2 = [y1_deriv_combined(:) y2_deriv_combined(:)];

% Weighted SINDy (STLS) using GP posterior samples
W = sqrt(weights(:));
Theta_path2_w = Theta_path2 .* W;
X_dot_path2_w = X_dot_path2 .* W;

lambda_base = 0.05;
lambda = lambda_base * std(X_dot_path2_w(:));
max_iter = 10;
epsguess_path2 = sindy_stls(Theta_path2_w, X_dot_path2_w, lambda, max_iter);

fprintf('Path 2 Model Coefficients (weighted SINDy, %d functions):\n', size(Theta_path2, 2));
fprintf('  X: [%s]\n', num2str(epsguess_path2(:,1)'));
fprintf('  Y: [%s]\n', num2str(epsguess_path2(:,2)'));

%% ========== Visualization ==========
fprintf('\n========== Generating Plots ==========\n');

% Plot 1: Original data and true solution
figure(1);
subplot(2,1,1);
plot(t_true, x_true, 'b-', 'LineWidth', 2, 'DisplayName', 'True X (prey)');
hold on;
scatter(t_sparse, x_sparse_noisy, 80, 'r', 'filled', 'DisplayName', 'Noisy sparse data (30 points)');
xlabel('Time');
ylabel('X (Prey)');
title('Lotka-Volterra: Prey Population (30 points: 3 replicates at 10 random time points)');
legend('Location', 'best');
grid on;

subplot(2,1,2);
plot(t_true, y_true, 'b-', 'LineWidth', 2, 'DisplayName', 'True Y (predator)');
hold on;
scatter(t_sparse, y_sparse_noisy, 80, 'r', 'filled', 'DisplayName', 'Noisy sparse data (30 points)');
xlabel('Time');
ylabel('Y (Predator)');
title('Lotka-Volterra: Predator Population (30 points: 3 replicates at 10 random time points)');
legend('Location', 'best');
grid on;

% Plot 2: GP fit and sampling
figure(2);
subplot(2,1,1);
plot(t_dense, y1pred, 'b-', 'LineWidth', 2, 'DisplayName', 'GP mean');
hold on;
scatter(t_sparse, x_sparse_noisy, 80, 'r', 'filled', 'DisplayName', 'Original 30 points');
scatter(t_dense, y1_sampled, 50, 'g', 'o', 'DisplayName', 'GP samples');
xlabel('Time');
ylabel('X (Prey)');
title('Path 2: GP Fit and Sampling (Prey)');
legend('Location', 'best');
grid on;

subplot(2,1,2);
plot(t_dense, y2pred, 'b-', 'LineWidth', 2, 'DisplayName', 'GP mean');
hold on;
scatter(t_sparse, y_sparse_noisy, 80, 'r', 'filled', 'DisplayName', 'Original 30 points');
scatter(t_dense, y2_sampled, 50, 'g', 'o', 'DisplayName', 'GP samples');
xlabel('Time');
ylabel('Y (Predator)');
title('Path 2: GP Fit and Sampling (Predator)');
legend('Location', 'best');
grid on;

% Plot 3: Model predictions comparison
figure(3);

% Compute averaged data points for Path 2 (same as Path 1 for consistency)
xdata_path2_avg = zeros(size(t_sparse_unique));
ydata_path2_avg = zeros(size(t_sparse_unique));
for i = 1:length(t_sparse_unique)
    idx_at_time = abs(t_sparse - t_sparse_unique(i)) < 1e-6;
    xdata_path2_avg(i) = mean(x_sparse_noisy(idx_at_time));
    ydata_path2_avg(i) = mean(y_sparse_noisy(idx_at_time));
end

% Simulate Path 1 model (continuous-time ODE)
% The model dX/dt = f(X) was learned using finite differences
% This is a direct continuous-time ODE, similar to Path 2

% Validate model before simulation
if any(isnan(epsguess_path1(:))) || any(isinf(epsguess_path1(:)))
    warning('Path 1: Model contains NaN or Inf coefficients. Using zero model.');
    epsguess_path1(isnan(epsguess_path1) | isinf(epsguess_path1)) = 0;
end

% Clip extreme coefficients to prevent instability
% Large coefficients can cause rapid divergence
max_coeff = 50;  % Maximum reasonable coefficient magnitude
num_clipped = sum(abs(epsguess_path1(:)) > max_coeff);
if num_clipped > 0
    warning('Path 1: Clipping %d coefficients with magnitude > %.1f to prevent instability.', num_clipped, max_coeff);
    epsguess_path1(abs(epsguess_path1) > max_coeff) = sign(epsguess_path1(abs(epsguess_path1) > max_coeff)) * max_coeff;
end

% Continuous-time ODE: dX/dt = f(X)
function_path1 = @(t, y) sindy_rhs_path1(t, y, epsguess_path1);

% Use true initial conditions
x0_path1 = x0;
y0_path1 = y0;

% Integrate with ODE solver (more stable than discrete stepping)
try
    % Try with standard tolerances first
    try
        [t_sim_path1, sol_path1] = ode45(function_path1, [0, 10], ...
                                         [x0_path1; y0_path1], ...
                                         odeset('RelTol', 1e-6, 'AbsTol', 1e-6, ...
                                                'MaxStep', 0.5, 'Refine', 2));
    catch
        % If ode45 fails, try ode15s (stiff solver)
        warning('Path 1: ode45 failed. Trying ode15s (stiff solver).');
        [t_sim_path1, sol_path1] = ode15s(function_path1, [0, 10], ...
                                          [x0_path1; y0_path1], ...
                                          odeset('RelTol', 1e-4, 'AbsTol', 1e-4, ...
                                                 'MaxStep', 0.5));
    end
    
    x_sim_path1 = sol_path1(:, 1);
    y_sim_path1 = sol_path1(:, 2);
    
    % Check for numerical issues
    if any(isnan(sol_path1(:))) || any(isinf(sol_path1(:)))
        warning('Path 1: ODE integration produced NaN or Inf. Model may be unstable.');
        % Truncate to valid portion
        valid_idx = ~(any(isnan(sol_path1), 2) | any(isinf(sol_path1), 2));
        if sum(valid_idx) > 1
            t_sim_path1 = t_sim_path1(valid_idx);
            x_sim_path1 = x_sim_path1(valid_idx);
            y_sim_path1 = y_sim_path1(valid_idx);
        else
            % If all invalid, create minimal output
            t_sim_path1 = [0; 10];
            x_sim_path1 = [x0_path1; x0_path1];
            y_sim_path1 = [y0_path1; y0_path1];
        end
    end
    
    % Check if solution stays in reasonable bounds
    if max(abs([x_sim_path1; y_sim_path1])) > 1e6
        warning('Path 1: ODE solution has very large values. Truncating to reasonable portion.');
        reasonable_idx = all(abs(sol_path1) < 1e6, 2);
        if sum(reasonable_idx) > 1
            t_sim_path1 = t_sim_path1(reasonable_idx);
            x_sim_path1 = x_sim_path1(reasonable_idx);
            y_sim_path1 = y_sim_path1(reasonable_idx);
        end
    end
    
catch ME
    % Fix MATLAB warning syntax - use format specifier
    warning('Path 1: ODE integration failed. Using initial condition only.');
    fprintf('  Error details: %s\n', ME.message);
    t_sim_path1 = [0; 10];
    x_sim_path1 = [x0_path1; x0_path1];
    y_sim_path1 = [y0_path1; y0_path1];
end

% Simulate Path 2 model (continuous-time ODE)
function_path2 = @(t, y) sindy_rhs_path2(t, y, epsguess_path2);

% Use true initial conditions for better comparison (or use mean of replicates)
% Option 1: Use true initial conditions
x0_path2 = x0;  % True initial condition
y0_path2 = y0;  % True initial condition

% Option 2: Use mean of replicates (uncomment to use)
% idx_first = abs(t_sparse - t_sparse_unique(1)) < 1e-6;
% x0_path2 = mean(x_sparse_noisy(idx_first));
% y0_path2 = mean(y_sparse_noisy(idx_first));

% Validate model before integration: check for NaN/Inf coefficients
if any(isnan(epsguess_path2(:))) || any(isinf(epsguess_path2(:)))
    warning('Path 2: Model contains NaN or Inf coefficients. Using zero model.');
    epsguess_path2(isnan(epsguess_path2) | isinf(epsguess_path2)) = 0;
end

% Clip extreme coefficients to prevent instability
% Large coefficients (especially in cubic terms) can cause rapid divergence
max_coeff = 50;  % Maximum reasonable coefficient magnitude (reduced from 100 for better stability)
num_clipped = sum(abs(epsguess_path2(:)) > max_coeff);
if num_clipped > 0
    warning('Path 2: Clipping %d coefficients with magnitude > %.1f to prevent instability.', num_clipped, max_coeff);
    epsguess_path2(abs(epsguess_path2) > max_coeff) = sign(epsguess_path2(abs(epsguess_path2) > max_coeff)) * max_coeff;
end

% Additional check: if cubic terms have large coefficients, they're likely problematic
cubic_indices = [7, 8, 9, 10];  % x^3, y^3, x^2*y, x*y^2 in expanded library
if any(abs(epsguess_path2(cubic_indices, :)) > 10)
    warning('Path 2: Large coefficients in cubic terms detected. These may cause instability.');
    % Further reduce cubic term coefficients
    for idx = cubic_indices
        if abs(epsguess_path2(idx, 1)) > 10
            epsguess_path2(idx, 1) = sign(epsguess_path2(idx, 1)) * min(10, abs(epsguess_path2(idx, 1)));
        end
        if abs(epsguess_path2(idx, 2)) > 10
            epsguess_path2(idx, 2) = sign(epsguess_path2(idx, 2)) * min(10, abs(epsguess_path2(idx, 2)));
        end
    end
end

% Integrate with error handling and reasonable tolerances
% Use looser tolerances to avoid step size issues, but still accurate enough
try
    % Try with standard tolerances first (ode45 for non-stiff systems)
    try
        [t_ode_path2, sol_path2] = ode45(function_path2, [0, 10], ...
                                         [x0_path2; y0_path2], ...
                                         odeset('RelTol', 1e-6, 'AbsTol', 1e-6, ...
                                                'MaxStep', 0.5, 'Refine', 2));
    catch
        % If ode45 fails (e.g., step size too small), try ode15s (stiff solver)
        warning('Path 2: ode45 failed. Trying ode15s (stiff solver).');
        [t_ode_path2, sol_path2] = ode15s(function_path2, [0, 10], ...
                                          [x0_path2; y0_path2], ...
                                          odeset('RelTol', 1e-4, 'AbsTol', 1e-4, ...
                                                 'MaxStep', 0.5));
    end
    
    % Check for numerical issues
    if any(isnan(sol_path2(:))) || any(isinf(sol_path2(:)))
        warning('Path 2: ODE integration produced NaN or Inf. Model may be unstable.');
        % Truncate to valid portion
        valid_idx = ~(any(isnan(sol_path2), 2) | any(isinf(sol_path2), 2));
        if sum(valid_idx) > 1
            t_ode_path2 = t_ode_path2(valid_idx);
            sol_path2 = sol_path2(valid_idx, :);
        else
            % If all invalid, create minimal output
            t_ode_path2 = [0; 10];
            sol_path2 = [x0_path2, y0_path2; x0_path2, y0_path2];
        end
    end
    
    % Check if solution stays in reasonable bounds
    if max(abs(sol_path2(:))) > 1e6
        warning('Path 2: ODE solution has very large values. Truncating to reasonable portion.');
        % Truncate to reasonable portion
        reasonable_idx = all(abs(sol_path2) < 1e6, 2);
        if sum(reasonable_idx) > 1
            t_ode_path2 = t_ode_path2(reasonable_idx);
            sol_path2 = sol_path2(reasonable_idx, :);
        end
    end
    
catch ME
    % Fix MATLAB warning syntax - use format specifier
    warning('Path 2: ODE integration failed. Using initial condition only.');
    fprintf('  Error details: %s\n', ME.message);
    t_ode_path2 = [0; 10];
    sol_path2 = [x0_path2, y0_path2; x0_path2, y0_path2];
end

% Path 1: SINDy predictions
subplot(2,2,1);
% Plot true curve first to establish axis limits
plot(t_true, x_true, 'k-', 'LineWidth', 2.5, 'DisplayName', 'True');
hold on;
% Then plot data points
scatter(t_sparse_unique, xdata_path1, 120, 'b', 'filled', 'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.5, 'DisplayName', 'Averaged data (10 pts)');
% Finally plot prediction
plot(t_sim_path1, x_sim_path1, 'r--', 'LineWidth', 2, 'DisplayName', 'Path 1: SINDy');
xlabel('Time');
ylabel('X (Prey)');
title('Path 1: SINDy Prediction');
legend('Location', 'best');
grid on;
% Ensure axis includes all data
axis tight;

subplot(2,2,2);
% Plot true curve first to establish axis limits
plot(t_true, y_true, 'k-', 'LineWidth', 2.5, 'DisplayName', 'True');
hold on;
% Then plot data points
scatter(t_sparse_unique, ydata_path1, 120, 'b', 'filled', 'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.5, 'DisplayName', 'Averaged data (10 pts)');
% Finally plot prediction
plot(t_sim_path1, y_sim_path1, 'r--', 'LineWidth', 2, 'DisplayName', 'Path 1: SINDy');
xlabel('Time');
ylabel('Y (Predator)');
title('Path 1: SINDy Prediction');
legend('Location', 'best');
grid on;
% Ensure axis includes all data
axis tight;

% Path 2: GP-enhanced SINDy predictions
subplot(2,2,3);
% Plot true curve first to establish axis limits
plot(t_true, x_true, 'k-', 'LineWidth', 2.5, 'DisplayName', 'True');
hold on;
% Then plot data points
scatter(t_sparse_unique, xdata_path2_avg, 120, 'b', 'filled', 'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.5, 'DisplayName', 'Averaged data (10 pts)');
% Finally plot prediction
plot(t_ode_path2, sol_path2(:,1), 'g--', 'LineWidth', 2, 'DisplayName', 'Path 2: GP+SINDy');
xlabel('Time');
ylabel('X (Prey)');
title('Path 2: GP-enhanced SINDy Prediction');
legend('Location', 'best');
grid on;
% Ensure axis includes all data
axis tight;

subplot(2,2,4);
% Plot true curve first to establish axis limits
plot(t_true, y_true, 'k-', 'LineWidth', 2.5, 'DisplayName', 'True');
hold on;
% Then plot data points
scatter(t_sparse_unique, ydata_path2_avg, 120, 'b', 'filled', 'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.5, 'DisplayName', 'Averaged data (10 pts)');
% Finally plot prediction
plot(t_ode_path2, sol_path2(:,2), 'g--', 'LineWidth', 2, 'DisplayName', 'Path 2: GP+SINDy');
xlabel('Time');
ylabel('Y (Predator)');
title('Path 2: GP-enhanced SINDy Prediction');
legend('Location', 'best');
grid on;
% Ensure axis includes all data
axis tight;

%% ========== Model Comparison and Error Analysis ==========
fprintf('\n========== Model Comparison ==========\n');

% Interpolate true solution to comparison time points
% Handle case where t_sim_path1 might have NaN or invalid values
valid_idx = ~isnan(t_sim_path1) & ~isnan(x_sim_path1) & ~isnan(y_sim_path1);
if sum(valid_idx) > 0
    x_true_path1 = interp1(t_sol, sol(:,1), t_sim_path1(valid_idx), 'linear', 'extrap');
    y_true_path1 = interp1(t_sol, sol(:,2), t_sim_path1(valid_idx), 'linear', 'extrap');
    
    % Only compute RMSE for valid points
    x_sim_valid = x_sim_path1(valid_idx);
    y_sim_valid = y_sim_path1(valid_idx);
else
    x_true_path1 = [];
    y_true_path1 = [];
    x_sim_valid = [];
    y_sim_valid = [];
end

x_true_path2 = interp1(t_sol, sol(:,1), t_ode_path2, 'linear', 'extrap');
y_true_path2 = interp1(t_sol, sol(:,2), t_ode_path2, 'linear', 'extrap');

% Compute RMSE for each path
if ~isempty(x_sim_valid) && length(x_sim_valid) == length(x_true_path1)
    rmse_path1_x = sqrt(mean((x_sim_valid - x_true_path1).^2));
    rmse_path1_y = sqrt(mean((y_sim_valid - y_true_path1).^2));
else
    rmse_path1_x = NaN;
    rmse_path1_y = NaN;
end
rmse_path2_x = sqrt(mean((sol_path2(:,1) - x_true_path2).^2));
rmse_path2_y = sqrt(mean((sol_path2(:,2) - y_true_path2).^2));

fprintf('Path 1 (Direct SINDy) RMSE:\n');
fprintf('  X (Prey): %.4f\n', rmse_path1_x);
fprintf('  Y (Predator): %.4f\n', rmse_path1_y);
fprintf('  Average: %.4f\n', (rmse_path1_x + rmse_path1_y)/2);

fprintf('\nPath 2 (GP-enhanced SINDy) RMSE:\n');
fprintf('  X (Prey): %.4f\n', rmse_path2_x);
fprintf('  Y (Predator): %.4f\n', rmse_path2_y);
fprintf('  Average: %.4f\n', (rmse_path2_x + rmse_path2_y)/2);

% True Lotka-Volterra coefficients (for reference)
% dx/dt = alpha*x - beta*x*y = 1.5*x - 1.0*x*y
% dy/dt = delta*x*y - gamma*y = 1.0*x*y - 3.0*y
% In library form [1, x, y, x^2, y^2, x*y]:
% dx/dt = 0*1 + 1.5*x + 0*y + 0*x^2 + 0*y^2 - 1.0*x*y
% dy/dt = 0*1 + 0*x - 3.0*y + 0*x^2 + 0*y^2 + 1.0*x*y
fprintf('\nTrue Lotka-Volterra coefficients (library: [1, x, y, x^2, y^2, x*y]):\n');
fprintf('  X: [0, %.2f, 0, 0, 0, %.2f]\n', alpha, -beta);
fprintf('  Y: [0, 0, %.2f, 0, 0, %.2f]\n', -gamma, delta);

fprintf('\nPath 1 recovered coefficients:\n');
fprintf('  X: [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', epsguess_path1(:,1)');
fprintf('  Y: [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', epsguess_path1(:,2)');

fprintf('\nPath 2 recovered coefficients (first 6 terms of expanded library):\n');
fprintf('  X: [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', epsguess_path2(1:6,1)');
fprintf('  Y: [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', epsguess_path2(1:6,2)');

fprintf('\n========== Analysis Complete ==========\n');
fprintf('Check figures 1-3 for visualization\n');

%% ========== Helper Functions ==========

function dydt = lotka_volterra_ode(t, y, alpha, beta, gamma, delta)
    % Lotka-Volterra ODE
    x = y(1);
    y_val = y(2);
    dydt = [alpha*x - beta*x*y_val;  % dx/dt
            delta*x*y_val - gamma*y_val];  % dy/dt
end

function dydt = sindy_rhs_path1(t, y, epsguess)
    % SINDy right-hand side for Path 1 (continuous-time ODE)
    % Model: dX/dt = f(X) learned using finite differences
    x = y(1);
    y_val = y(2);
    
    % Library functions: [1, x, y, x^2, y^2, x*y]
    library = [1; x; y_val; x^2; y_val^2; x*y_val];
    
    % Direct continuous-time ODE: dX/dt = f(X)
    dxdt = epsguess(:,1)' * library;
    dydt_val = epsguess(:,2)' * library;
    
    dydt = [dxdt; dydt_val];
end

function dydt = sindy_rhs_path2(t, y, epsguess)
    % SINDy right-hand side for Path 2 (expanded library)
    x = y(1);
    y_val = y(2);
    
    % Library functions: [1, x, y, x^2, y^2, x*y, x^3, y^3, x^2*y, x*y^2]
    library = [1; x; y_val; x^2; y_val^2; x*y_val; x^3; y_val^3; x^2*y_val; x*y_val^2];
    
    dxdt = epsguess(:,1)' * library;
    dydt_val = epsguess(:,2)' * library;
    
    dydt = [dxdt; dydt_val];
end

function Xi = sindy_stls(Theta, dXdt, lambda, max_iter)
    % Sequential Thresholded Least Squares (STLS) for SINDy
    % Theta: library matrix, dXdt: derivative targets
    % lambda: threshold, max_iter: number of thresholding iterations
    Xi = Theta \ dXdt;  % Initial least-squares fit
    for k = 1:max_iter
        small = abs(Xi) < lambda;
        Xi(small) = 0;
        for i = 1:size(dXdt, 2)
            big = ~small(:, i);
            if any(big)
                Xi(big, i) = Theta(:, big) \ dXdt(:, i);
            end
        end
    end
end
