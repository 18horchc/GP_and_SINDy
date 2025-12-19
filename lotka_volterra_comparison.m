% Comparison of ESINDy on sparse noisy data (Path 1) vs GP-enhanced ESINDy (Path 2)
% For Lotka-Volterra equations
%
% Path 1: Direct ESINDy on 6 noisy sparse points
% Path 2: Fit GP to 6 points -> sample more data -> use GP derivative -> ESINDy with expanded library

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

% Select 6 sparse time points
t_sparse_unique = linspace(0, 10, 6);  % 6 evenly spaced time points
x_sparse_true_unique = interp1(t_sol, sol(:,1), t_sparse_unique);
y_sparse_true_unique = interp1(t_sol, sol(:,2), t_sparse_unique);

% Generate 3 independent noisy measurements at each time point
% This gives us 18 total data points (3 at each of 6 time points)
num_replicates = 3;  % Number of measurements per time point
noise_level = 0.1;  % 10% noise

% Initialize arrays for all 18 points
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
    for j = 1:num_replicates
        idx = (i-1)*num_replicates + j;
        noise_x = normrnd(0, noise_level * std(x_sparse_true_unique));
        noise_y = normrnd(0, noise_level * std(y_sparse_true_unique));
        
        x_sparse_noisy(idx) = x_true_i + noise_x;
        y_sparse_noisy(idx) = y_true_i + noise_y;
    end
end

% Ensure non-negative (Lotka-Volterra should be positive)
x_sparse_noisy = max(0.01, x_sparse_noisy);
y_sparse_noisy = max(0.01, y_sparse_noisy);

fprintf('Generated %d sparse noisy points from Lotka-Volterra\n', length(t_sparse));
fprintf('  - %d unique time points: [%s]\n', length(t_sparse_unique), num2str(t_sparse_unique));
fprintf('  - %d replicates per time point\n', num_replicates);
fprintf('  - Total: %d data points\n', length(t_sparse));

%% ========== PATH 1: Direct ESINDy on averaged data ==========
fprintf('\n========== PATH 1: Direct ESINDy ==========\n');

% For Path 1, we use AVERAGED data at each unique time point
% This gives us 6 points (one average per time point) instead of 18 raw points
% Average replicates at each unique time point to get 6 sequential points
xdata_path1 = zeros(size(t_sparse_unique));
ydata_path1 = zeros(size(t_sparse_unique));

for i = 1:length(t_sparse_unique)
    idx_at_time = abs(t_sparse - t_sparse_unique(i)) < 1e-6;
    xdata_path1(i) = mean(x_sparse_noisy(idx_at_time));
    ydata_path1(i) = mean(y_sparse_noisy(idx_at_time));
end

xdata_path1 = xdata_path1';
ydata_path1 = ydata_path1';

fprintf('Path 1 uses averaged data: %d points (averages at %d unique time points)\n', ...
        length(xdata_path1), length(t_sparse_unique));

% ESINDy uses discrete-time formulation: X(t+1) = f(X(t))
% NOTE: No derivatives are computed in Path 1!
% Instead, we directly use the next state X(t+1) as the target
% Since we only have 6 points, we have 5 transitions
%
% Theta_path1: Library functions evaluated at X(t) [current state]
% X2_path1:    X(t+1) [next state] - this is what we're trying to predict
Theta_path1 = [ones(length(xdata_path1)-1, 1) ...
               xdata_path1(1:end-1) ...
               ydata_path1(1:end-1) ...
               xdata_path1(1:end-1).^2 ...
               ydata_path1(1:end-1).^2 ...
               xdata_path1(1:end-1).*ydata_path1(1:end-1)];

X2_path1 = [xdata_path1(2:end) ydata_path1(2:end)];

% Run ESINDy ensemble
N = 1000;  % Number of ensemble members
epsguessstored_path1 = zeros(6, 2, N);
choices = [1:1:6];

for i = 1:N
    l = 5;  % Use 5 out of 6 candidate functions
    total_cols = 6;
    p = randsample(total_cols, l);
    p = sort(p);
    g = ismember(choices, p);
    
    % Find which function was excluded
    stored_idx = find(g == 0);
    if isempty(stored_idx)
        stored_idx = 1;  % Fallback
    else
        stored_idx = stored_idx(1);
    end
    
    % Create reduced Theta matrix
    Thetait = Theta_path1(:, p);
    X2it = X2_path1;
    
    lambda = 0.01;
    
    % STLS (Sequential Thresholded Least Squares)
    epsguess = Thetait \ X2it;
    
    for c = 1:10
        smallinds = (abs(epsguess) < lambda);
        epsguess(smallinds) = 0;
        biginds1 = ~smallinds(:,1);
        biginds2 = ~smallinds(:,2);
        epsguess(biginds1,1) = Thetait(:,biginds1) \ X2it(:,1);
        epsguess(biginds2,2) = Thetait(:,biginds2) \ X2it(:,2);
    end
    
    % Store coefficients (map back to full 6-function library)
    epsguessstored_path1(stored_idx, :, i) = NaN;
    if stored_idx > 1
        epsguessstored_path1(1:stored_idx-1, :, i) = epsguess(1:stored_idx-1, :);
    end
    if stored_idx < 6
        epsguessstored_path1(stored_idx+1:end, :, i) = epsguess(stored_idx:end, :);
    end
end

% Find most common model
epsguessstored_path1(isnan(epsguessstored_path1)) = 0;

% Get the first model as representative (user can modify this)
epsguess_path1 = epsguessstored_path1(:, :, 1);

fprintf('Path 1 Model Coefficients:\n');
fprintf('  X: [%s]\n', num2str(epsguess_path1(:,1)'));
fprintf('  Y: [%s]\n', num2str(epsguess_path1(:,2)'));

%% ========== PATH 2: GP-enhanced ESINDy ==========
fprintf('\n========== PATH 2: GP-enhanced ESINDy ==========\n');

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

% Using Squared Exponential kernel (as in Wang & Barber 2014)
% This kernel is well-suited for smooth functions and has stable analytical derivatives
% The paper uses GP gradient matching for ODE parameter estimation, similar to our approach
% Note: The paper uses constant basis functions for simplicity
fprintf('Fitting GP models with Squared Exponential kernel (Wang & Barber approach)...\n');
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

% Generate more time points for sampling
t_dense = linspace(0, 10, 50);  % 50 points instead of 6
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

% Sample ONE GP curve from the posterior distribution per run
% This gives us a full curve that respects the GP covariance structure
% Each run will sample a different curve, providing variability
num_gp_samples = 1;  % One curve per run (as requested)
y1_gp_samples = sample_gp_posterior(gprMdl_x, X_dense, num_gp_samples);
y2_gp_samples = sample_gp_posterior(gprMdl_y, X_dense, num_gp_samples);

% Use the sampled curve
y1_sampled = y1_gp_samples(:, 1);
y2_sampled = y2_gp_samples(:, 1);

% Ensure non-negative
y1_sampled = max(0.01, y1_sampled);
y2_sampled = max(0.01, y2_sampled);

% Combine GP mean AND samples for data usage in Path 2
% Strategy: Use GP mean at some points, GP samples at others
% This gives us both the smooth mean estimate and the variability from sampling
t_combined = t_dense;  % Use dense time grid

% Interleave GP mean and samples: use mean at even indices, samples at odd indices
% This gives us a mix of both
x_combined = zeros(size(t_dense));
y_combined = zeros(size(t_dense));
for i = 1:length(t_dense)
    if mod(i, 2) == 0
        % Even indices: use GP mean
        x_combined(i) = y1pred(i);
        y_combined(i) = y2pred(i);
    else
        % Odd indices: use GP sample
        x_combined(i) = y1_sampled(i);
        y_combined(i) = y2_sampled(i);
    end
end

% Replace with original noisy data at original time points (use mean of replicates)
% This ensures we keep the actual observed data
t_sparse_unique = unique(t_sparse);
for i = 1:length(t_sparse_unique)
    t_i = t_sparse_unique(i);
    [~, closest_idx] = min(abs(t_combined - t_i));
    
    % Find all points at this time and use their mean
    idx_at_time = abs(t_sparse - t_i) < 1e-6;
    x_combined(closest_idx) = mean(x_sparse_noisy(idx_at_time));
    y_combined(closest_idx) = mean(y_sparse_noisy(idx_at_time));
end

fprintf('Combined dataset: %d points (%d original + %d GP-sampled)\n', ...
        length(t_combined), length(t_sparse_unique), length(t_combined) - length(t_sparse_unique));

% Compute GP derivative analytically
% The derivative of a GP is itself a GP! We can compute it by differentiating
% the covariance kernel directly, which is more principled than finite differences.
%
% Mathematical foundation:
% For a squared exponential kernel k(x, x') = σ² exp(-(x-x')²/(2l²)):
% - The derivative GP has covariance: k'(x, x') = σ²/l² * (1 - (x-x')²/l²) * exp(-(x-x')²/(2l²))
% - The mean of the derivative GP: μ' = K'(X_new, X) * K(X,X)^(-1) * y
%
% How this fits into ESINDy workflow:
% 1. We fit GP to sparse noisy data → get smooth function estimate
% 2. We compute analytical derivative of GP → get smooth derivative estimate
% 3. We use these derivatives in continuous-time SINDy: dX/dt = f(X)
% 4. The GP derivative provides denoised, uncertainty-quantified derivatives
%    that are more reliable than finite differences, especially with sparse data
%
% This gives us the analytical derivative at any point, with proper uncertainty
% quantification from the GP framework.
[y1_deriv_gp, y1_deriv_var] = gp_derivative_analytical(gprMdl_x, X_dense);
[y2_deriv_gp, y2_deriv_var] = gp_derivative_analytical(gprMdl_y, X_dense);

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

% Derivatives are already computed on t_dense, which matches t_combined
y1_deriv_combined = y1_deriv_gp;
y2_deriv_combined = y2_deriv_gp;

% For ESINDy with derivatives, we'll use a continuous-time approach
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

% Run ESINDy ensemble with expanded library
num_functions_path2 = size(Theta_path2, 2);
epsguessstored_path2 = zeros(num_functions_path2, 2, N);
choices_path2 = [1:1:num_functions_path2];

for i = 1:N
    l = num_functions_path2 - 1;  % Use all but one function
    total_cols = num_functions_path2;
    p = randsample(total_cols, l);
    p = sort(p);
    g = ismember(choices_path2, p);
    
    stored_idx = find(g == 0);
    if isempty(stored_idx)
        stored_idx = 1;
    else
        stored_idx = stored_idx(1);
    end
    
    Thetait = Theta_path2(:, p);
    X2it = X_dot_path2;
    
    lambda = 0.01;
    
    epsguess = Thetait \ X2it;
    
    for c = 1:10
        smallinds = (abs(epsguess) < lambda);
        epsguess(smallinds) = 0;
        biginds1 = ~smallinds(:,1);
        biginds2 = ~smallinds(:,2);
        epsguess(biginds1,1) = Thetait(:,biginds1) \ X2it(:,1);
        epsguess(biginds2,2) = Thetait(:,biginds2) \ X2it(:,2);
    end
    
    epsguessstored_path2(stored_idx, :, i) = NaN;
    if stored_idx > 1
        epsguessstored_path2(1:stored_idx-1, :, i) = epsguess(1:stored_idx-1, :);
    end
    if stored_idx < num_functions_path2
        epsguessstored_path2(stored_idx+1:end, :, i) = epsguess(stored_idx:end, :);
    end
end

epsguessstored_path2(isnan(epsguessstored_path2)) = 0;
epsguess_path2 = epsguessstored_path2(:, :, 1);

fprintf('Path 2 Model Coefficients (expanded library with %d functions):\n', num_functions_path2);
fprintf('  X: [%s]\n', num2str(epsguess_path2(:,1)'));
fprintf('  Y: [%s]\n', num2str(epsguess_path2(:,2)'));

%% ========== Visualization ==========
fprintf('\n========== Generating Plots ==========\n');

% Plot 1: Original data and true solution
figure(1);
subplot(2,1,1);
plot(t_true, x_true, 'b-', 'LineWidth', 2, 'DisplayName', 'True X (prey)');
hold on;
scatter(t_sparse, x_sparse_noisy, 80, 'r', 'filled', 'DisplayName', 'Noisy sparse data (18 points)');
xlabel('Time');
ylabel('X (Prey)');
title('Lotka-Volterra: Prey Population (18 points: 3 replicates at 6 time points)');
legend('Location', 'best');
grid on;

subplot(2,1,2);
plot(t_true, y_true, 'b-', 'LineWidth', 2, 'DisplayName', 'True Y (predator)');
hold on;
scatter(t_sparse, y_sparse_noisy, 80, 'r', 'filled', 'DisplayName', 'Noisy sparse data (18 points)');
xlabel('Time');
ylabel('Y (Predator)');
title('Lotka-Volterra: Predator Population (18 points: 3 replicates at 6 time points)');
legend('Location', 'best');
grid on;

% Plot 2: GP fit and sampling
figure(2);
subplot(2,1,1);
plot(t_dense, y1pred, 'b-', 'LineWidth', 2, 'DisplayName', 'GP mean');
hold on;
scatter(t_sparse, x_sparse_noisy, 80, 'r', 'filled', 'DisplayName', 'Original 18 points');
scatter(t_dense, y1_sampled, 50, 'g', 'o', 'DisplayName', 'GP samples');
xlabel('Time');
ylabel('X (Prey)');
title('Path 2: GP Fit and Sampling (Prey)');
legend('Location', 'best');
grid on;

subplot(2,1,2);
plot(t_dense, y2pred, 'b-', 'LineWidth', 2, 'DisplayName', 'GP mean');
hold on;
scatter(t_sparse, y_sparse_noisy, 80, 'r', 'filled', 'DisplayName', 'Original 18 points');
scatter(t_dense, y2_sampled, 50, 'g', 'o', 'DisplayName', 'GP samples');
xlabel('Time');
ylabel('Y (Predator)');
title('Path 2: GP Fit and Sampling (Predator)');
legend('Location', 'best');
grid on;

% Plot 3: Model predictions comparison
figure(3);

% Simulate Path 1 model (discrete-time)
% CRITICAL: Use unique time points for simulation, not all 18 points!
% The model was trained on averaged data at 6 unique time points,
% so we must simulate at those same 6 unique time points.
% Using all 18 points would apply the discrete map multiple times at the same time,
% causing numerical blowup (values like 10^166).
t_sim_path1 = t_sparse_unique;  % Use 6 unique time points, not 18 with replicates
x_sim_path1 = zeros(size(t_sim_path1));
y_sim_path1 = zeros(size(t_sim_path1));

% Use mean of replicates at first time point for initial condition
idx_first = abs(t_sparse - t_sim_path1(1)) < 1e-6;
x_sim_path1(1) = mean(x_sparse_noisy(idx_first));
y_sim_path1(1) = mean(y_sparse_noisy(idx_first));

% Apply discrete-time map: X(t+1) = f(X(t))
% This advances from one unique time point to the next
for i = 2:length(t_sim_path1)
    x_prev = x_sim_path1(i-1);
    y_prev = y_sim_path1(i-1);
    x_sim_path1(i) = epsguess_path1(1,1) + epsguess_path1(2,1)*x_prev + ...
                     epsguess_path1(3,1)*y_prev + epsguess_path1(4,1)*x_prev^2 + ...
                     epsguess_path1(5,1)*y_prev^2 + epsguess_path1(6,1)*x_prev*y_prev;
    y_sim_path1(i) = epsguess_path1(1,2) + epsguess_path1(2,2)*x_prev + ...
                     epsguess_path1(3,2)*y_prev + epsguess_path1(4,2)*x_prev^2 + ...
                     epsguess_path1(5,2)*y_prev^2 + epsguess_path1(6,2)*x_prev*y_prev;
end

% Simulate Path 2 model (continuous-time ODE)
function_path2 = @(t, y) sindy_rhs_path2(t, y, epsguess_path2);
% Use mean of replicates at first time point for initial condition
idx_first = abs(t_sparse - t_sparse_unique(1)) < 1e-6;
x0_path2 = mean(x_sparse_noisy(idx_first));
y0_path2 = mean(y_sparse_noisy(idx_first));
[t_ode_path2, sol_path2] = ode45(function_path2, [0, 10], ...
                                 [x0_path2; y0_path2], ...
                                 odeset('RelTol', 1e-6));

subplot(2,2,1);
plot(t_true, x_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True');
hold on;
plot(t_sim_path1, x_sim_path1, 'r--', 'LineWidth', 2, 'DisplayName', 'Path 1');
% Path 1 uses averaged data, so show averaged points
scatter(t_sparse_unique, xdata_path1, 100, 'b', 'filled', 'DisplayName', 'Averaged data (6 pts)');
xlabel('Time');
ylabel('X (Prey)');
title('Path 1: ESINDy Prediction (uses averaged data)');
legend('Location', 'best');
grid on;

subplot(2,2,2);
plot(t_true, y_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True');
hold on;
plot(t_sim_path1, y_sim_path1, 'r--', 'LineWidth', 2, 'DisplayName', 'Path 1');
% Path 1 uses averaged data, so show averaged points
scatter(t_sparse_unique, ydata_path1, 100, 'b', 'filled', 'DisplayName', 'Averaged data (6 pts)');
xlabel('Time');
ylabel('Y (Predator)');
title('Path 1: ESINDy Prediction (uses averaged data)');
legend('Location', 'best');
grid on;

subplot(2,2,3);
plot(t_true, x_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True');
hold on;
plot(t_ode_path2, sol_path2(:,1), 'g--', 'LineWidth', 2, 'DisplayName', 'Path 2');
scatter(t_sparse, x_sparse_noisy, 80, 'b', 'filled', 'DisplayName', 'Data (18 pts)');
xlabel('Time');
ylabel('X (Prey)');
title('Path 2: GP-enhanced ESINDy Prediction');
legend('Location', 'best');
grid on;

subplot(2,2,4);
plot(t_true, y_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True');
hold on;
plot(t_ode_path2, sol_path2(:,2), 'g--', 'LineWidth', 2, 'DisplayName', 'Path 2');
scatter(t_sparse, y_sparse_noisy, 80, 'b', 'filled', 'DisplayName', 'Data (18 pts)');
xlabel('Time');
ylabel('Y (Predator)');
title('Path 2: GP-enhanced ESINDy Prediction');
legend('Location', 'best');
grid on;

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

fprintf('Path 1 (Direct ESINDy) RMSE:\n');
fprintf('  X (Prey): %.4f\n', rmse_path1_x);
fprintf('  Y (Predator): %.4f\n', rmse_path1_y);
fprintf('  Average: %.4f\n', (rmse_path1_x + rmse_path1_y)/2);

fprintf('\nPath 2 (GP-enhanced ESINDy) RMSE:\n');
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

