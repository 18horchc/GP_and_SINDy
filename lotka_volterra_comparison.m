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

% Parameters for Lotka-Volterra
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
t_sparse = linspace(0, 10, 6);  % 6 evenly spaced points
x_sparse_true = interp1(t_sol, sol(:,1), t_sparse);
y_sparse_true = interp1(t_sol, sol(:,2), t_sparse);

% Add Gaussian noise
noise_level = 0.1;  % 10% noise
noise_x = normrnd(0, noise_level * std(x_sparse_true), size(x_sparse_true));
noise_y = normrnd(0, noise_level * std(y_sparse_true), size(y_sparse_true));

x_sparse_noisy = x_sparse_true + noise_x;
y_sparse_noisy = y_sparse_true + noise_y;

% Ensure non-negative (Lotka-Volterra should be positive)
x_sparse_noisy = max(0.01, x_sparse_noisy);
y_sparse_noisy = max(0.01, y_sparse_noisy);

fprintf('Generated 6 sparse noisy points from Lotka-Volterra\n');
fprintf('Time points: [%s]\n', num2str(t_sparse));
fprintf('X values (noisy): [%s]\n', num2str(x_sparse_noisy));
fprintf('Y values (noisy): [%s]\n', num2str(y_sparse_noisy));

%% ========== PATH 1: Direct ESINDy on 6 noisy points ==========
fprintf('\n========== PATH 1: Direct ESINDy ==========\n');

% For ESINDy, we need discrete-time data
% We'll use the 6 points directly
xdata_path1 = x_sparse_noisy';
ydata_path1 = y_sparse_noisy';

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

% Fit GP models directly (no square root transform needed for Lotka-Volterra)
% 
% Why no square root transform?
% - The square root transform was likely for count data with variance proportional
%   to mean. For Lotka-Volterra, direct fitting works fine.
%
% Why Standardize = false?
% - When Standardize = true, fitrgp standardizes both predictors (time) and response
% - This complicates derivative computation: derivatives of standardized variables
%   are not the same as derivatives of original variables
% - For time series with derivatives, it's cleaner to work in original units
% - If needed, you can manually standardize, but then must account for it in derivatives
gprMdl_x = fitrgp(X_gp, y1_gp, ...
    'KernelFunction', 'squaredexponential', ...
    'BasisFunction', 'constant', ...
    'FitMethod', 'exact', ...
    'PredictMethod', 'exact', ...
    'Standardize', false);

gprMdl_y = fitrgp(X_gp, y2_gp, ...
    'KernelFunction', 'squaredexponential', ...
    'BasisFunction', 'constant', ...
    'FitMethod', 'exact', ...
    'PredictMethod', 'exact', ...
    'Standardize', false);

% Generate more time points for sampling
t_dense = linspace(0, 10, 50);  % 50 points instead of 6
X_dense = t_dense(:);

% Predict from GP (directly, no transform needed)
[y1pred, ~, ~] = predict(gprMdl_x, X_dense);
[y2pred, ~, ~] = predict(gprMdl_y, X_dense);

% Ensure non-negative
y1pred = max(0.01, y1pred);
y2pred = max(0.01, y2pred);

% Sample from GP (add some uncertainty)
% We can use the prediction intervals or add noise
y1_sampled = y1pred + 0.05 * std(y1pred) * randn(size(y1pred));
y2_sampled = y2pred + 0.05 * std(y2pred) * randn(size(y2pred));
y1_sampled = max(0.01, y1_sampled);
y2_sampled = max(0.01, y2_sampled);

% Combine original 6 points with new simulated data
% Use original points where available, GP samples elsewhere
t_combined = t_dense;  % Use dense time grid
x_combined = y1_sampled;  % Start with GP samples
y_combined = y2_sampled;

% Replace with original noisy data at original time points
for i = 1:length(t_sparse)
    [~, closest_idx] = min(abs(t_combined - t_sparse(i)));
    x_combined(closest_idx) = x_sparse_noisy(i);
    y_combined(closest_idx) = y_sparse_noisy(i);
end

fprintf('Combined dataset: %d points (6 original + %d GP-sampled)\n', ...
        length(t_combined), length(t_combined) - 6);

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
scatter(t_sparse, x_sparse_noisy, 100, 'r', 'filled', 'DisplayName', 'Noisy sparse data');
xlabel('Time');
ylabel('X (Prey)');
title('Lotka-Volterra: Prey Population');
legend('Location', 'best');
grid on;

subplot(2,1,2);
plot(t_true, y_true, 'b-', 'LineWidth', 2, 'DisplayName', 'True Y (predator)');
hold on;
scatter(t_sparse, y_sparse_noisy, 100, 'r', 'filled', 'DisplayName', 'Noisy sparse data');
xlabel('Time');
ylabel('Y (Predator)');
title('Lotka-Volterra: Predator Population');
legend('Location', 'best');
grid on;

% Plot 2: GP fit and sampling
figure(2);
subplot(2,1,1);
plot(t_dense, y1pred, 'b-', 'LineWidth', 2, 'DisplayName', 'GP mean');
hold on;
scatter(t_sparse, x_sparse_noisy, 100, 'r', 'filled', 'DisplayName', 'Original 6 points');
scatter(t_dense, y1_sampled, 50, 'g', 'o', 'DisplayName', 'GP samples');
xlabel('Time');
ylabel('X (Prey)');
title('Path 2: GP Fit and Sampling (Prey)');
legend('Location', 'best');
grid on;

subplot(2,1,2);
plot(t_dense, y2pred, 'b-', 'LineWidth', 2, 'DisplayName', 'GP mean');
hold on;
scatter(t_sparse, y_sparse_noisy, 100, 'r', 'filled', 'DisplayName', 'Original 6 points');
scatter(t_dense, y2_sampled, 50, 'g', 'o', 'DisplayName', 'GP samples');
xlabel('Time');
ylabel('Y (Predator)');
title('Path 2: GP Fit and Sampling (Predator)');
legend('Location', 'best');
grid on;

% Plot 3: Model predictions comparison
figure(3);

% Simulate Path 1 model (discrete-time)
t_sim_path1 = t_sparse;
x_sim_path1 = zeros(size(t_sim_path1));
y_sim_path1 = zeros(size(t_sim_path1));
x_sim_path1(1) = x_sparse_noisy(1);
y_sim_path1(1) = y_sparse_noisy(1);

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
[t_ode_path2, sol_path2] = ode45(function_path2, [0, 10], ...
                                 [x_sparse_noisy(1); y_sparse_noisy(1)], ...
                                 odeset('RelTol', 1e-6));

subplot(2,2,1);
plot(t_true, x_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True');
hold on;
plot(t_sim_path1, x_sim_path1, 'r--', 'LineWidth', 2, 'DisplayName', 'Path 1');
scatter(t_sparse, x_sparse_noisy, 100, 'b', 'filled', 'DisplayName', 'Data');
xlabel('Time');
ylabel('X (Prey)');
title('Path 1: ESINDy Prediction');
legend('Location', 'best');
grid on;

subplot(2,2,2);
plot(t_true, y_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True');
hold on;
plot(t_sim_path1, y_sim_path1, 'r--', 'LineWidth', 2, 'DisplayName', 'Path 1');
scatter(t_sparse, y_sparse_noisy, 100, 'b', 'filled', 'DisplayName', 'Data');
xlabel('Time');
ylabel('Y (Predator)');
title('Path 1: ESINDy Prediction');
legend('Location', 'best');
grid on;

subplot(2,2,3);
plot(t_true, x_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True');
hold on;
plot(t_ode_path2, sol_path2(:,1), 'g--', 'LineWidth', 2, 'DisplayName', 'Path 2');
scatter(t_sparse, x_sparse_noisy, 100, 'b', 'filled', 'DisplayName', 'Data');
xlabel('Time');
ylabel('X (Prey)');
title('Path 2: GP-enhanced ESINDy Prediction');
legend('Location', 'best');
grid on;

subplot(2,2,4);
plot(t_true, y_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True');
hold on;
plot(t_ode_path2, sol_path2(:,2), 'g--', 'LineWidth', 2, 'DisplayName', 'Path 2');
scatter(t_sparse, y_sparse_noisy, 100, 'b', 'filled', 'DisplayName', 'Data');
xlabel('Time');
ylabel('Y (Predator)');
title('Path 2: GP-enhanced ESINDy Prediction');
legend('Location', 'best');
grid on;

%% ========== Model Comparison and Error Analysis ==========
fprintf('\n========== Model Comparison ==========\n');

% Interpolate true solution to comparison time points
x_true_path1 = interp1(t_sol, sol(:,1), t_sim_path1, 'linear', 'extrap');
y_true_path1 = interp1(t_sol, sol(:,2), t_sim_path1, 'linear', 'extrap');

x_true_path2 = interp1(t_sol, sol(:,1), t_ode_path2, 'linear', 'extrap');
y_true_path2 = interp1(t_sol, sol(:,2), t_ode_path2, 'linear', 'extrap');

% Compute RMSE for each path
rmse_path1_x = sqrt(mean((x_sim_path1 - x_true_path1).^2));
rmse_path1_y = sqrt(mean((y_sim_path1 - y_true_path1).^2));
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

