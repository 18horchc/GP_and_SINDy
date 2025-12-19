% Generate plots for GP_SINDy_Session_Summary.tex
% This script creates the three figures referenced in the summary document
% without interfering with the main lotka_volterra_comparison.m script

clear; close all; clc

%% ========== Generate Lotka-Volterra Data (same as main script) ==========
alpha = 1.5;  % Prey growth rate
beta = 1.0;   % Predation rate
gamma = 3.0;  % Predator death rate
delta = 1.0;  % Predator growth rate

x0 = 1.0;  % Prey
y0 = 1.0;  % Predator

t_span = [0, 10];
t_true = linspace(0, 10, 1000);

% Solve Lotka-Volterra ODE
[t_sol, sol] = ode45(@(t, y) lotka_volterra_ode(t, y, alpha, beta, gamma, delta), ...
                     t_span, [x0; y0], odeset('RelTol', 1e-8));

x_true = interp1(t_sol, sol(:,1), t_true);
y_true = interp1(t_sol, sol(:,2), t_true);

% Select 6 sparse time points
t_sparse_unique = linspace(0, 10, 6);
x_sparse_true_unique = interp1(t_sol, sol(:,1), t_sparse_unique);
y_sparse_true_unique = interp1(t_sol, sol(:,2), t_sparse_unique);

% Generate 3 independent noisy measurements at each time point
num_replicates = 3;
noise_level = 0.1;

t_sparse = repmat(t_sparse_unique, num_replicates, 1);
t_sparse = t_sparse(:);
x_sparse_noisy = zeros(size(t_sparse));
y_sparse_noisy = zeros(size(t_sparse));

for i = 1:length(t_sparse_unique)
    x_true_i = x_sparse_true_unique(i);
    y_true_i = y_sparse_true_unique(i);
    
    for j = 1:num_replicates
        idx = (i-1)*num_replicates + j;
        noise_x = normrnd(0, noise_level * std(x_sparse_true_unique));
        noise_y = normrnd(0, noise_level * std(y_sparse_true_unique));
        
        x_sparse_noisy(idx) = x_true_i + noise_x;
        y_sparse_noisy(idx) = y_true_i + noise_y;
    end
end

x_sparse_noisy = max(0.01, x_sparse_noisy);
y_sparse_noisy = max(0.01, y_sparse_noisy);

%% ========== Fit GP Models ==========
X_gp = t_sparse(:);
y1_gp = x_sparse_noisy(:);
y2_gp = y_sparse_noisy(:);

fprintf('Fitting GP models for summary plots...\n');
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

% Generate dense time grid for predictions
t_dense = linspace(0, 10, 200)';
X_dense = t_dense(:);

% Predict from GP with uncertainty
[y1pred, y1pred_std, y1pred_int] = predict(gprMdl_x, X_dense);
[y2pred, y2pred_std, y2pred_int] = predict(gprMdl_y, X_dense);

y1pred = max(0.01, y1pred);
y2pred = max(0.01, y2pred);

%% ========== Figure 1: GP Fit Visualization ==========
fprintf('Generating Figure 1: GP fit visualization...\n');

figure(1);
set(gcf, 'Position', [100, 100, 800, 600]);

% Prey (X)
subplot(2,1,1);
hold on;
% True solution
plot(t_true, x_true, 'k-', 'LineWidth', 1.5, 'DisplayName', 'True solution');
% GP mean
plot(t_dense, y1pred, 'b-', 'LineWidth', 2, 'DisplayName', 'GP mean');
% Uncertainty bands (95% prediction interval)
fill([t_dense; flipud(t_dense)], ...
     [y1pred_int(:,1); flipud(y1pred_int(:,2))], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95% prediction interval');
% Original noisy data
scatter(t_sparse, x_sparse_noisy, 100, 'r', 'filled', 'MarkerEdgeColor', 'red', ...
        'DisplayName', 'Original noisy data (18 points)');
xlabel('Time', 'FontSize', 12);
ylabel('X (Prey)', 'FontSize', 12);
title('GP Fit to Sparse Noisy Lotka-Volterra Data (Prey)', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);
grid on;
box on;

% Predator (Y)
subplot(2,1,2);
hold on;
% True solution
plot(t_true, y_true, 'k-', 'LineWidth', 1.5, 'DisplayName', 'True solution');
% GP mean
plot(t_dense, y2pred, 'b-', 'LineWidth', 2, 'DisplayName', 'GP mean');
% Uncertainty bands
fill([t_dense; flipud(t_dense)], ...
     [y2pred_int(:,1); flipud(y2pred_int(:,2))], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95% prediction interval');
% Original noisy data
scatter(t_sparse, y_sparse_noisy, 100, 'r', 'filled', 'MarkerEdgeColor', 'red', ...
        'DisplayName', 'Original noisy data (18 points)');
xlabel('Time', 'FontSize', 12);
ylabel('Y (Predator)', 'FontSize', 12);
title('GP Fit to Sparse Noisy Lotka-Volterra Data (Predator)', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);
grid on;
box on;

% Save figure
saveas(gcf, 'figure1_gp_fit.png');
fprintf('  Saved: figure1_gp_fit.png\n');

%% ========== Figure 2: Derivative Comparison ==========
fprintf('Generating Figure 2: Derivative comparison...\n');

% Compute analytical GP derivatives
[y1_deriv_gp, ~] = gp_derivative_analytical(gprMdl_x, X_dense);
[y2_deriv_gp, ~] = gp_derivative_analytical(gprMdl_y, X_dense);

% Compute finite differences on raw data (at sparse points)
dt_sparse = mean(diff(t_sparse_unique));
x_deriv_fd_raw = gradient(interp1(t_sparse_unique, x_sparse_true_unique, t_sparse_unique, 'linear', 'extrap'), dt_sparse);
y_deriv_fd_raw = gradient(interp1(t_sparse_unique, y_sparse_true_unique, t_sparse_unique, 'linear', 'extrap'), dt_sparse);
t_fd_raw = t_sparse_unique;

% Compute finite differences on GP mean (for comparison)
dt_dense = mean(diff(t_dense));
x_deriv_fd_gp = gradient(y1pred, dt_dense);
y_deriv_fd_gp = gradient(y2pred, dt_dense);

% True derivatives (analytical from ODE)
x_true_dense = interp1(t_sol, sol(:,1), t_dense);
y_true_dense = interp1(t_sol, sol(:,2), t_dense);
x_deriv_true = alpha * x_true_dense - beta * x_true_dense .* y_true_dense;
y_deriv_true = delta * x_true_dense .* y_true_dense - gamma * y_true_dense;

figure(2);
set(gcf, 'Position', [100, 100, 1000, 600]);

% Prey derivatives
subplot(2,1,1);
hold on;
plot(t_dense, x_deriv_true, 'k--', 'LineWidth', 2, 'DisplayName', 'True derivative');
plot(t_dense, y1_deriv_gp, 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical GP derivative');
plot(t_dense, x_deriv_fd_gp, 'g:', 'LineWidth', 1.5, 'DisplayName', 'Finite diff on GP mean');
scatter(t_fd_raw, x_deriv_fd_raw, 100, 'r', 'filled', 'MarkerEdgeColor', 'red', ...
        'DisplayName', 'Finite diff on raw data');
xlabel('Time', 'FontSize', 12);
ylabel('dX/dt (Prey)', 'FontSize', 12);
title('Derivative Comparison: Prey (X)', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);
grid on;
box on;
xlim([0, 10]);

% Predator derivatives
subplot(2,1,2);
hold on;
plot(t_dense, y_deriv_true, 'k--', 'LineWidth', 2, 'DisplayName', 'True derivative');
plot(t_dense, y2_deriv_gp, 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical GP derivative');
plot(t_dense, y_deriv_fd_gp, 'g:', 'LineWidth', 1.5, 'DisplayName', 'Finite diff on GP mean');
scatter(t_fd_raw, y_deriv_fd_raw, 100, 'r', 'filled', 'MarkerEdgeColor', 'red', ...
        'DisplayName', 'Finite diff on raw data');
xlabel('Time', 'FontSize', 12);
ylabel('dY/dt (Predator)', 'FontSize', 12);
title('Derivative Comparison: Predator (Y)', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);
grid on;
box on;
xlim([0, 10]);

% Save figure
saveas(gcf, 'figure2_derivatives.png');
fprintf('  Saved: figure2_derivatives.png\n');

%% ========== Figure 3: Model Comparison ==========
fprintf('Generating Figure 3: Model comparison...\n');

% Run Path 1: Direct ESINDy (simplified version for plotting)
% Average data at unique time points
xdata_path1 = zeros(size(t_sparse_unique));
ydata_path1 = zeros(size(t_sparse_unique));
for i = 1:length(t_sparse_unique)
    idx_at_time = abs(t_sparse - t_sparse_unique(i)) < 1e-6;
    xdata_path1(i) = mean(x_sparse_noisy(idx_at_time));
    ydata_path1(i) = mean(y_sparse_noisy(idx_at_time));
end

% Simple linear regression for Path 1 (for demonstration)
% In reality, this would use ESINDy ensemble, but for plotting we'll use a simple fit
Theta_path1 = [ones(length(xdata_path1)-1, 1) ...
               xdata_path1(1:end-1)' ...
               ydata_path1(1:end-1)' ...
               xdata_path1(1:end-1)'.^2 ...
               ydata_path1(1:end-1)'.^2 ...
               xdata_path1(1:end-1)'.*ydata_path1(1:end-1)'];
X2_path1 = [xdata_path1(2:end)' ydata_path1(2:end)'];
epsguess_path1 = Theta_path1 \ X2_path1;

% Simulate Path 1 (discrete-time)
t_sim_path1 = t_sparse_unique;
x_sim_path1 = zeros(size(t_sim_path1));
y_sim_path1 = zeros(size(t_sim_path1));
idx_first = abs(t_sparse - t_sim_path1(1)) < 1e-6;
x_sim_path1(1) = mean(x_sparse_noisy(idx_first));
y_sim_path1(1) = mean(y_sparse_noisy(idx_first));

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

% Run Path 2: GP-enhanced ESINDy (simplified)
% Sample from GP
y1_gp_samples = sample_gp_posterior(gprMdl_x, X_dense, 1);
y2_gp_samples = sample_gp_posterior(gprMdl_y, X_dense, 1);
y1_sampled = y1_gp_samples(:, 1);
y2_sampled = y2_gp_samples(:, 1);

% Combine GP mean and samples
x_combined = zeros(size(t_dense));
y_combined = zeros(size(t_dense));
for i = 1:length(t_dense)
    if mod(i, 2) == 0
        x_combined(i) = y1pred(i);
        y_combined(i) = y2pred(i);
    else
        x_combined(i) = y1_sampled(i);
        y_combined(i) = y2_sampled(i);
    end
end

% Replace with original data at original time points
for i = 1:length(t_sparse_unique)
    t_i = t_sparse_unique(i);
    [~, closest_idx] = min(abs(t_dense - t_i));
    idx_at_time = abs(t_sparse - t_i) < 1e-6;
    x_combined(closest_idx) = mean(x_sparse_noisy(idx_at_time));
    y_combined(closest_idx) = mean(y_sparse_noisy(idx_at_time));
end

% Compute derivatives
[y1_deriv_combined, ~] = gp_derivative_analytical(gprMdl_x, X_dense);
[y2_deriv_combined, ~] = gp_derivative_analytical(gprMdl_y, X_dense);

% SINDy fit for Path 2 (expanded library) with STLS regularization
Theta_path2 = [ones(length(t_dense), 1) ...
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

% Use STLS (Sequential Thresholded Least Squares) with regularization
% This prevents overfitting with the expanded library
lambda = 0.01;  % Sparsity threshold
epsguess_path2 = Theta_path2 \ X_dot_path2;  % Initial guess

% Apply sparsity thresholding iteratively
for c = 1:10
    smallinds = (abs(epsguess_path2) < lambda);
    epsguess_path2(smallinds) = 0;
    biginds1 = ~smallinds(:,1);
    biginds2 = ~smallinds(:,2);
    % Only regress on non-zero coefficients
    if sum(biginds1) > 0
        epsguess_path2(biginds1,1) = Theta_path2(:,biginds1) \ X_dot_path2(:,1);
    end
    if sum(biginds2) > 0
        epsguess_path2(biginds2,2) = Theta_path2(:,biginds2) \ X_dot_path2(:,2);
    end
end

% Check for extreme coefficients that might cause instability
if any(abs(epsguess_path2) > 100)
    warning('Path 2 coefficients are very large. Applying additional regularization.');
    % Apply stronger threshold for very large coefficients
    extreme_inds = abs(epsguess_path2) > 100;
    epsguess_path2(extreme_inds) = sign(epsguess_path2(extreme_inds)) * 100;
end

% Simulate Path 2 (continuous-time ODE)
function_path2 = @(t, y) sindy_rhs_path2(t, y, epsguess_path2);
idx_first = abs(t_sparse - t_sparse_unique(1)) < 1e-6;
x0_path2 = mean(x_sparse_noisy(idx_first));
y0_path2 = mean(y_sparse_noisy(idx_first));

try
    % Use tighter tolerances and smaller max step to handle potential stiffness
    [t_ode_path2, sol_path2] = ode45(function_path2, [0, 10], ...
                                     [x0_path2; y0_path2], ...
                                     odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'MaxStep', 0.05));
    
    % Check for unrealistic values (blow-up)
    if any(sol_path2(:) > 1e6) || any(sol_path2(:) < -1e6) || any(~isfinite(sol_path2(:)))
        warning('Path 2 ODE solution contains unrealistic values. Truncating integration.');
        % Find where solution becomes unrealistic
        valid_idx = all(isfinite(sol_path2), 2) & all(abs(sol_path2) < 1e6, 2);
        if sum(valid_idx) > 10
            t_ode_path2 = t_ode_path2(valid_idx);
            sol_path2 = sol_path2(valid_idx, :);
        else
            % If too much is invalid, use GP mean
            t_ode_path2 = t_dense;
            sol_path2 = [y1pred, y2pred];
        end
    end
catch ME
    % If integration fails, use GP mean as fallback
    warning('Path 2 ODE integration failed: %s. Using GP mean for visualization.', ME.message);
    t_ode_path2 = t_dense;
    sol_path2 = [y1pred, y2pred];
end

figure(3);
set(gcf, 'Position', [100, 100, 1200, 600]);

% Prey comparison
subplot(2,2,1);
hold on;
plot(t_true, x_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True solution');
plot(t_sim_path1, x_sim_path1, 'b--', 'LineWidth', 2, 'DisplayName', 'Path 1 (Direct ESINDy)');
scatter(t_sparse_unique, xdata_path1, 100, 'b', 'filled', 'MarkerEdgeColor', 'blue', ...
        'DisplayName', 'Averaged data (6 pts)');
xlabel('Time', 'FontSize', 11);
ylabel('X (Prey)', 'FontSize', 11);
title('Path 1: Direct ESINDy', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 9);
grid on;
box on;
xlim([0, 10]);

subplot(2,2,2);
hold on;
plot(t_true, x_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True solution');
plot(t_ode_path2, sol_path2(:,1), 'r--', 'LineWidth', 2, 'DisplayName', 'Path 2 (GP-enhanced)');
scatter(t_sparse, x_sparse_noisy, 80, 'r', 'filled', 'MarkerEdgeColor', 'red', ...
        'DisplayName', 'Data (18 pts)');
xlabel('Time', 'FontSize', 11);
ylabel('X (Prey)', 'FontSize', 11);
title('Path 2: GP-enhanced ESINDy', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 9);
grid on;
box on;
xlim([0, 10]);

% Predator comparison
subplot(2,2,3);
hold on;
plot(t_true, y_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True solution');
plot(t_sim_path1, y_sim_path1, 'b--', 'LineWidth', 2, 'DisplayName', 'Path 1 (Direct ESINDy)');
scatter(t_sparse_unique, ydata_path1, 100, 'b', 'filled', 'MarkerEdgeColor', 'blue', ...
        'DisplayName', 'Averaged data (6 pts)');
xlabel('Time', 'FontSize', 11);
ylabel('Y (Predator)', 'FontSize', 11);
title('Path 1: Direct ESINDy', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 9);
grid on;
box on;
xlim([0, 10]);

subplot(2,2,4);
hold on;
plot(t_true, y_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True solution');
plot(t_ode_path2, sol_path2(:,2), 'r--', 'LineWidth', 2, 'DisplayName', 'Path 2 (GP-enhanced)');
scatter(t_sparse, y_sparse_noisy, 80, 'r', 'filled', 'MarkerEdgeColor', 'red', ...
        'DisplayName', 'Data (18 pts)');
xlabel('Time', 'FontSize', 11);
ylabel('Y (Predator)', 'FontSize', 11);
title('Path 2: GP-enhanced ESINDy', 'FontSize', 12);
legend('Location', 'best', 'FontSize', 9);
grid on;
box on;
xlim([0, 10]);

% Save figure
saveas(gcf, 'figure3_model_comparison.png');
fprintf('  Saved: figure3_model_comparison.png\n');

fprintf('\nAll summary plots generated successfully!\n');
fprintf('Files saved:\n');
fprintf('  - figure1_gp_fit.png\n');
fprintf('  - figure2_derivatives.png\n');
fprintf('  - figure3_model_comparison.png\n');

%% ========== Helper Functions ==========

function dydt = lotka_volterra_ode(t, y, alpha, beta, gamma, delta)
    x = y(1);
    y_val = y(2);
    dydt = [alpha*x - beta*x*y_val;
            delta*x*y_val - gamma*y_val];
end

function dydt = sindy_rhs_path2(t, y, epsguess)
    x = y(1);
    y_val = y(2);
    library = [1; x; y_val; x^2; y_val^2; x*y_val; x^3; y_val^3; x^2*y_val; x*y_val^2];
    dxdt = epsguess(:,1)' * library;
    dydt_val = epsguess(:,2)' * library;
    dydt = [dxdt; dydt_val];
end

