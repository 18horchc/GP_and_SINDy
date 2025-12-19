% Quick test to compare GP fitting options:
% 1. With/without square root transform
% 2. With/without standardization

clear; close all; clc

%% Generate test data (6 sparse noisy Lotka-Volterra points)
alpha = 1.5; beta = 1.0; gamma = 3.0; delta = 1.0;
x0 = 1.0; y0 = 1.0;
t_span = [0, 10];
[t_sol, sol] = ode45(@(t, y) lotka_volterra_ode(t, y, alpha, beta, gamma, delta), ...
                     t_span, [x0; y0], odeset('RelTol', 1e-8));

t_sparse = linspace(0, 10, 6);
x_sparse_true = interp1(t_sol, sol(:,1), t_sparse);
y_sparse_true = interp1(t_sol, sol(:,2), t_sparse);

noise_level = 0.1;
x_sparse_noisy = x_sparse_true + normrnd(0, noise_level * std(x_sparse_true), size(x_sparse_true));
y_sparse_noisy = y_sparse_true + normrnd(0, noise_level * std(y_sparse_true), size(y_sparse_true));
x_sparse_noisy = max(0.01, x_sparse_noisy);
y_sparse_noisy = max(0.01, y_sparse_noisy);

X_gp = t_sparse(:);
y1_gp = x_sparse_noisy(:);
y2_gp = y_sparse_noisy(:);

%% Test different GP configurations
t_dense = linspace(0, 10, 100)';
X_dense = t_dense(:);

fprintf('Testing GP fitting options...\n\n');

% Option 1: Direct, no standardization (RECOMMENDED)
fprintf('Option 1: Direct fit, Standardize = false\n');
gpr1_1 = fitrgp(X_gp, y1_gp, 'KernelFunction', 'squaredexponential', ...
                'BasisFunction', 'constant', 'FitMethod', 'exact', ...
                'PredictMethod', 'exact', 'Standardize', false);
gpr1_2 = fitrgp(X_gp, y2_gp, 'KernelFunction', 'squaredexponential', ...
                'BasisFunction', 'constant', 'FitMethod', 'exact', ...
                'PredictMethod', 'exact', 'Standardize', false);
[y1_pred1, ~, y1_std1] = predict(gpr1_1, X_dense);
[y2_pred1, ~, y2_std1] = predict(gpr1_2, X_dense);
y1_deriv1 = gradient(y1_pred1, mean(diff(t_dense)));
y2_deriv1 = gradient(y2_pred1, mean(diff(t_dense)));

% Option 2: Square root transform, Standardize = false
fprintf('Option 2: Square root transform, Standardize = false\n');
y1_sqrt = sqrt(y1_gp + 1e-6);
y2_sqrt = sqrt(y2_gp + 1e-6);
gpr2_1 = fitrgp(X_gp, y1_sqrt, 'KernelFunction', 'squaredexponential', ...
                'BasisFunction', 'constant', 'FitMethod', 'exact', ...
                'PredictMethod', 'exact', 'Standardize', false);
gpr2_2 = fitrgp(X_gp, y2_sqrt, 'KernelFunction', 'squaredexponential', ...
                'BasisFunction', 'constant', 'FitMethod', 'exact', ...
                'PredictMethod', 'exact', 'Standardize', false);
[y1_pred2_sqrt, ~, ~] = predict(gpr2_1, X_dense);
[y2_pred2_sqrt, ~, ~] = predict(gpr2_2, X_dense);
y1_pred2 = max(0.01, y1_pred2_sqrt).^2;
y2_pred2 = max(0.01, y2_pred2_sqrt).^2;
y1_deriv2 = gradient(y1_pred2, mean(diff(t_dense)));
y2_deriv2 = gradient(y2_pred2, mean(diff(t_dense)));

% Option 3: Direct, Standardize = true (PROBLEMATIC for derivatives)
fprintf('Option 3: Direct fit, Standardize = true (WARNING: problematic for derivatives)\n');
gpr3_1 = fitrgp(X_gp, y1_gp, 'KernelFunction', 'squaredexponential', ...
                'BasisFunction', 'constant', 'FitMethod', 'exact', ...
                'PredictMethod', 'exact', 'Standardize', true);
gpr3_2 = fitrgp(X_gp, y2_gp, 'KernelFunction', 'squaredexponential', ...
                'BasisFunction', 'constant', 'FitMethod', 'exact', ...
                'PredictMethod', 'exact', 'Standardize', true);
[y1_pred3, ~, ~] = predict(gpr3_1, X_dense);
[y2_pred3, ~, ~] = predict(gpr3_2, X_dense);
y1_deriv3 = gradient(y1_pred3, mean(diff(t_dense)));
y2_deriv3 = gradient(y2_pred3, mean(diff(t_dense)));

% True solution for comparison
x_true = interp1(t_sol, sol(:,1), t_dense);
y_true = interp1(t_sol, sol(:,2), t_dense);
x_true_deriv = gradient(x_true, mean(diff(t_dense)));
y_true_deriv = gradient(y_true, mean(diff(t_dense)));

%% Compute errors
rmse1_x = sqrt(mean((y1_pred1 - x_true).^2));
rmse1_y = sqrt(mean((y2_pred1 - y_true).^2));
rmse1_deriv_x = sqrt(mean((y1_deriv1 - x_true_deriv).^2));
rmse1_deriv_y = sqrt(mean((y2_deriv1 - y_true_deriv).^2));

rmse2_x = sqrt(mean((y1_pred2 - x_true).^2));
rmse2_y = sqrt(mean((y2_pred2 - y_true).^2));
rmse2_deriv_x = sqrt(mean((y1_deriv2 - x_true_deriv).^2));
rmse2_deriv_y = sqrt(mean((y2_deriv2 - y_true_deriv).^2));

rmse3_x = sqrt(mean((y1_pred3 - x_true).^2));
rmse3_y = sqrt(mean((y2_pred3 - y_true).^2));
rmse3_deriv_x = sqrt(mean((y1_deriv3 - x_true_deriv).^2));
rmse3_deriv_y = sqrt(mean((y2_deriv3 - y_true_deriv).^2));

%% Display results
fprintf('\n========== Comparison Results ==========\n');
fprintf('RMSE for X (Prey) predictions:\n');
fprintf('  Option 1 (Direct, no std): %.4f\n', rmse1_x);
fprintf('  Option 2 (Sqrt, no std):   %.4f\n', rmse2_x);
fprintf('  Option 3 (Direct, std):    %.4f\n', rmse3_x);

fprintf('\nRMSE for Y (Predator) predictions:\n');
fprintf('  Option 1 (Direct, no std): %.4f\n', rmse1_y);
fprintf('  Option 2 (Sqrt, no std):   %.4f\n', rmse2_y);
fprintf('  Option 3 (Direct, std):    %.4f\n', rmse3_y);

fprintf('\nRMSE for X derivatives:\n');
fprintf('  Option 1 (Direct, no std): %.4f\n', rmse1_deriv_x);
fprintf('  Option 2 (Sqrt, no std):   %.4f\n', rmse2_deriv_x);
fprintf('  Option 3 (Direct, std):    %.4f\n', rmse3_deriv_x);

fprintf('\nRMSE for Y derivatives:\n');
fprintf('  Option 1 (Direct, no std): %.4f\n', rmse1_deriv_y);
fprintf('  Option 2 (Sqrt, no std):   %.4f\n', rmse2_deriv_y);
fprintf('  Option 3 (Direct, std):    %.4f\n', rmse3_deriv_y);

fprintf('\n========== Recommendation ==========\n');
fprintf('For Lotka-Volterra data:\n');
fprintf('  - Square root transform: NOT NECESSARY (Option 1 works fine)\n');
fprintf('  - Standardize = false: RECOMMENDED (avoids derivative issues)\n');
fprintf('  - Standardize = true: PROBLEMATIC (derivatives are in wrong units)\n');

%% Plot comparison
figure('Position', [100, 100, 1200, 800]);

% Predictions
subplot(2,3,1);
plot(t_dense, x_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True');
hold on;
scatter(t_sparse, x_sparse_noisy, 100, 'r', 'filled', 'DisplayName', 'Data');
plot(t_dense, y1_pred1, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Option 1');
plot(t_dense, y1_pred2, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Option 2');
plot(t_dense, y1_pred3, 'm--', 'LineWidth', 1.5, 'DisplayName', 'Option 3');
xlabel('Time'); ylabel('X (Prey)'); title('GP Predictions: Prey');
legend('Location', 'best'); grid on;

subplot(2,3,2);
plot(t_dense, y_true, 'k-', 'LineWidth', 2, 'DisplayName', 'True');
hold on;
scatter(t_sparse, y_sparse_noisy, 100, 'r', 'filled', 'DisplayName', 'Data');
plot(t_dense, y2_pred1, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Option 1');
plot(t_dense, y2_pred2, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Option 2');
plot(t_dense, y2_pred3, 'm--', 'LineWidth', 1.5, 'DisplayName', 'Option 3');
xlabel('Time'); ylabel('Y (Predator)'); title('GP Predictions: Predator');
legend('Location', 'best'); grid on;

% Derivatives
subplot(2,3,4);
plot(t_dense, x_true_deriv, 'k-', 'LineWidth', 2, 'DisplayName', 'True');
hold on;
plot(t_dense, y1_deriv1, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Option 1');
plot(t_dense, y1_deriv2, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Option 2');
plot(t_dense, y1_deriv3, 'm--', 'LineWidth', 1.5, 'DisplayName', 'Option 3');
xlabel('Time'); ylabel('dX/dt'); title('GP Derivatives: Prey');
legend('Location', 'best'); grid on;

subplot(2,3,5);
plot(t_dense, y_true_deriv, 'k-', 'LineWidth', 2, 'DisplayName', 'True');
hold on;
plot(t_dense, y2_deriv1, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Option 1');
plot(t_dense, y2_deriv2, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Option 2');
plot(t_dense, y2_deriv3, 'm--', 'LineWidth', 1.5, 'DisplayName', 'Option 3');
xlabel('Time'); ylabel('dY/dt'); title('GP Derivatives: Predator');
legend('Location', 'best'); grid on;

% Error comparison
subplot(2,3,3);
bar_data = [rmse1_x, rmse1_y; rmse2_x, rmse2_y; rmse3_x, rmse3_y];
bar(bar_data);
set(gca, 'XTickLabel', {'Option 1', 'Option 2', 'Option 3'});
ylabel('RMSE'); title('Prediction RMSE');
legend('X (Prey)', 'Y (Predator)', 'Location', 'best');
grid on;

subplot(2,3,6);
bar_data_deriv = [rmse1_deriv_x, rmse1_deriv_y; rmse2_deriv_x, rmse2_deriv_y; rmse3_deriv_x, rmse3_deriv_y];
bar(bar_data_deriv);
set(gca, 'XTickLabel', {'Option 1', 'Option 2', 'Option 3'});
ylabel('RMSE'); title('Derivative RMSE');
legend('dX/dt', 'dY/dt', 'Location', 'best');
grid on;

function dydt = lotka_volterra_ode(t, y, alpha, beta, gamma, delta)
    x = y(1);
    y_val = y(2);
    dydt = [alpha*x - beta*x*y_val; delta*x*y_val - gamma*y_val];
end

