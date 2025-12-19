%% Visualization Generator: Ground Truth vs. Path 1 (Discrete) vs. Path 2 (GP Continuous)
clear; clc; close all;

%% 1. Generate Sparse Noisy Lotka-Volterra Data
t_fine = linspace(0, 15, 100); 
params = [1.1, 0.4, 0.1, 0.4]; % alpha, beta, delta, gamma
lv_ode = @(t, x) [params(1)*x(1) - params(2)*x(1)*x(2); ...
                  params(3)*x(1)*x(2) - params(4)*x(2)];
[~, X_true] = ode45(lv_ode, t_fine, [10; 5]);

% Sample exactly 6 points (Noisy and Sparse)
sample_idx = round(linspace(1, 100, 6)); 
t_obs = t_fine(sample_idx)';
X_obs = X_true(sample_idx, :) + normrnd(0, 0.4, [6, 2]);

%% 2. PATH 1: ESINDy on Sparse Raw Data (Discrete Map - Author Style)
% Matches ESINDy.m logic: X_{k+1} = f(X_k)
N = 1000;
lambda_sindy = 0.1;

% Build Library from first 5 points to predict the next 5 points
M1_raw = X_obs(:,1); M2_raw = X_obs(:,2);
Theta1 = [ones(5, 1), M1_raw(1:5), M2_raw(1:5), M1_raw(1:5).^2, M2_raw(1:5).^2, M1_raw(1:5).*M2_raw(1:5)];
Y1 = [M1_raw(2:6), M2_raw(2:6)]; % Target is the "Next" state

% Run Ensemble STLS (Using Median for robustness)
Xi_Path1 = run_Ensemble_STLS(Theta1, Y1, N, lambda_sindy);

%% 3. PATH 2: GP-Augmented Continuous ESINDy (Research Alternative)
% Step A: Fit GP to 6 points and "hallucinate" 100 simulated points
n_sim = 100;
t_sim_gp = linspace(min(t_obs), max(t_obs), n_sim)';
[X_sim, dX_sim] = fit_GP_and_Deriv(t_obs, X_obs, t_sim_gp);

% Step B: Build an EXPANDED Library (Continuous Time: dX/dt = f(X))
Theta2 = [ones(n_sim, 1), X_sim(:,1), X_sim(:,2), ...
          X_sim(:,1).^2, X_sim(:,2).^2, X_sim(:,1).*X_sim(:,2)];

Xi_Path2 = run_Ensemble_STLS(Theta2, dX_sim, N, lambda_sindy);

%% 4. Model Reconstruction
t_recon = linspace(0, 15, 200);
x0 = X_obs(1,:)';

% --- Path 1: Discrete Map Iteration (Author Style Plotting) ---
X_rec1 = zeros(length(t_recon), 2);
X_rec1(1,:) = x0';
for k = 1:length(t_recon)-1
    curr = X_rec1(k,:);
    lib_vec = [1, curr(1), curr(2), curr(1)^2, curr(2)^2, curr(1)*curr(2)];
    X_rec1(k+1,:) = lib_vec * Xi_Path1; % Direct discrete step
end

% --- Path 2: Continuous ODE Integration (GP Style) ---
rhs_path2 = @(t,x) ([1, x(1), x(2), x(1)^2, x(2)^2, x(1)*x(2)] * Xi_Path2)';
[~, X_rec2] = ode45(rhs_path2, t_recon, x0);

%% 5. Comparative Plotting
figure('Color', 'w', 'Position', [100, 100, 900, 600]);

% Subplot 1: Path 1 (Discrete Map)
subplot(2,1,1);
plot(t_recon, X_rec1(:,1), 'b--', 'LineWidth', 1.5); hold on;
plot(t_recon, X_rec1(:,2), 'r--', 'LineWidth', 1.5);
scatter(t_obs, X_obs(:,1), 80, 'b', 'filled', 'MarkerEdgeColor', 'k');
scatter(t_obs, X_obs(:,2), 80, 'r', 'filled', 'MarkerEdgeColor', 'k');
title('Path 1: Discrete Map Prediction (Author Style)');
legend('Prey (Rec)', 'Pred (Rec)', 'Raw 6 Samples'); grid on;

% Subplot 2: Path 2 (GP Continuous)
subplot(2,1,2);
plot(t_recon, X_rec2(:,1), 'b-', 'LineWidth', 2); hold on;
plot(t_recon, X_rec2(:,2), 'r-', 'LineWidth', 2);
scatter(t_obs, X_obs(:,1), 80, 'b', 'filled', 'MarkerEdgeColor', 'k');
scatter(t_obs, X_obs(:,2), 80, 'r', 'filled', 'MarkerEdgeColor', 'k');
title('Path 2: Continuous ODE Recovery (GP Augmented)');
legend('Prey (Rec)', 'Pred (Rec)', 'Raw 6 Samples'); grid on;

%% Helper Functions
function Xi = run_Ensemble_STLS(Theta, Y, N, lambda)
    num_lib = size(Theta, 2); num_vars = size(Y, 2);
    stored = zeros(num_lib, num_vars, N);
    for i = 1:N
        subset = sort(randsample(num_lib, num_lib-1)); 
        T_sub = Theta(:, subset);
        Xi_sub = T_sub \ Y;
        for k = 1:10
            small = abs(Xi_sub) < lambda; Xi_sub(small) = 0;
            for v = 1:num_vars
                big = ~small(:,v); Xi_sub(big,v) = T_sub(:,big) \ Y(:,v);
            end
        end
        temp = zeros(num_lib, num_vars); temp(subset, :) = Xi_sub;
        stored(:,:,i) = temp;
    end
    Xi = median(stored, 3); %
end

function [X_gp, dX_gp] = fit_GP_and_Deriv(t, X, t_new)
    X_gp = zeros(length(t_new), 2); dX_gp = zeros(length(t_new), 2);
    for i = 1:2
        mdl = fitrgp(t, X(:,i), 'KernelFunction', 'squaredexponential', ...
                     'Sigma', 0.4, 'Standardize', true, 'BasisFunction', 'constant');
        X_gp(:,i) = predict(mdl, t_new);
        L = mdl.KernelInformation.KernelParameters(1); sf = mdl.KernelInformation.KernelParameters(2);
        dist = t_new - mdl.X'; K = (sf^2) * exp(-0.5 * (dist.^2) / L^2);
        dK = K .* (-dist / L^2); dX_gp(:,i) = dK * mdl.Alpha;
    end
end