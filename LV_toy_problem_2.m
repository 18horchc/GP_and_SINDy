clear; close all; clc;

%% 1. Generate Sparse Noisy Lotka-Volterra Data
t_fine = linspace(0, 15, 100); 
params = [1.1, 0.4, 0.1, 0.4]; % Ground truth alpha, beta, delta, gamma
lv_ode = @(t, x) [params(1)*x(1) - params(2)*x(1)*x(2); ...
                  params(3)*x(1)*x(2) - params(4)*x(2)];
[~, X_true] = ode45(lv_ode, t_fine, [10; 5]);

% Sample exactly 6 points
sample_idx = round(linspace(1, 100, 6)); 
t_obs = t_fine(sample_idx)';
X_obs = X_true(sample_idx, :) + normrnd(0, 0.4, [6, 2]);

%% 2. PATH 1: ESINDy on Sparse Raw Data (Discrete Map)
% Using the style from ESINDy.m: discrete steps
N = 1000;
M1_raw = X_obs(:,1); M2_raw = X_obs(:,2);
Theta1 = [ones(5, 1), M1_raw(1:5), M2_raw(1:5), M1_raw(1:5).^2, M2_raw(1:5).^2, M1_raw(1:5).*M2_raw(1:5)];
Y1 = [M1_raw(2:6), M2_raw(2:6)];

Xi_Path1 = run_Ensemble_STLS(Theta1, Y1, N, 0.1); % Returns discrete-time coefficients

%% 3. PATH 2: GP-Augmented Continuous ESINDy
% Step A: Fit GP to the 6 points and "hallucinate" 100 more points
n_sim = 100;
t_sim = linspace(min(t_obs), max(t_obs), n_sim)';
[X_sim, dX_sim] = fit_GP_and_Deriv(t_obs, X_obs, t_sim);

% Step B: Build an EXPANDED Library (3rd Order)
% [1, x, y, x^2, y^2, xy, x^3, y^3, x^2y, xy^2] - 10 terms
Theta2 = [ones(n_sim, 1), X_sim(:,1), X_sim(:,2), ...
          X_sim(:,1).^2, X_sim(:,2).^2, X_sim(:,1).*X_sim(:,2), ...
          X_sim(:,1).^3, X_sim(:,2).^3, (X_sim(:,1).^2).*X_sim(:,2), X_sim(:,1).*(X_sim(:,2).^2)];

Xi_Path2 = run_Ensemble_STLS(Theta2, dX_sim, N, 0.1); % Returns continuous-time coefficients

%% 4. Model Reconstruction and Results
fprintf('Path 1 (Discrete) Coefficients for M1:\n'); disp(Xi_Path1(:,1));
fprintf('Path 2 (Continuous) Coefficients for M1:\n'); disp(Xi_Path2(:,1));

% Plot Path 2 results (the focus of the GP improvement)
figure('Color','w');
subplot(2,1,1); hold on;
scatter(t_obs, X_obs(:,1), 100, 'ko', 'filled');
plot(t_sim, X_sim(:,1), 'r-', 'LineWidth', 2);
title('Path 2: GP Denoising & Upsampling (Species 1)'); grid on;

subplot(2,1,2); hold on;
plot(t_sim, dX_sim(:,1), 'b-', 'LineWidth', 2);
title('Path 2: Analytical Derivative derived from 6 points'); grid on;


%% 4. Model Reconstruction
t_sim_plot = linspace(0, 15, 300);
x0 = X_true(1,:)';

% --- Path 1: Discrete Map Iteration ---
% x_{k+1} = Theta * Xi
X_rec1 = zeros(length(t_sim_plot), 2);
X_rec1(1,:) = x0';
dt_step = t_sim_plot(2) - t_sim_plot(1);
for k = 1:length(t_sim_plot)-1
    x_curr = X_rec1(k,:);
    theta_curr = [1, x_curr(1), x_curr(2), x_curr(1)^2, x_curr(2)^2, x_curr(1)*x_curr(2)];
    X_rec1(k+1,:) = theta_curr * Xi_Path1; % Discrete step
end

% --- Path 2: Continuous ODE Integration ---
% dx/dt = Theta * Xi
rhs_path2 = @(t,x) ([ones(1,1), x(1), x(2), x(1)^2, x(2)^2, x(1)*x(2), ...
                    x(1)^3, x(2)^3, (x(1)^2)*x(2), x(1)*(x(2)^2)] * Xi_Path2)';
try
    [t_rec2, X_rec2] = ode45(rhs_path2, t_sim_plot, x0);
catch
    X_rec2 = nan(length(t_sim_plot), 2); t_rec2 = t_sim_plot;
end

%% 5. Plotting Results

% Colors & Styling
c_prey = [0 0.4470 0.7410]; c_pred = [0.8500 0.3250 0.0980];

% --- Figure 1: Ground Truth ---
figure('Name', '1. Ground Truth', 'Color', 'w', 'Position', [100 600 800 250]);
plot(t_fine, X_true(:,1), 'Color', c_prey, 'LineWidth', 2); hold on;
plot(t_fine, X_true(:,2), 'Color', c_pred, 'LineWidth', 2);
title('Ground Truth Lotka-Volterra'); legend('Prey', 'Predator'); grid on;

% --- Figure 2: Path 1 (Sparse Raw Data) ---
figure('Name', '2. Path 1 Results', 'Color', 'w', 'Position', [100 350 800 250]);
plot(t_sim_plot, X_rec1(:,1), '--', 'Color', c_prey, 'LineWidth', 1.5); hold on;
plot(t_sim_plot, X_rec1(:,2), '--', 'Color', c_pred, 'LineWidth', 1.5);
scatter(t_obs, X_obs(:,1), 80, c_prey, 'filled', 'MarkerEdgeColor', 'k');
scatter(t_obs, X_obs(:,2), 80, c_pred, 'filled', 'MarkerEdgeColor', 'k');
title('Path 1: Discrete eSINDy (6 Noisy Samples Only)');
legend('Path 1 Prey', 'Path 1 Pred', 'Raw Samples'); grid on;

% --- Figure 3: Path 2 (GP Augmented Data) ---
figure('Name', '3. Path 2 Results', 'Color', 'w', 'Position', [100 100 800 250]);
% Plot integrated model
plot(t_rec2, X_rec2(:,1), '-', 'Color', c_prey, 'LineWidth', 2); hold on;
plot(t_rec2, X_rec2(:,2), '-', 'Color', c_pred, 'LineWidth', 2);
% Plot simulated GP points (the "more data" you generated)
scatter(t_sim, X_sim(:,1), 15, c_prey, 'MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
scatter(t_sim, X_sim(:,2), 15, c_pred, 'MarkerFaceAlpha', 0.3, 'HandleVisibility', 'off');
% Plot original 6 points
scatter(t_obs, X_obs(:,1), 80, c_prey, 'filled', 'MarkerEdgeColor', 'k');
scatter(t_obs, X_obs(:,2), 80, c_pred, 'filled', 'MarkerEdgeColor', 'k');
title('Path 2: Continuous eSINDy (GP Augmented - 100 Sim Points)');
legend('Path 2 Prey', 'Path 2 Pred', 'Raw Samples'); grid on;

%% --- HELPER FUNCTIONS ---

function Xi_final = run_Ensemble_STLS(Theta, Y, N, lambda)
    % Matches logic of ESINDy.m ensemble
    num_lib = size(Theta, 2);
    num_vars = size(Y, 2);
    stored = zeros(num_lib, num_vars, N);
    
    for i = 1:N
        % Randomly sample a subset of candidate functions
        subset = sort(randsample(num_lib, num_lib-1)); 
        Theta_it = Theta(:, subset);
        Xi_it = Theta_it \ Y; % Least Squares
        
        % STLS refinement
        for k = 1:10
            small = abs(Xi_it) < lambda;
            Xi_it(small) = 0;
            for v = 1:num_vars
                big = ~small(:,v);
                Xi_it(big, v) = Theta_it(:, big) \ Y(:, v);
            end
        end
        % Store result in the correct library index
        temp = zeros(num_lib, num_vars);
        temp(subset, :) = Xi_it;
        stored(:,:,i) = temp;
    end
    Xi_final = median(stored, 3); % More robust than manual sorting
end

function [X_gp, dX_gp] = fit_GP_and_Deriv(t, X, t_new)
    X_gp = zeros(length(t_new), 2); dX_gp = zeros(length(t_new), 2);
    for i = 1:2
        % Standardize helps the GP stabilize on sparse data
        mdl = fitrgp(t, X(:,i), 'KernelFunction', 'squaredexponential', ...
                     'Sigma', 0.4, 'Standardize', true);
        X_gp(:,i) = predict(mdl, t_new);
        
        % Analytical Derivative Calculation
        L = mdl.KernelInformation.KernelParameters(1);
        sf = mdl.KernelInformation.KernelParameters(2);
        dist = t_new - mdl.X';
        K = (sf^2) * exp(-0.5 * (dist.^2) / L^2);
        dK = K .* (-dist / L^2);
        dX_gp(:,i) = dK * mdl.Alpha;
    end
end