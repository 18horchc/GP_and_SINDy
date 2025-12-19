clear; close all; clc;

%% 1. Generate Noisy Lotka-Volterra Data
t_span = linspace(0, 15, 100); 
params = [1.1, 0.4, 0.1, 0.4]; % alpha, beta, delta, gamma
lv_ode = @(t, x) [params(1)*x(1) - params(2)*x(1)*x(2); ...
                  params(3)*x(1)*x(2) - params(4)*x(2)];
[t_true, X_true] = ode45(lv_ode, t_span, [10; 5]);

% Create "Sparse & Noisy" Observed Data
sample_idx = 1:5:length(t_true); 
t_obs = t_true(sample_idx);
X_obs = X_true(sample_idx, :) + normrnd(0, 0.5, [length(sample_idx), 2]);

%% 2. PATH 1: Polynomial Interpolation (Original Style)
t_proxy = linspace(0, 15, 6)'; % Original code used ~15-30 points
X_poly_proxy = zeros(length(t_proxy), 2);

for i = 1:2
    p = polyfit(t_obs, X_obs(:,i), 4); % 4th order fit
    X_poly_proxy(:,i) = polyval(p, t_proxy);
end

% Build Library for Path 1 (2nd Order)
% Theta = [1, x, y, x^2, y^2, xy]
Theta1 = [ones(size(X_poly_proxy,1),1), X_poly_proxy(:,1), X_poly_proxy(:,2), ...
          X_poly_proxy(:,1).^2, X_poly_proxy(:,2).^2, X_poly_proxy(:,1).*X_poly_proxy(:,2)];
% Finite Difference for Path 1
dX1 = diff(X_poly_proxy) ./ diff(t_proxy);
Theta1_reduced = Theta1(1:end-1, :);

%% 3. PATH 2: GP Interpolation (Enhanced Method)
% Here we upsample to 200 points to support a larger library
[t_gp, X_gp, dX_gp] = fitAndPlotGP_Multi(t_obs, X_obs, 200);

% Build Library for Path 2 
% [1, x, y, x^2, y^2, xy]
Theta2 = [ones(size(X_gp,1),1), X_gp(:,1), X_gp(:,2), ...
          X_gp(:,1).^2, X_gp(:,2).^2, X_gp(:,1).*X_gp(:,2)];

         % X_gp(:,1).^3, X_gp(:,2).^3, (X_gp(:,1).^2).*X_gp(:,2), X_gp(:,1).*(X_gp(:,2).^2)];

%% 4. Ensemble SINDy (Style Matched to your Code)
N = 500; 
lambda = 0.1;
[Xi_Path1] = run_ESINDy_Logic(Theta1_reduced, dX1, N, lambda);
[Xi_Path2] = run_ESINDy_Logic(Theta2, dX_gp, N, lambda);

fprintf('Path 1 (Poly) Discovered Coefficients:\n'); disp(Xi_Path1);
fprintf('Path 2 (GP) Discovered Coefficients:\n'); disp(Xi_Path2);

%% 5. INTEGRATE DISCOVERED MODELS
% We define a function to use the Xi coefficients with our library structure
t_sim = linspace(0, 15, 200);
x0 = X_true(1,:)'; % Ensure initial condition is a column vector

% Path 1 Reconstruction (2nd Order Library)
% Added ( )' to convert the row vector result into a column vector
rhs_path1 = @(t,x) ([ones(1,1), x(1), x(2), x(1)^2, x(2)^2, x(1)*x(2)] * Xi_Path1)';

[t_rec1, X_rec1] = ode45(rhs_path1, t_sim, x0);

% Path 2 Reconstruction (3rd Order Library)
% Added ( )' to convert the row vector result into a column vector
rhs_path2 = @(t,x) ([ones(1,1), x(1), x(2), x(1)^2, x(2)^2, x(1)*x(2)] * Xi_Path2)' 

                    %x(1)^3, x(2)^3, (x(1)^2)*x(2), x(1)*(x(2)^2)] * Xi_Path2)';

[t_rec2, X_rec2] = ode45(rhs_path2, t_sim, x0);

%% 6. PLOTTING RESULTS

% --- FIGURE 1: GROUND TRUTH ---
figure('Name', 'Ground Truth', 'Color', 'w');
plot(t_true, X_true(:,1), 'b-', 'LineWidth', 2); hold on;
plot(t_true, X_true(:,2), 'r-', 'LineWidth', 2);
title('Original Lotka-Volterra Dynamics (Noiseless)');
xlabel('Time'); ylabel('Population');
legend('Prey (True)', 'Predator (True)'); grid on;

% --- FIGURE 2: PATH 1 (Polynomial + Central Diff) ---
figure('Name', 'Path 1 Results', 'Color', 'w');
plot(t_rec1, X_rec1(:,1), 'b--', 'LineWidth', 1.5); hold on;
plot(t_rec1, X_rec1(:,2), 'r--', 'LineWidth', 1.5);
scatter(t_obs, X_obs(:,1), 40, 'b', 'filled', 'MarkerFaceAlpha', 0.4);
scatter(t_obs, X_obs(:,2), 40, 'r', 'filled', 'MarkerFaceAlpha', 0.4);
title('Path 1: Polynomial Proxy Recovery');
xlabel('Time'); ylabel('Population');
legend('Prey (Path 1)', 'Predator (Path 1)', 'Observed Prey', 'Observed Pred'); grid on;

% --- FIGURE 3: PATH 2 (GP + Analytical Derivative) ---
figure('Name', 'Path 2 Results', 'Color', 'w');
plot(t_rec2, X_rec2(:,1), 'b-', 'LineWidth', 2); hold on;
plot(t_rec2, X_rec2(:,2), 'r-', 'LineWidth', 2);
scatter(t_obs, X_obs(:,1), 40, 'b', 'filled', 'MarkerFaceAlpha', 0.4);
scatter(t_obs, X_obs(:,2), 40, 'r', 'filled', 'MarkerFaceAlpha', 0.4);
title('Path 2: GP + Analytical Derivative Recovery');
xlabel('Time'); ylabel('Population');
legend('Prey (Path 2)', 'Predator (Path 2)', 'Observed Prey', 'Observed Pred'); grid on;

%% Helper Functions (Style Matched)
function [Xi_final] = run_ESINDy_Logic(Theta, dX, N, lambda)
    num_lib = size(Theta, 2);
    num_vars = size(dX, 2);
    stored_coeffs = zeros(num_lib, num_vars, N);
    
    for i = 1:N
        % Randomly drop one candidate function (Ensemble logic)
        subset = randsample(num_lib, num_lib-1);
        Theta_sub = Theta(:, subset);
        
        % Sequentially Thresholded Least Squares (STLS)
        Xi_sub = Theta_sub \ dX;
        for k = 1:10
            small_idx = abs(Xi_sub) < lambda;
            Xi_sub(small_idx) = 0;
            for v = 1:num_vars
                big_idx = ~small_idx(:,v);
                Xi_sub(big_idx, v) = Theta_sub(:, big_idx) \ dX(:, v);
            end
        end
        % Re-map back to full library size
        temp = zeros(num_lib, num_vars);
        temp(subset, :) = Xi_sub;
        stored_coeffs(:,:,i) = temp;
    end
    Xi_final = median(stored_coeffs, 3); % Use median to get most likely model
end

function [t_new, X_gp, dX_gp] = fitAndPlotGP_Multi(t, X, n_points)
    t_new = linspace(min(t), max(t), n_points)';
    X_gp = zeros(n_points, 2); dX_gp = zeros(n_points, 2);
    for i = 1:2
        % Using the analytical derivative logic from our previous step
        gprMdl = fitrgp(t, X(:,i), 'KernelFunction', 'squaredexponential', 'Standardize', true);
        [mu, ~, ~] = predict(gprMdl, t_new);
        X_gp(:,i) = mu;
        
        % Analytical Derivative calc (simplified for brevity)
        L = gprMdl.KernelInformation.KernelParameters(1);
        sf = gprMdl.KernelInformation.KernelParameters(2);
        dist = t_new - gprMdl.X';
        K = (sf^2) * exp(-0.5 * (dist.^2) / L^2);
        dK = K .* (-dist / L^2);
        dX_gp(:,i) = dK * gprMdl.Alpha;
    end
end