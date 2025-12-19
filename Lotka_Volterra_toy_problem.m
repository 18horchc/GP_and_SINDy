%% Visualization Generator: Ground Truth vs. Path A vs. Path B
clear; clc; close all;

% Add paths to your helper functions if they are in subfolders
% addpath('./blocks'); 

%% 1. Setup a Specific "Representative" Scenario
% We pick parameters where noise is significant enough that Central Diff
% struggles, but the underlying signal is strong enough for GP to recover.
p_true.alpha = 1.1; p_true.beta = 0.4; p_true.delta = 0.1; p_true.gamma = 0.4;
t_span_true = 0:0.05:30; 
x0 = [10; 5];

% "Difficult" Data Parameters
noise_sigma = 0.05;  % Significant noise
n_samples = 20;      % Sparse sampling
lambda_sindy = 0.1;  % SINDy threshold

% --- Generate Data ---
[t_true, X_true, ~, Xi_true] = LV_ground_truth(t_span_true, x0, p_true);
[t_m, X_m] = degrade_data(t_true, X_true, n_samples, noise_sigma);

%% 2. Run Path A: Raw Data + Central Difference + ESINDy
fprintf('Running Path A (Central Difference)...\n');
% 2a. Central Difference Derivative
dt_m = t_m(2) - t_m(1);
dX_central = zeros(size(X_m));
% Simple central difference for interior points
dX_central(2:end-1, :) = (X_m(3:end, :) - X_m(1:end-2, :)) / (2*dt_m);
% Trim edges (SINDy needs X and dXdt to match size)
X_A_in = X_m(2:end-1, :);
dX_A_in = dX_central(2:end-1, :);

% 2b. Run ESINDy
[Xi_A, ~] = run_ESINDy(X_A_in,3, dX_A_in, lambda_sindy, 50, 0.8);

% 2c. Integrate Discovered Model A
t_span_recon = [min(t_true), max(t_true)];
odefun_A = @(t,x) reconstruct_dynamics(t, x, Xi_A);
% Use Try-Catch because bad models might diverge to infinity quickly
try
    [tA_recon, XA_recon] = ode45(odefun_A, t_span_recon, x0);
catch
    tA_recon = t_span_recon(1); XA_recon = x0';
    warning('Model A diverged during integration.');
end

%% 3. Run Path B: GP Interpolation + ESINDy
fprintf('Running Path B (Gaussian Process)...\n');
% 3a. GP Interpolation & Analytical Derivative
% Upsample significantly for smooth derivatives (e.g., 200 points)
[t_gp, X_gp, dX_gp] = fitAndPlotGP(t_m, X_m, 200);


[mdl1, x1, p1] = fitAndPlotGP(t_1, y_1, 'GP of x', 'x', false);



% 3b. Run ESINDy on the GP-augmented data
[Xi_B, ~] = run_ESINDy(X_gp, dX_gp, lambda_sindy, 50, 0.9);
disp('Path B Coefficients:');
disp(Xi_B);

% 3c. Integrate Discovered Model B
odefun_B = @(t,x) reconstruct_dynamics(t, x, Xi_B);
try
    [tB_recon, XB_recon] = ode45(odefun_B, t_span_recon, x0);
catch
    tB_recon = t_span_recon(1); XB_recon = x0';
    warning('Model B diverged during integration.');
end

%% 4. Generate the Three Requested Plots

% --- Shared Plotting Settings ---
c_prey = [0 0.4470 0.7410];    % Blue
c_pred = [0.8500 0.3250 0.0980]; % Red
marker_alpha = 0.6;

% --- Plot 1: True Model + Noisy Data ---
figure('Color', 'w', 'Position', [100, 100, 800, 300]);
plot(t_true, X_true(:,1), '-', 'Color', c_prey, 'LineWidth', 2, 'DisplayName', 'True Prey (x)'); hold on;
plot(t_true, X_true(:,2), '-', 'Color', c_pred, 'LineWidth', 2, 'DisplayName', 'True Predator (y)');
% Noisy dots
scatter(t_m, X_m(:,1), 30, c_prey, 'filled', 'MarkerFaceAlpha', marker_alpha, 'DisplayName', 'Prey Samples');
scatter(t_m, X_m(:,2), 30, c_pred, 'filled', 'MarkerFaceAlpha', marker_alpha, 'DisplayName', 'Predator Samples');
title('1. Ground Truth Model & Noisy Measurements');
xlabel('Time'); ylabel('Population'); legend('Location', 'northeastoutside'); grid on; ylim([0 max(X_true(:))*1.2]);

% --- Plot 2: Path A Recovered Model + Noisy Data ---
figure('Color', 'w', 'Position', [100, 450, 800, 300]);
% Reconstructed Trajectory A
plot(tA_recon, XA_recon(:,1), '--', 'Color', c_prey, 'LineWidth', 2, 'DisplayName', 'Path A Model (Prey)'); hold on;
plot(tA_recon, XA_recon(:,2), '--', 'Color', c_pred, 'LineWidth', 2, 'DisplayName', 'Path A Model (Pred)');
% Noisy dots (kept same for reference)
scatter(t_m, X_m(:,1), 30, c_prey, 'filled', 'MarkerFaceAlpha', marker_alpha, 'HandleVisibility', 'off');
scatter(t_m, X_m(:,2), 30, c_pred, 'filled', 'MarkerFaceAlpha', marker_alpha, 'HandleVisibility', 'off');
title('2. Path A Recovery (Central Diff + ESINDy)');
xlabel('Time'); ylabel('Population'); legend('Location', 'northeastoutside'); grid on; ylim([0 max(X_true(:))*1.2]);

% --- Plot 3: Path B Recovered Model + Noisy Data ---
figure('Color', 'w', 'Position', [100, 800, 800, 300]);
% Reconstructed Trajectory B
plot(tB_recon, XB_recon(:,1), '-.', 'Color', c_prey, 'LineWidth', 2, 'DisplayName', 'Path B Model (Prey)'); hold on;
plot(tB_recon, XB_recon(:,2), '-.', 'Color', c_pred, 'LineWidth', 2, 'DisplayName', 'Path B Model (Pred)');
% Noisy dots (kept same for reference)
scatter(t_m, X_m(:,1), 30, c_prey, 'filled', 'MarkerFaceAlpha', marker_alpha, 'HandleVisibility', 'off');
scatter(t_m, X_m(:,2), 30, c_pred, 'filled', 'MarkerFaceAlpha', marker_alpha, 'HandleVisibility', 'off');
title('3. Path B Recovery (GP + ESINDy)');
xlabel('Time'); ylabel('Population'); legend('Location', 'northeastoutside'); grid on; ylim([0 max(X_true(:))*1.2]);


%% Helper Function: Reconstruct Dynamics for ODE45
function dxdt = reconstruct_dynamics(t, x, Xi)
    % This function converts the SINDy coefficient matrix Xi back into
    % derivative values given a current state x.
    % Assumes PolyOrder = 2 library: [1, x(1), x(2), x(1)^2, x(1)x(2), x(2)^2]
    
    X_state = x'; % Transpose to row vector for library build
    
    % Build library vector for this single state
    Theta_vec = [1, X_state(1), X_state(2), ...
                 X_state(1)^2, X_state(1)*X_state(2), X_state(2)^2];
             
    % Calculate derivatives: dXdt = Theta * Xi
    dxdt_row = Theta_vec * Xi;
    dxdt = dxdt_row'; % Transpose back to column vector for ode45
end