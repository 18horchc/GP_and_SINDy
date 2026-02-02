
clear; clc; close all; 

%% Toy problem Equations

%Log Growth ---------------------
% 1. Define parameters
r = 0.5;
K = 1000;
tspan = [0 20];
P0 = 10;

% 2. Define Logistic ODE: dP/dt = r*P*(1 - P/K) 
logisticODE = @(t, P) r * P * (1 - P / K);

% 3. Solve using ODE45 
[t, P] = ode45(logisticODE, tspan, P0);

% 4. Plot 
figure;
plot(t, P, '-o', 'LineWidth', 2);
title('Logistic Growth (ODE45 Simulation)');
xlabel('Time');
ylabel('Population');
grid on;

%Lotka Volterra ---------------------
function dPdt = lotka_volterra(t, P)
%LOTKA_VOLTERRA Defines the Lotka-Volterra predator-prey model
% P(1) is the prey population (x)
% P(2) is the predator population (y)

    % Define the parameters (these can be made input arguments for more flexibility)
    alpha = 1.1;  % Prey growth rate
    beta = 0.4;   % Predation rate
    delta = 0.4;  % Predator growth rate per prey eaten
    gamma = 0.4;  % Predator death rate

    x = P(1);
    y = P(2);

    dPdt = zeros(2,1);    % Pre-allocate dPdt
    dPdt(1) = alpha*x - beta*x*y;
    dPdt(2) = delta*x*y - gamma*y;
end


% 1. Define parameters and initial conditions (arbitrary values can be changed here)
P0 = [20; 5];      % Initial populations: [Prey; Predators]
tspan = [0 50];    % Time span for the simulation (e.g., 0 to 50 time units)

% 2. Solve the ODE using ode45
% The solver calls the function defined in lotka_volterra.m
[t, P] = ode45(@lotka_volterra, tspan, P0);

% Extract populations
prey = P(:, 1);
predators = P(:, 2);


% 3. Plot the results

% Plot population dynamics over time
figure;
subplot(2,1,1);
plot(t, prey, 'b-', 'LineWidth', 2);
hold on;
plot(t, predators, 'r-', 'LineWidth', 2);
xlabel('Time');
ylabel('Population');
legend('Prey (x)', 'Predators (y)');
title('Lotka-Volterra Population Dynamics Over Time');
grid on;

% Plot the phase space (predators vs. prey)
subplot(2,1,2);
plot(prey, predators, 'k-', 'LineWidth', 2);
xlabel('Prey Population');
ylabel('Predator Population');
title('Phase Space Plot (Predators vs. Prey)');
grid on;

%Reaction kinetic ---------------------

function stiff_kinetics()
    % Initial concentrations [A, B, C]
    C0 = [1.0, 0, 0]; 
    tspan = [0 100];
    
    % Use stiff solver ode15s
    [t, C] = ode15s(@(t, y) stiff_system(t, y), tspan, C0);
    
    % Plotting
    plot(t, C(:,1), '-r', t, C(:,2), '-g', t, C(:,3), '-b', 'LineWidth', 2);
    xlabel('Time'); ylabel('Concentration');
    legend('A', 'B', 'C');
    title('Stiff Kinetic Curve');
end

function dCdt = stiff_system(~, C)
    % Reaction rates: A -> B (fast) -> C (slow)
    k1 = 1e3;  % Fast
    k2 = 1e-2; % Slow
    
    % Rate equations
    dCdt = zeros(3,1);
    dCdt(1) = -k1 * C(1);      % dA/dt
    dCdt(2) =  k1 * C(1) - k2 * C(2); % dB/dt
    dCdt(3) =  k2 * C(2);      % dC/dt
end

%% ========================================================================
%  Microglial Cell Dynamics Models (Amato & Arnold, 2025)
%  Reference: https://github.com/anarnold/Data-Driven-Modeling-Microglia
%  Paper: https://arxiv.org/abs/2404.10915
%  ========================================================================

% Experimental data from MCAO-induced stroke studies (Table 1 in paper)
% Time points where experimental data was observed
exp_time = [0, 1, 2, 3, 5, 7, 14];
exp_M1 = [5, 27.5, 122.5, 139.8, 325, 445, 816.67];  % M1 cell counts
exp_M2 = [5, 78.33, 179.5, 126.4, 800, 319, 136.67]; % M2 cell counts

% Additional validation data point at Day 35 (Suenaga et al., 2015)
exp_time_day35 = 35;
exp_M1_day35 = 550;
exp_M2_day35 = 50;

%% Method 1: SINDy + MCMC (Continuous-Time ODE Model)
% -------------------------------------------------------------------------
% From equations (3.1)-(3.2) with MCMC posterior mean estimates (eq 3.3):
%   dM1/dt = 32.9650 - 0.1377*M1 + 0.3150*M2
%   dM2/dt = 78.5130 - 0.1550*M1 + 0.0217*M2
%
% The MCMC step estimates θ3, θ4, θ5 (most sensitive parameters from GSS)

% Parameters with MCMC posterior means (from Figure 8 and eq 3.3)
A_sindy = [-0.1377,  0.3150;    % [θ2, θ3]
           -0.1550,  0.0217];   % [θ5, θ6]
b_sindy = [32.9650; 78.5130];   % [θ1; θ4]

% Define ODE: dM/dt = A*M + b
sindy_ode = @(t, M) A_sindy * M + b_sindy;

% Time span and initial conditions
tspan_sindy = [0 50];
M0 = [5; 5];  % Initial conditions from experimental data at Day 0

% Solve ODE
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[t_sindy, M_sindy] = ode45(sindy_ode, tspan_sindy, M0, opts);

M1_sindy = M_sindy(:, 1);
M2_sindy = M_sindy(:, 2);

% Compute equilibrium point
M_eq_sindy = -A_sindy \ b_sindy;
fprintf('SINDy+MCMC Equilibrium: M1* = %.2f, M2* = %.2f cells/mm^2\n', ...
        M_eq_sindy(1), M_eq_sindy(2));

%% Method 2: ESINDy + DBN (Discrete-Time Model)
% -------------------------------------------------------------------------
% From equations (3.4)-(3.5):
%   M1(t+1) = 30.5432 + 0.8960*M1(t) + 0.3242*M2(t)
%   M2(t+1) = 83.4946 - 0.1096*M1(t) + 0.8692*M2(t)
%
% Coefficients averaged across 10 proxy datasets (Figure 11-12)

% ESINDy coefficients (θ1-θ6 from equations 3.4-3.5)
theta_esindy = [30.5432, 0.8960, 0.3242;   % [θ1, θ2, θ3] for M1
                83.4946, -0.1096, 0.8692]; % [θ4, θ5, θ6] for M2

% Time points (discrete, daily)
t_esindy = 0:1:50;
n_steps = length(t_esindy);

% Initialize arrays
M1_esindy = zeros(n_steps, 1);
M2_esindy = zeros(n_steps, 1);

% Initial conditions
M1_esindy(1) = 5;
M2_esindy(1) = 5;

% Iterate the discrete-time model
for i = 2:n_steps
    M1_esindy(i) = theta_esindy(1,1) + theta_esindy(1,2)*M1_esindy(i-1) + theta_esindy(1,3)*M2_esindy(i-1);
    M2_esindy(i) = theta_esindy(2,1) + theta_esindy(2,2)*M1_esindy(i-1) + theta_esindy(2,3)*M2_esindy(i-1);
end

%% Plot: SINDy+MCMC Model (Figure 10 style)
figure;
plot(t_sindy, M1_sindy, 'k-', 'LineWidth', 2); hold on;
plot(t_sindy, M2_sindy, 'r-', 'LineWidth', 2);
xlabel('Time (Days)', 'FontSize', 14);
ylabel('cells/mm^2', 'FontSize', 14);
title('SINDy+MCMC Model: Microglial Cell Dynamics', 'FontSize', 14);
legend('M1', 'M2', 'Location', 'best');
xlim([0 50]); ylim([0 1500]);
grid on;
set(gca, 'FontSize', 12);

%% Plot: ESINDy+DBN Model (Figure 13 style)
figure;
plot(t_esindy, M1_esindy, 'k-', 'LineWidth', 2); hold on;
plot(t_esindy, M2_esindy, 'r-', 'LineWidth', 2);
xlabel('Time (Days)', 'FontSize', 14);
ylabel('cells/mm^2', 'FontSize', 14);
title('ESINDy+DBN Model: Microglial Cell Dynamics', 'FontSize', 14);
legend('M1', 'M2', 'Location', 'best');
xlim([0 50]); ylim([0 1500]);
grid on;
set(gca, 'FontSize', 12);

%% Plot: Comparison of Both Methods (Figure 14 style)
figure;
subplot(1,2,1);
plot(t_sindy, M1_sindy, 'k-', 'LineWidth', 2); hold on;
plot(t_sindy, M2_sindy, 'r-', 'LineWidth', 2);
xlabel('Time (Days)', 'FontSize', 12);
ylabel('cells/mm^2', 'FontSize', 12);
title('SINDy+MCMC', 'FontSize', 14);
legend('M1', 'M2', 'Location', 'best');
xlim([0 50]); ylim([0 1500]);
grid on;

subplot(1,2,2);
plot(t_esindy, M1_esindy, 'k-', 'LineWidth', 2); hold on;
plot(t_esindy, M2_esindy, 'r-', 'LineWidth', 2);
xlabel('Time (Days)', 'FontSize', 12);
ylabel('cells/mm^2', 'FontSize', 12);
title('ESINDy+DBN', 'FontSize', 14);
legend('M1', 'M2', 'Location', 'best');
xlim([0 50]); ylim([0 1500]);
grid on;

sgtitle('Comparison: Data-Driven Models of Microglial Cell Dynamics (Amato & Arnold 2025)', 'FontSize', 14);

%% Display model equations
fprintf('\n=== Model Equations ===\n');
fprintf('SINDy+MCMC (continuous-time ODE):\n');
fprintf('  dM1/dt = %.4f + %.4f*M1 + %.4f*M2\n', b_sindy(1), A_sindy(1,1), A_sindy(1,2));
fprintf('  dM2/dt = %.4f + %.4f*M1 + %.4f*M2\n', b_sindy(2), A_sindy(2,1), A_sindy(2,2));
fprintf('\nESINDy+DBN (discrete-time):\n');
fprintf('  M1(t+1) = %.4f + %.4f*M1(t) + %.4f*M2(t)\n', theta_esindy(1,:));
fprintf('  M2(t+1) = %.4f + %.4f*M1(t) + %.4f*M2(t)\n', theta_esindy(2,:));
