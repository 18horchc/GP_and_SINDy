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

%Sara's model ---------------------
% SINDy + MCMC
%% Linear ODE system for M1(t), M2(t)
%   M1' = 32.9650 - 0.1377*M1 + 0.3157*M2
%   M2' = 70.3415 - 0.1569*M1 + 0.0217*M2

% Parameters
A = [-0.1377,  0.3157;
     -0.1569,  0.0217];
b = [32.9650; 70.3415];

% RHS: dM/dt = A*M + b
rhs = @(t,M) A*M + b;

% Time span and initial conditions (EDIT THESE)
tspan = [0 50];      % time range
M0    = [5; 5];      % [M1(0); M2(0)]

% Solve
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t,M] = ode45(rhs, tspan, M0, opts);

M1 = M(:,1);
M2 = M(:,2);

% Plot
figure;
plot(t, M1, 'LineWidth', 2); hold on;
plot(t, M2, 'LineWidth', 2);
grid on;
xlabel('t');
ylabel('State value');
legend('M1(t)','M2(t)','Location','best');
title('ODE solution for M1 and M2');

% (Optional) steady state (equilibrium) where M' = 0: A*M + b = 0
M_ss = -A\b;
fprintf('Steady state: M1* = %.4f, M2* = %.4f\n', M_ss(1), M_ss(2));



%ESINDy + DBN



