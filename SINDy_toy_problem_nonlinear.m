% 2D nonlinear

% x is prey and y is pred
% x' = ax - bxy
% y' = dxy - gy




% Parameters
alpha = 1.1; beta = 0.4;
delta = 0.1; gamma = 0.4;
params = [alpha, beta, delta, gamma];

% Time span and initial conditions [Prey, Predators]
tspan = linspace(0, 50, 500);
y0 = [10, 5]; 

% Solve the system
[t, sol] = ode45(@(t, y) [params(1)*y(1) - params(2)*y(1)*y(2); ...
                          params(3)*y(1)*y(2) - params(4)*y(2)], tspan, y0);

% Plot the "Analytical Truth"
plot(t, sol(:,1), 'b', t, sol(:,2), 'r', 'LineWidth', 2);
legend('Prey (x)', 'Predators (y)');
xlabel('Time'); ylabel('Population');
title('Lotka-Volterra Toy Problem');