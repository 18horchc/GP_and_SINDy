
%ODE (Exponential Growth Model): dP/dt = rP
%Analytical Solution: P(t) = P0*e^(rt)


% Parameters
r = 0.1;  %growth rate
P0 = 10; %Initial Population
t = linspace(0, 50, 25)'; % 25 data points

% Analytical Solution (Clean Data)
P_clean = P0 * exp(r * t);

% Add a tiny bit of noise (to see if SINDy/GP can handle it)
P_noisy = P_clean + 0.5 * randn(size(P_clean));

% Plot to verify
plot(t, P_clean, 'k-', 'LineWidth', 2); hold on;
plot(t, P_noisy, 'ro', 'MarkerFaceColor', 'r');
legend('Analytical Truth', 'Noisy Data Points');
xlabel('Time (t)'); ylabel('Population (P)');
title('Toy Problem: Exponential Growth');