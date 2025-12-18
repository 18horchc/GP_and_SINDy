clear; close all; clc;

%Think about adding second eqution in same to model M1 and M2
%% 1. Generate Noisy Logistic Data
r1 = 0.6; K1 = 100; x1_0 = 5;
r2 = 1.2; K2 = 500; x2_0 = 5;

t_dense = linspace(0, 20, 500);
t_sampled = linspace(0, 19, 20);
dt = t_sampled(2) - t_sampled(1); % Time step

%Real curves
x1_real = (K1 * x1_0 * exp(r1 * t_dense)) ./ (K1 + x1_0 * (exp(r1 * t_dense) - 1));
x2_real = (K2 * x2_0 * exp(r2 * t_dense)) ./ (K2 + x2_0 * (exp(r2 * t_dense) - 1));

% Logistic growth: x(t) = (K*x0*e^rt) / (K + x0(e^rt - 1))
x1_clean = (K1 * x1_0 * exp(r1 * t_sampled)) ./ (K1 + x1_0 * (exp(r1 * t_sampled) - 1));
x2_clean = (K2 * x2_0 * exp(r2 * t_sampled)) ./ (K2 + x2_0 * (exp(r2 * t_sampled) - 1));

% Add 5% Noise
noise_level = 0.05;
x1_noisy = x1_clean + noise_level * std(x1_clean) * randn(size(x1_clean));
x2_noisy = x2_clean + noise_level * std(x2_clean) * randn(size(x2_clean));

figure;
plot(t_sampled, x1_noisy, 'bx', 'MarkerFaceColor', 'b'); hold on;
plot(t_sampled, x2_noisy, 'ro', 'MarkerFaceColor', 'r'); 
grid on;
xlabel('Time'); ylabel('Population');
title('Logistic Growth: 20 Sampled Points');
legend('x1 Sampled Data', 'x2 Sampled Data');



figure;
% Plot the Real Curves (Lines)
plot(t_dense, x1_real, 'b-', 'LineWidth', 1.5); hold on;
plot(t_dense, x2_real, 'r-', 'LineWidth', 1.5); 
% Plot the Sampled Data (Markers)
plot(t_sampled, x1_noisy, 'bx', 'MarkerSize', 8); 
plot(t_sampled, x2_noisy, 'ro', 'MarkerSize', 6); 
grid on;
xlabel('Time'); ylabel('Population');
title('Logistic Growth: Analytical Truth vs. Noisy Samples');
legend('x1 Real Curve', 'x2 Real Curve', 'x1 Sampled Data', 'x2 Sampled Data', 'Location', 'best');

% Combine into a matrix (matching the xobs structure: points x variables)
xobs_noisy = [x1_noisy', x2_noisy'];

%% 2. Compute Derivatives for x1 and x2
dxobs = zeros(size(xobs_noisy));

% Endpoints (Indices 1-2): 3rd order forward difference
dxobs(1,:) = (-11/6*xobs_noisy(1,:) + 3*xobs_noisy(2,:) - 3/2*xobs_noisy(3,:) + 1/3*xobs_noisy(4,:)) / dt;
dxobs(2,:) = (-11/6*xobs_noisy(2,:) + 3*xobs_noisy(3,:) - 3/2*xobs_noisy(4,:) + 1/3*xobs_noisy(5,:)) / dt;

% Middle Section (Indices 3-18): 4th order centered difference
% Logic: (-1/12*x(i+2) + 2/3*x(i+1) - 2/3*x(i-1) + 1/12*x(i-2)) / dt
dxobs(3:18,:) = (-1/12*xobs_noisy(5:20,:) + 2/3*xobs_noisy(4:19,:) - 2/3*xobs_noisy(2:17,:) + 1/12*xobs_noisy(1:16,:)) / dt;

% Endpoints (Indices 19-20): 3rd order backward difference
dxobs(19,:) = (11/6*xobs_noisy(19,:) - 3*xobs_noisy(18,:) + 3/2*xobs_noisy(17,:) - 1/3*xobs_noisy(16,:)) / dt;
dxobs(20,:) = (11/6*xobs_noisy(20,:) - 3*xobs_noisy(19,:) + 3/2*xobs_noisy(18,:) - 1/3*xobs_noisy(17,:)) / dt;
%% SINDy

X_dot = dxobs; %Derivative information; rhs of 2.1


n = 2;
polys = 0:2;
gamma = 0; 
common_params = {polys,[]};

tol_ode = 1e-10;         % set tolerance (abs and rel) of ode45
options = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,n));

Theta = build_theta(xobs_noisy,common_params);

lambda = 0.01; %threshhold value

%STLS
epsguess = Theta\X_dot; % initial guess: Least-squares

for k=1:100
smallinds = (abs(epsguess)<lambda); % find small coefficients
epsguess(smallinds)=0; % and threshold (set small coeff = 0)
for ind = 1:2 % n is state dimension (in this case we have M1 and M2)
biginds = ~smallinds(:,ind);
% Regress dynamics onto remaining terms to find sparse Xi
epsguess(biginds,ind) = Theta(:,biginds)\X_dot(:,ind); %compute LS for non-small coeff
end
end

[time, sol] = ode15s(@SINDyfwd, [0:0.1:19], [x1_noisy(1), x2_noisy(1)], [], epsguess); %Run fwd model

%Plot (used to obtain figure 6)
figure(199)
plot(time, sol(:,1), 'k', 'LineWidth',2.0)
hold on
s = scatter(t_sampled, x1_noisy, 'k', 'filled');
s.Marker = 'hexagram';
s.SizeData = 150;
hold on
plot(time, sol(:,2), 'r', 'LineWidth',2.0)
hold on
s = scatter(t_sampled, x2_noisy, 'r', 'filled');
s.Marker = 'hexagram';
s.SizeData = 150;
hold off
xlabel('Time (Days)', 'fontsize',20) 
ylim([0 1500])
ylabel('cells/mm^2', 'fontsize',20)
set(gca, 'fontsize', 20)

epsguess % display resulting coefficients


%% Function for fwd model
function rhs = SINDyfwd(t, inits, epsguess)
M1 = inits(1); %initial condition for M1
M2 = inits(2); %initial condition for M2

%coefficient values
theta1 = epsguess(1,1);
theta2 = epsguess(2,1);
theta3 = epsguess(3,1);
theta4 = epsguess(4,1);
theta5 = epsguess(5,1);
theta6 = epsguess(6,1);
theta7 = epsguess(1, 2);
theta8 = epsguess(2, 2);
theta9 = epsguess(3, 2);
theta10 = epsguess(4, 2);
theta11 = epsguess(5,2); 
theta12 = epsguess(6,2);

%Equations
dM1 = theta1 + theta2*M1 + theta3*M2 + theta4*M1^2 + theta5*M2^2 +theta6*M1*M2;
dM2 = theta7 + theta8*M1 + theta9*M2 + theta10*M1^2 + theta11*M2^2 + theta12*M1*M2;

rhs = [dM1; dM2];
end

