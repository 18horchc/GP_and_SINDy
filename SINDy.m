% Authors: Sara Amato, Andrea Arnold <anarnold@wpi.edu>
% Last updated: August 29, 2025
%
% Corresponds with the manuscript:
% 
% S. Amato and A. Arnold (2025) Data-driven modeling and prediction of 
% microglial cell dynamics in the ischemic penumbra. 
% Preprint available on arXiv: https://arxiv.org/abs/2404.10915


clear; close all; clc

%% Data 

newtime = [0, 1, 2, 3, 7, 14]; %These are the points where we observe experimental data
datapointsM2 = [5, 78.33, 179.5, 126.4, 319, 136.67]; %Experimental measurements of M2
datapointsM1 = [5, 27.5, 122.5, 139.8, 445, 816.67]; %Experimental measurements of M1

%% Obtain Derivative Information 

t = [0:0.1:14]; %time that we want to fit the polynomial over

pM1four = polyfit(newtime, datapointsM1, 4); %obtain coefficients for a fourth order polynomial fit to M1 data
polyM1four = polyval(pM1four, t); %Compute y values of polynomial for t

%analytically compute derivative for all values of t
for i = 1:length(t)
M1estderiv(i) = 0 + pM1four(4) + 2*pM1four(3)*t(i) + 3*pM1four(2)*t(i)^2 + 4*pM1four(1)*t(i)^3; %
end

%Extract derivative information for time points where experimental data is
%present
derivM1 = [M1estderiv(1), M1estderiv(11), M1estderiv(21), M1estderiv(31), M1estderiv(71), M1estderiv(141)];

pM2quad = polyfit(newtime, datapointsM2, 2); %obtain coefficients for a second order polynomial fit to M2 data
polyM2quad = polyval(pM2quad, t); %Compute y values of polynomial for t

%analytically compute derivative for all values of t
for i = 1:length(t)
M2estderiv(i) = 0 + pM2quad(2) + 2*pM2quad(1)*t(i);%
end

%Extract derivative information for time points where experimental data is
%present
derivM2 = [M2estderiv(1), M2estderiv(11), M2estderiv(21), M2estderiv(31), M2estderiv(71), M2estderiv(141)];


%% SINDy

X_dot = [derivM1; derivM2]'; %Derivative information; rhs of 2.1

Theta = [ones(length(datapointsM1),1)'; datapointsM1(1:end); datapointsM2(1:end); datapointsM1(1:end).^2;...
   datapointsM2(1:end).^2; datapointsM1(1:end).*datapointsM2(1:end)]'; %candidate function matrix

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

[time, sol] = ode15s(@SINDyfwd, [0:0.1:14], [datapointsM1(1), datapointsM2(1)], [], epsguess); %Run fwd model

%Plot (used to obtain figure 6)
figure(199)
plot(time, sol(:,1), 'k', 'LineWidth',2.0)
hold on
s = scatter(newtime, datapointsM1, 'k', 'filled');
s.Marker = 'hexagram';
s.SizeData = 150;
hold on
plot(time, sol(:,2), 'r', 'LineWidth',2.0)
hold on
s = scatter(newtime, datapointsM2, 'r', 'filled');
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

