% Authors: Sara Amato, Andrea Arnold <anarnold@wpi.edu>
% Last updated: August 29, 2025
%
% Corresponds with the manuscript:
% 
% S. Amato and A. Arnold (2025) Data-driven modeling and prediction of 
% microglial cell dynamics in the ischemic penumbra. 
% Preprint available on arXiv: https://arxiv.org/abs/2404.10915


%% MCMC for ODE model

% Note: This code requires use of MCMC Toolbox for MATLAB, available at:
% https://github.com/mjlaine/mcmcstat

clear; close all; clc;

tic %timer

%% Data

tpoints = [0, 1, 2, 3, 5, 7, 14]; %time points where experimental data is observed
datapointsM1 = [5, 27.5, 122.5, 139.8, 325, 445, 816.67]; %M1 experimental data
datapointsM2 = [5, 78.33, 179.5, 126.4, 800, 319, 136.67]; %M2 experimental data


data(:, 1) = datapointsM1;
data(:, 2) = datapointsM2;

data1.x = [0:0.1:14];
data1.y = data(:, :);


%% Specify Parameters

params = {
{'theta_4', normrnd(70.3415, 0.2*70.3415), -Inf, Inf, 70.3415, 0.2*70.3415}
{'theta_3', normrnd(0.3157, 0.2*abs(0.3157)), -Inf, Inf, 0.3157, 0.2*abs(0.3157)}
{'theta_5', normrnd(-0.1569, 0.2*abs(-0.1569)), -Inf, Inf, -0.1569, 0.2*abs(-0.1569)}
};

%% Options

model.ssfun = @model1av
model.N = length(data1.y);
model.S20 = [1];
model.N0 = [2];

options.updatesigma = 1;

%% Run toolbox

options.nsimu = 200000;
[results, chain, s2chain] = mcmcrun(model, data1, params, options);

%Run again starting at ending point of previous
options.nsimu = 200000;
[results, chain1, s2chain] = mcmcrun(model, data1, params, options, results);

toc

%% Results

chainstats(chain1, results) 
theta4mean = mean(chain1(:, 1)) %posterior means
theta3mean = mean(chain1(:, 2))
theta5mean = mean(chain1(:, 3))


figure(5); clf
mcmcplot(chain1,[],results,'pairs');
figure(6); clf
mcmcplot(chain1,[],results,'denspanel',2);

%% Forward Model

[time, sol] = ode15s(@inflammation, [0:0.1:14], [5,5], [], [theta4mean, theta3mean, theta5mean]);

%figure 8 rightmost plot
figure(12)
plot(time, sol(:,1), 'k', 'LineWidth',2.0)
hold on
s = scatter(tpoints, datapointsM1, 'k', 'filled',  'SizeData', 150);
s.Marker = 'hexagram';
s.SizeData = 150;
hold on
plot(time, sol(:,2), 'r', 'LineWidth',2.0)
hold on
s = scatter(tpoints, datapointsM2, 'r', 'filled',  'SizeData', 150);
s.Marker = 'hexagram';
s.SizeData = 150;
xlabel('Time (Days)', 'fontsize',20) 
ylim([0 1500])
ylabel('cells/mm^2', 'fontsize',20)
set(gca, 'fontsize', 20)

%% Prediction plots for UQ- Figure 10

N = 1000; %Number of times to sample
predtime = [0:0.1:50];
solpred = zeros(length(predtime), 2, N); 
for i = 1:N %Sample N times from chain and run fwd model
   m = randsample([1:1:N], 1, true);
    theta4pred = chain(m, 1);
    theta3pred = chain(m, 2);
    theta5pred = chain(m, 3);
    [time, solpred(:, :, i)] = ode15s(@inflammationfwd, [0:0.1:50], [5,5], [], [theta4pred, theta5pred, theta3pred]);
end

%get mean and std in each case
for j = 1:length(predtime)
    CIM1mean(j) = mean(solpred(j, 1, :));
    CIM1std(j) = std(solpred(j, 1, :));
    z = 2;
    sqrtn = 1;
    upperCIM1(j) = CIM1mean(j) + (z*(CIM1std(j)/sqrtn));
    lowerCIM1(j) = CIM1mean(j) - (z*(CIM1std(j)/sqrtn));
    
    CIM2mean(j) = mean(solpred(j, 2, :));
    CIM2std(j) = std(solpred(j, 2, :));
    upperCIM2(j) = CIM2mean(j) + z*(CIM2std(j)/sqrtn);
    lowerCIM2(j) = CIM2mean(j) - z*(CIM2std(j)/sqrtn);
end

figure(199)
x = predtime;
fill([x fliplr(x)],[upperCIM1 fliplr(lowerCIM1)],'k','EdgeColor','none','FaceAlpha',0.2);
hold on
fill([x fliplr(x)],[upperCIM2 fliplr(lowerCIM2)],'r','EdgeColor','none','FaceAlpha',0.2);
hold on
plot(predtime, CIM1mean, 'k', 'LineWidth',2.0)
hold on
s = scatter(tpoints, datapointsM1, 'k', 'filled');
s.Marker = 'hexagram';
s.SizeData = 150;
hold off
hold on
plot(predtime, CIM2mean, 'r', 'LineWidth',2.0)
hold on
s =  scatter(tpoints, datapointsM2, 'r', 'filled',  'SizeData', 150);
s.Marker = 'hexagram';
s.SizeData = 150;
hold off
hold on
scatter(35, 550,  'k', 'filled',  'SizeData', 150)
hold on
scatter(35, 50,  'r', 'filled',  'SizeData', 150)
xlabel('Time (Days)', 'fontsize',20) 
ylim([0 1500])
xlim([0 50])
ylabel('cells/mm^2', 'fontsize',20)
set(gca, 'fontsize', 20)


%% Functions

function ss = model1av(theta, data1) %%SS function used in MCMC
time = data1.x;
ydata = data1.y;
init = [5,5];
ymodel = inflammationfun(time, theta, init);
tpoints = [0, 1, 2, 3, 5, 7, 14];
for i = 1:length(time)
    if time(i)== tpoints(1)
        ycomp(1, :) = ymodel(1, :);
    elseif time(i) == tpoints(2)
        ycomp(2, :) = ymodel(11, :);
    elseif time(i) == tpoints(3)
        ycomp(3, :) = ymodel(21, :);
    elseif time(i) == tpoints(4)
        ycomp(4, :) = ymodel(31, :);
   elseif time(i) == tpoints(5)
        ycomp(5, :) = ymodel(51, :);
    elseif time(i) == tpoints(6)
        ycomp(6, :) = ymodel(71, :);
    elseif time(i) == tpoints(7)
        ycomp(7, :) = ymodel(141, :);
    end
end

ss = sum(sum((ycomp - (ydata(:, :))).^2));
end



function y = inflammationfun(time, theta, inits) %Model used in mcmc
[t, y] = ode15s(@inflammation, time, inits, [], theta);
end


function rhs = inflammation(t, inits, theta) %model used in MCMC
M1 = inits(1);
M2 = inits(2);

theta1 = 32.9650;
theta2 = -0.1377;
theta3 = theta(2);
theta4 = 0;
theta5 = 0;
theta6 = 0;
theta7 = theta(1);
theta8 = theta(3);
theta9 = 0.0217;
theta10 = 0;
theta11 = 0; 
theta12 = 0;

% Equations
dM1 = theta1 + theta2*M1 + theta3*M2 + theta4*M1^2 + theta5*M2^2 + theta6*M1*M2;
dM2 = theta7 + theta8*M1 + theta9*M2 + theta10*M1^2 + theta11*M2^2 + theta12*M1*M2;

rhs = [dM1; dM2];
end


function rhs = inflammationfwd(t, inits, theta) %FWD Model
M1 = inits(1);
M2 = inits(2);

theta1 = 32.9650;
theta2 = -0.1377;
theta3 = theta(3);
theta4 = 0;
theta5 = 0;
theta6 = 0;
theta7 = theta(1);
theta8 = theta(2);
theta9 = 0.0217;
theta10 = 0;
theta11 = 0; 
theta12 = 0;

% Equations
dM1 = theta1 + theta2*M1 + theta3*M2 + theta4*M1^2 + theta5*M2^2 + theta6*M1*M2;
dM2 = theta7 + theta8*M1 + theta9*M2 + theta10*M1^2 + theta11*M2^2 + theta12*M1*M2;

rhs = [dM1; dM2];
end

