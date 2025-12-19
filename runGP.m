clear; clc; close all;

%% Data
t_1 = [0,0,1,1,2,2,3,3,3,3,3,5,7,7,7,7,14,14,14]';
y_1 = [0,10,5,50,120,125,375,62,102,100,60,325,600,55,900,225,750,400,1300]';

t_2 = [0,0, 1, 1,1, 2, 2, 3,3,3,3,3, 5, 7, 7, 7, 7, 14,14,14]';
y_2 = [0,10,170,15,50, 90,269,300,15,57,100,160, 800, 600,6,400,270,200,110,100]';

%% Run Analysis
figure;

% Call for dataset 1
subplot(2,1,1);
[mdl1, x1, p1] = fitAndPlotGP(t_1, y_1, 'GP of x', 'x', false);

% Call for dataset 2 
subplot(2,1,2);
[mdl2, x2, p2] = fitAndPlotGP(t_2, y_2, 'GP of y', 'y', false);


%% With derivative
% Example for M1 data
[mdl1, t_gp, x_gp, dx_gp] = fitAndPlotGP(t_1, y_1, 'M1 Analysis', 'Counts', true);

% Your SINDy "Data" matrix (X) and "Derivative" matrix (Xdot)
% If you have multiple species (M1 and M2), concatenate them side-by-side
X_sindy = x_gp;      % This is your denoised state
Xdot_sindy = dx_gp;  % This is your analytical derivative