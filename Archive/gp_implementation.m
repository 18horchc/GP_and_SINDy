clear; clc; close all;
%% Data
tpoints_M1 = [0,0,1,1,2,2,3,3,3,3,3,5,7,7,7,7,14,14,14]'; %time points where experimental data is observed
datapointsM1 = [0,10,5,50,120,125,375,62,102,100,60,325,600,55,900,225,750,400,1300]'; %M1 experimental data

tpoints_M2 = [0,0, 1, 1,1, 2, 2, 3,3,3,3,3, 5, 7, 7, 7, 7, 14,14,14]'; %time points where experimental data is observed
datapointsM2 = [0,10,170,15,50, 90,269,300,15,57,100,160, 800, 600,6,400,270,200,110,100]'; %M2 experimental data


%% Fitting GP
X_1 = tpoints_M1(:);
y_1 = datapointsM1(:);
y1_sqrt = sqrt(y_1 + 1e-6);  


X_2 = tpoints_M2(:);
y_2 = datapointsM2(:);
y2_sqrt = sqrt(y_2 + 1e-6); 


gprMdl_1 = fitrgp(X_1, y1_sqrt, ...
    'KernelFunction', 'squaredexponential', ...
    'BasisFunction', 'constant', ...
    'FitMethod', 'exact', ...
    'PredictMethod', 'exact', ...
    'Standardize', true);

gprMdl_2 = fitrgp(X_2, y2_sqrt, ...
    'KernelFunction', 'squaredexponential', ...
    'BasisFunction', 'constant', ...
    'FitMethod', 'exact', ...
    'PredictMethod', 'exact', ...
    'Standardize', true);


X1new = linspace(min(X_1), max(X_1), 300)';
[y1pred_sqrt, ~, y1int_sqrt] = predict(gprMdl_1, X1new);

X2new = linspace(min(X_2), max(X_2), 300)';
[y2pred_sqrt, ~, y2int_sqrt] = predict(gprMdl_2, X2new);


y1pred = max(0, y1pred_sqrt).^2;
y1low  = max(0, y1int_sqrt(:,1)).^2;
y1high = max(0, y1int_sqrt(:,2)).^2;

y2pred = max(0, y2pred_sqrt).^2;
y2low  = max(0, y2int_sqrt(:,1)).^2;
y2high = max(0, y2int_sqrt(:,2)).^2;


%% Plotting
figure; 
subplot(2,1,1);
hold on;
scatter(X_1, y_1, 50, 'r', 'filled');
plot(X1new, y1pred, 'b-', 'LineWidth', 2);
plot(X1new, y1low,  'k--', 'LineWidth', 1.4);
plot(X1new, y1high, 'k--', 'LineWidth', 1.4);
xlabel('Time (days)');
ylabel('M1 cell counts');
title('GP of M1');
legend('Observed','GP mean','95% prediction interval','Location','best');


subplot(2,1,2)
hold on;
scatter(X_2, y_2, 50, 'r', 'filled');
plot(X2new, y2pred, 'b-', 'LineWidth', 2);
plot(X2new, y2low,  'k--', 'LineWidth', 1.4);
plot(X2new, y2high, 'k--', 'LineWidth', 1.4);
xlabel('Time (days)');
ylabel('M2 cell counts');
title('GP of M2');
legend('Observed','GP mean','95% prediction interval','Location','best');
hold off;

