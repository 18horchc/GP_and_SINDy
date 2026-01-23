%% 
% I think microglia data most likely follow a log growth curve. I want to model 
% my populations with this in mind and compare kernels.
% 
% I also want to think about M1 and M2 together as predator-prey
% 
% 
% 
% Questions
%% 
% * What does the optimizer options in fitgrp() do? Note this is different than 
% the hyperparameter optimization
%% Data

clear; clc; close all; 
%samples
tpoints_M1 = [0,0,1,1,2,2,3,3,3,3,3,5,7,7,7,7,14,14,14]'; %time points where experimental data is observed
datapointsM1 = [0,10,5,50,120,125,375,62,102,100,60,325,600,55,900,225,750,400,1300]'; %M1 experimental data

tpoints_M2 = [0,0, 1, 1,1, 2, 2, 3,3,3,3,3, 5, 7, 7, 7, 7, 14,14,14]'; %time points where experimental data is observed
datapointsM2 = [0,10,170,15,50, 90,269,300,15,57,100,160, 800, 600,6,400,270,200,110,100]'; %M2 experimental data


%for testing
xp = linspace(0, 14, 300)'; 

%% Pre-processing
% Log transform to ensure positive: y_log = log(y_raw + 1)
% The +1 handles zero counts
datapointsM1_log = log(datapointsM1 + 1);

% sqrt transform to ensure positive (for later use if needed)
% Square root transform: y_sqrt = sqrt(y_raw)
% No +1 needed since sqrt(0) = 0
datapointsM1_sqrt = sqrt(datapointsM1);

% Standardize so I can use mean 0 (try with and try without)
% Calculate mean and standard deviation of log-transformed data
mu_train_log = mean(datapointsM1_log);
sigma_train_log = std(datapointsM1_log);

% Z-score standardization: z = (y_log - mu_train) / sigma_train
datapointsM1_standardized_log = (datapointsM1_log - mu_train_log) / sigma_train_log;

% Calculate mean and standard deviation of sqrt-transformed data
mu_train_sqrt = mean(datapointsM1_sqrt);
sigma_train_sqrt = std(datapointsM1_sqrt);

% Z-score standardization for sqrt data: z = (y_sqrt - mu_train) / sigma_train
datapointsM1_standardized_sqrt = (datapointsM1_sqrt - mu_train_sqrt) / sigma_train_sqrt;

fprintf('Original M1 data: min = %.2f, max = %.2f, mean = %.2f, std = %.2f\n', ...
    min(datapointsM1), max(datapointsM1), mean(datapointsM1), std(datapointsM1));
fprintf('Log-transformed M1 data: min = %.4f, max = %.4f, mean = %.4f, std = %.4f\n', ...
    min(datapointsM1_log), max(datapointsM1_log), mean(datapointsM1_log), std(datapointsM1_log));
fprintf('Standardized (log) M1 data: min = %.4f, max = %.4f, mean = %.6f, std = %.6f\n', ...
    min(datapointsM1_standardized_log), max(datapointsM1_standardized_log), ...
    mean(datapointsM1_standardized_log), std(datapointsM1_standardized_log));
fprintf('(Note: Mean should be ~0 and std should be ~1)\n\n');
fprintf('Sqrt-transformed M1 data: min = %.4f, max = %.4f, mean = %.4f, std = %.4f\n', ...
    min(datapointsM1_sqrt), max(datapointsM1_sqrt), mean(datapointsM1_sqrt), std(datapointsM1_sqrt));
fprintf('Standardized (sqrt) M1 data: min = %.4f, max = %.4f, mean = %.6f, std = %.6f\n', ...
    min(datapointsM1_standardized_sqrt), max(datapointsM1_standardized_sqrt), ...
    mean(datapointsM1_standardized_sqrt), std(datapointsM1_standardized_sqrt));
fprintf('(Note: Mean should be ~0 and std should be ~1)\n\n');






%% *Comparing kernels* 
% *Without Bayesian optimizer*
% More info about kernels: <https://www.mathworks.com/help/stats/kernel-covariance-function-options.html 
% https://www.mathworks.com/help/stats/kernel-covariance-function-options.html> 
% 
% Not setting any inital values for kernel parameters
% 
% Using default constant basis function for all kernels
%
% Other settings:
% * Using FitMethod = exact (default) _(Exact Gaussian process regression. This 
% value is the default if n ≤ 2000, where n is the number of observations.)_
% * Always using predict = exact _(Exact Gaussian process regression method. 
% This value is the default if n ≤ 10,000.)_
% * Standardize = false (default) - Data is already log-transformed and standardized,
%   so we don't want fitrgp() to standardize it again

%Squared exponental (hyperparams: lengthsale (l))
% Fit GP on standardized log-transformed data (z)
gprMdl_SE_M1_log = fitrgp(tpoints_M1, datapointsM1_standardized_log, ...
    'KernelFunction', 'squaredexponential', 'Standardize', false);
[ypred_SE_M1_log, std_SE_M1_log] = predict(gprMdl_SE_M1_log,xp);


%Exponential (aka Matern 1/2)
gprMdl_exp_M1_log = fitrgp(tpoints_M1, datapointsM1_standardized_log, ...
    'KernelFunction', 'exponential', 'Standardize', false);
[ypred_exp_M1_log, std_exp_M1_log] = predict(gprMdl_exp_M1_log,xp);


%Matern 3/2 [v=3/2] (hyperparams: v, l)
gprMdl_M32_M1_log = fitrgp(tpoints_M1, datapointsM1_standardized_log, ...
    'KernelFunction', 'matern32', 'Standardize', false);
[ypred_M32_M1_log, std_M32_M1_log] = predict(gprMdl_M32_M1_log,xp);

%Matern 5/2 [v=5/2] (hyperparams: v, l)
gprMdl_M52_M1_log = fitrgp(tpoints_M1, datapointsM1_standardized_log, ...
    'KernelFunction', 'matern52', 'Standardize', false);
[ypred_M52_M1_log, std_M52_M1_log] = predict(gprMdl_M52_M1_log,xp);

%Rational Quadratic (hyperparams: alpha, l)
gprMdl_RQ_M1_log = fitrgp(tpoints_M1, datapointsM1_standardized_log, ...
    'KernelFunction', 'rationalquadratic', 'Standardize', false);
[ypred_RQ_M1_log, std_RQ_M1_log] = predict(gprMdl_RQ_M1_log,xp);

%ARD = automatic relevance determination
%ARD squared exponential
gprMdl_ardSE_M1_log = fitrgp(tpoints_M1, datapointsM1_standardized_log, ...
    'KernelFunction', 'ardsquaredexponential', 'Standardize', false);
[ypred_ardSE_M1_log, std_ardSE_M1_log] = predict(gprMdl_ardSE_M1_log,xp);

%ARD exponential
gprMdl_ardExp_M1_log = fitrgp(tpoints_M1, datapointsM1_standardized_log, ...
    'KernelFunction', 'ardexponential', 'Standardize', false);
[ypred_ardExp_M1_log, std_ardExp_M1_log] = predict(gprMdl_ardExp_M1_log,xp);

%ARD Matern 3/2
gprMdl_ardM32_M1_log = fitrgp(tpoints_M1, datapointsM1_standardized_log, ...
    'KernelFunction', 'ardmatern32', 'Standardize', false);
[ypred_ardM32_M1_log, std_ardM32_M1_log] = predict(gprMdl_ardM32_M1_log,xp);


%ARD Matern 5/2
gprMdl_ardM52_M1_log = fitrgp(tpoints_M1, datapointsM1_standardized_log, ...
    'KernelFunction', 'ardmatern52', 'Standardize', false);
[ypred_ardM52_M1_log, std_ardM52_M1_log] = predict(gprMdl_ardM52_M1_log,xp);

%ARD rational quadratic
gprMdl_ardRQ_M1_log = fitrgp(tpoints_M1, datapointsM1_standardized_log, ...
    'KernelFunction', 'ardrationalquadratic', 'Standardize', false);
[ypred_ardRQ_M1_log, std_ardRQ_M1_log] = predict(gprMdl_ardRQ_M1_log,xp);

%% GP Models with Square Root Transform
% Fit GP models on standardized sqrt-transformed data

%Squared exponental (hyperparams: lengthsale (l))
% Fit GP on standardized sqrt-transformed data (z)
gprMdl_SE_M1_sqrt = fitrgp(tpoints_M1, datapointsM1_standardized_sqrt, ...
    'KernelFunction', 'squaredexponential', 'Standardize', false);
[ypred_SE_M1_sqrt, std_SE_M1_sqrt] = predict(gprMdl_SE_M1_sqrt,xp);

%Exponential (aka Matern 1/2)
gprMdl_exp_M1_sqrt = fitrgp(tpoints_M1, datapointsM1_standardized_sqrt, ...
    'KernelFunction', 'exponential', 'Standardize', false);
[ypred_exp_M1_sqrt, std_exp_M1_sqrt] = predict(gprMdl_exp_M1_sqrt,xp);

%Matern 3/2 [v=3/2] (hyperparams: v, l)
gprMdl_M32_M1_sqrt = fitrgp(tpoints_M1, datapointsM1_standardized_sqrt, ...
    'KernelFunction', 'matern32', 'Standardize', false);
[ypred_M32_M1_sqrt, std_M32_M1_sqrt] = predict(gprMdl_M32_M1_sqrt,xp);

%Matern 5/2 [v=5/2] (hyperparams: v, l)
gprMdl_M52_M1_sqrt = fitrgp(tpoints_M1, datapointsM1_standardized_sqrt, ...
    'KernelFunction', 'matern52', 'Standardize', false);
[ypred_M52_M1_sqrt, std_M52_M1_sqrt] = predict(gprMdl_M52_M1_sqrt,xp);

%Rational Quadratic (hyperparams: alpha, l)
gprMdl_RQ_M1_sqrt = fitrgp(tpoints_M1, datapointsM1_standardized_sqrt, ...
    'KernelFunction', 'rationalquadratic', 'Standardize', false);
[ypred_RQ_M1_sqrt, std_RQ_M1_sqrt] = predict(gprMdl_RQ_M1_sqrt,xp);

%ARD = automatic relevance determination
%ARD squared exponential
gprMdl_ardSE_M1_sqrt = fitrgp(tpoints_M1, datapointsM1_standardized_sqrt, ...
    'KernelFunction', 'ardsquaredexponential', 'Standardize', false);
[ypred_ardSE_M1_sqrt, std_ardSE_M1_sqrt] = predict(gprMdl_ardSE_M1_sqrt,xp);

%ARD exponential
gprMdl_ardExp_M1_sqrt = fitrgp(tpoints_M1, datapointsM1_standardized_sqrt, ...
    'KernelFunction', 'ardexponential', 'Standardize', false);
[ypred_ardExp_M1_sqrt, std_ardExp_M1_sqrt] = predict(gprMdl_ardExp_M1_sqrt,xp);

%ARD Matern 3/2
gprMdl_ardM32_M1_sqrt = fitrgp(tpoints_M1, datapointsM1_standardized_sqrt, ...
    'KernelFunction', 'ardmatern32', 'Standardize', false);
[ypred_ardM32_M1_sqrt, std_ardM32_M1_sqrt] = predict(gprMdl_ardM32_M1_sqrt,xp);

%ARD Matern 5/2
gprMdl_ardM52_M1_sqrt = fitrgp(tpoints_M1, datapointsM1_standardized_sqrt, ...
    'KernelFunction', 'ardmatern52', 'Standardize', false);
[ypred_ardM52_M1_sqrt, std_ardM52_M1_sqrt] = predict(gprMdl_ardM52_M1_sqrt,xp);

%ARD rational quadratic
gprMdl_ardRQ_M1_sqrt = fitrgp(tpoints_M1, datapointsM1_standardized_sqrt, ...
    'KernelFunction', 'ardrationalquadratic', 'Standardize', false);
[ypred_ardRQ_M1_sqrt, std_ardRQ_M1_sqrt] = predict(gprMdl_ardRQ_M1_sqrt,xp);

%%To explore later %%
%non-stationary kernels
%combining kernels through sums or products
%% Reverse Transformations: Convert predictions back to original count space
% Step 1: Reverse Z-Score (Un-standardize) to get back to log-space
% μ_log = (μ_z × σ_train) + μ_train
% σ²_log = σ²_z × σ²_train (note: std_SE_M1_log is σ_z, so σ²_z = std_SE_M1_log²)

% Helper function to reverse transform predictions (for log transform)
% Input: ypred_z (mean in standardized space), std_z (std in standardized space)
% Output: struct with predictions in count space
reverseTransformLog = @(ypred_z, std_z) struct(...
    'mu_log', ypred_z * sigma_train_log + mu_train_log, ...
    'sigma_log', std_z * sigma_train_log, ...
    'mu_count_mean', exp(ypred_z * sigma_train_log + mu_train_log + 0.5 * (std_z * sigma_train_log).^2) - 1, ...
    'mu_count_median', exp(ypred_z * sigma_train_log + mu_train_log) - 1, ...
    'lower_log', (ypred_z - 1.96 * std_z) * sigma_train_log + mu_train_log, ...
    'upper_log', (ypred_z + 1.96 * std_z) * sigma_train_log + mu_train_log, ...
    'lower_count', exp((ypred_z - 1.96 * std_z) * sigma_train_log + mu_train_log) - 1, ...
    'upper_count', exp((ypred_z + 1.96 * std_z) * sigma_train_log + mu_train_log) - 1);

% Helper function to reverse transform predictions (for sqrt transform)
% Input: ypred_z (mean in standardized space), std_z (std in standardized space)
% Output: struct with predictions in count space
reverseTransformSqrt = @(ypred_z, std_z) struct(...
    'mu_sqrt', ypred_z * sigma_train_sqrt + mu_train_sqrt, ...
    'sigma_sqrt', std_z * sigma_train_sqrt, ...
    'mu_count_mean', (ypred_z * sigma_train_sqrt + mu_train_sqrt).^2 + (std_z * sigma_train_sqrt).^2, ...
    'mu_count_median', (ypred_z * sigma_train_sqrt + mu_train_sqrt).^2, ...
    'lower_sqrt', (ypred_z - 1.96 * std_z) * sigma_train_sqrt + mu_train_sqrt, ...
    'upper_sqrt', (ypred_z + 1.96 * std_z) * sigma_train_sqrt + mu_train_sqrt, ...
    'lower_count', max(0, (ypred_z - 1.96 * std_z) * sigma_train_sqrt + mu_train_sqrt).^2, ...
    'upper_count', ((ypred_z + 1.96 * std_z) * sigma_train_sqrt + mu_train_sqrt).^2);

% Transform all predictions back to count space (log transform)
pred_SE_log = reverseTransformLog(ypred_SE_M1_log, std_SE_M1_log);
pred_exp_log = reverseTransformLog(ypred_exp_M1_log, std_exp_M1_log);
pred_M32_log = reverseTransformLog(ypred_M32_M1_log, std_M32_M1_log);
pred_M52_log = reverseTransformLog(ypred_M52_M1_log, std_M52_M1_log);
pred_RQ_log = reverseTransformLog(ypred_RQ_M1_log, std_RQ_M1_log);
pred_ardSE_log = reverseTransformLog(ypred_ardSE_M1_log, std_ardSE_M1_log);
pred_ardExp_log = reverseTransformLog(ypred_ardExp_M1_log, std_ardExp_M1_log);
pred_ardM32_log = reverseTransformLog(ypred_ardM32_M1_log, std_ardM32_M1_log);
pred_ardM52_log = reverseTransformLog(ypred_ardM52_M1_log, std_ardM52_M1_log);
pred_ardRQ_log = reverseTransformLog(ypred_ardRQ_M1_log, std_ardRQ_M1_log);

% Transform all predictions back to count space (sqrt transform)
pred_SE_sqrt = reverseTransformSqrt(ypred_SE_M1_sqrt, std_SE_M1_sqrt);
pred_exp_sqrt = reverseTransformSqrt(ypred_exp_M1_sqrt, std_exp_M1_sqrt);
pred_M32_sqrt = reverseTransformSqrt(ypred_M32_M1_sqrt, std_M32_M1_sqrt);
pred_M52_sqrt = reverseTransformSqrt(ypred_M52_M1_sqrt, std_M52_M1_sqrt);
pred_RQ_sqrt = reverseTransformSqrt(ypred_RQ_M1_sqrt, std_RQ_M1_sqrt);
pred_ardSE_sqrt = reverseTransformSqrt(ypred_ardSE_M1_sqrt, std_ardSE_M1_sqrt);
pred_ardExp_sqrt = reverseTransformSqrt(ypred_ardExp_M1_sqrt, std_ardExp_M1_sqrt);
pred_ardM32_sqrt = reverseTransformSqrt(ypred_ardM32_M1_sqrt, std_ardM32_M1_sqrt);
pred_ardM52_sqrt = reverseTransformSqrt(ypred_ardM52_M1_sqrt, std_ardM52_M1_sqrt);
pred_ardRQ_sqrt = reverseTransformSqrt(ypred_ardRQ_M1_sqrt, std_ardRQ_M1_sqrt);


% 
%%  Plot all GP models

% Squared Exponential
% Plot original count data and reverse-transformed predictions
figure;
% Top subplot: Log transform
subplot(2,1,1);
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_SE_log.mu_count_median, 'g', 'LineWidth', 1.5);  % Median (most likely path)
plot(xp, pred_SE_log.lower_count, 'g--', 'LineWidth', 1);     % Lower 95% CI
plot(xp, pred_SE_log.upper_count, 'g--', 'LineWidth', 1);      % Upper 95% CI
xlabel('time');
ylabel('M1 cell count');
title('GPR with Squared Exponential Kernel (Log Transform)');
legend('Data', 'GPR predictions (SE kernel, median)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;
% Bottom subplot: Sqrt transform
subplot(2,1,2);
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_SE_sqrt.mu_count_mean, 'r', 'LineWidth', 1.5);  % Mean (corrected)
plot(xp, pred_SE_sqrt.lower_count, 'r--', 'LineWidth', 1);     % Lower 95% CI
plot(xp, pred_SE_sqrt.upper_count, 'r--', 'LineWidth', 1);      % Upper 95% CI
xlabel('time');
ylabel('M1 cell count');
title('GPR with Squared Exponential Kernel (Sqrt Transform)');
legend('Data', 'GPR predictions (SE kernel, mean)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;

% Exponential
figure;
% Top subplot: Log transform
subplot(2,1,1);
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_exp_log.mu_count_median, 'g', 'LineWidth', 1.5);
plot(xp, pred_exp_log.lower_count, 'g--', 'LineWidth', 1);
plot(xp, pred_exp_log.upper_count, 'g--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with Exponential Kernel (Log Transform)');
legend('Data', 'GPR predictions (Exponential kernel, median)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;
% Bottom subplot: Sqrt transform
subplot(2,1,2);
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_exp_sqrt.mu_count_mean, 'r', 'LineWidth', 1.5);
plot(xp, pred_exp_sqrt.lower_count, 'r--', 'LineWidth', 1);
plot(xp, pred_exp_sqrt.upper_count, 'r--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with Exponential Kernel (Sqrt Transform)');
legend('Data', 'GPR predictions (Exponential kernel, mean)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;

% Matern 3/2
figure;
% Top subplot: Log transform
subplot(2,1,1);
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_M32_log.mu_count_median, 'g', 'LineWidth', 1.5);
plot(xp, pred_M32_log.lower_count, 'g--', 'LineWidth', 1);
plot(xp, pred_M32_log.upper_count, 'g--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with Matern 3/2 Kernel (Log Transform)');
legend('Data', 'GPR predictions (Matern 3/2 kernel, median)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;
% Bottom subplot: Sqrt transform
subplot(2,1,2);
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_M32_sqrt.mu_count_mean, 'r', 'LineWidth', 1.5);
plot(xp, pred_M32_sqrt.lower_count, 'r--', 'LineWidth', 1);
plot(xp, pred_M32_sqrt.upper_count, 'r--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with Matern 3/2 Kernel (Sqrt Transform)');
legend('Data', 'GPR predictions (Matern 3/2 kernel, mean)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;

% Matern 5/2
figure;
% Top subplot: Log transform
subplot(2,1,1);
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_M52_log.mu_count_median, 'g', 'LineWidth', 1.5);
plot(xp, pred_M52_log.lower_count, 'g--', 'LineWidth', 1);
plot(xp, pred_M52_log.upper_count, 'g--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with Matern 5/2 Kernel (Log Transform)');
legend('Data', 'GPR predictions (Matern 5/2 kernel, median)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;
% Bottom subplot: Sqrt transform
subplot(2,1,2);
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_M52_sqrt.mu_count_mean, 'r', 'LineWidth', 1.5);
plot(xp, pred_M52_sqrt.lower_count, 'r--', 'LineWidth', 1);
plot(xp, pred_M52_sqrt.upper_count, 'r--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with Matern 5/2 Kernel (Sqrt Transform)');
legend('Data', 'GPR predictions (Matern 5/2 kernel, mean)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;

% Rational Quadratic
figure;
% Top subplot: Log transform
subplot(2,1,1);
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_RQ_log.mu_count_median, 'g', 'LineWidth', 1.5);
plot(xp, pred_RQ_log.lower_count, 'g--', 'LineWidth', 1);
plot(xp, pred_RQ_log.upper_count, 'g--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with Rational Quadratic Kernel (Log Transform)');
legend('Data', 'GPR predictions (Rational Quadratic kernel, median)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;
% Bottom subplot: Sqrt transform
subplot(2,1,2);
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_RQ_sqrt.mu_count_mean, 'r', 'LineWidth', 1.5);
plot(xp, pred_RQ_sqrt.lower_count, 'r--', 'LineWidth', 1);
plot(xp, pred_RQ_sqrt.upper_count, 'r--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with Rational Quadratic Kernel (Sqrt Transform)');
legend('Data', 'GPR predictions (Rational Quadratic kernel, mean)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;

% ARD Squared Exponential
figure;
% Top subplot: Log transform
subplot(2,1,1);
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_ardSE_log.mu_count_median, 'g', 'LineWidth', 1.5);
plot(xp, pred_ardSE_log.lower_count, 'g--', 'LineWidth', 1);
plot(xp, pred_ardSE_log.upper_count, 'g--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with ARD Squared Exponential Kernel (Log Transform)');
legend('Data', 'GPR predictions (ARD SE kernel, median)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;
% Bottom subplot: Sqrt transform
subplot(2,1,2);
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_ardSE_sqrt.mu_count_mean, 'r', 'LineWidth', 1.5);
plot(xp, pred_ardSE_sqrt.lower_count, 'r--', 'LineWidth', 1);
plot(xp, pred_ardSE_sqrt.upper_count, 'r--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with ARD Squared Exponential Kernel (Sqrt Transform)');
legend('Data', 'GPR predictions (ARD SE kernel, mean)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;

% ARD Exponential
figure;
% Top subplot: Log transform
subplot(2,1,1);
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_ardExp_log.mu_count_median, 'g', 'LineWidth', 1.5);
plot(xp, pred_ardExp_log.lower_count, 'g--', 'LineWidth', 1);
plot(xp, pred_ardExp_log.upper_count, 'g--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with ARD Exponential Kernel (Log Transform)');
legend('Data', 'GPR predictions (ARD Exponential kernel, median)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;
% Bottom subplot: Sqrt transform
subplot(2,1,2);
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_ardExp_sqrt.mu_count_mean, 'r', 'LineWidth', 1.5);
plot(xp, pred_ardExp_sqrt.lower_count, 'r--', 'LineWidth', 1);
plot(xp, pred_ardExp_sqrt.upper_count, 'r--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with ARD Exponential Kernel (Sqrt Transform)');
legend('Data', 'GPR predictions (ARD Exponential kernel, mean)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;

% ARD Matern 3/2
figure;
% Top subplot: Log transform
subplot(2,1,1);
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_ardM32_log.mu_count_median, 'g', 'LineWidth', 1.5);
plot(xp, pred_ardM32_log.lower_count, 'g--', 'LineWidth', 1);
plot(xp, pred_ardM32_log.upper_count, 'g--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with ARD Matern 3/2 Kernel (Log Transform)');
legend('Data', 'GPR predictions (ARD Matern 3/2 kernel, median)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;
% Bottom subplot: Sqrt transform
subplot(2,1,2);
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_ardM32_sqrt.mu_count_mean, 'r', 'LineWidth', 1.5);
plot(xp, pred_ardM32_sqrt.lower_count, 'r--', 'LineWidth', 1);
plot(xp, pred_ardM32_sqrt.upper_count, 'r--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with ARD Matern 3/2 Kernel (Sqrt Transform)');
legend('Data', 'GPR predictions (ARD Matern 3/2 kernel, mean)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;

% ARD Matern 5/2
figure;
% Top subplot: Log transform
subplot(2,1,1);
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_ardM52_log.mu_count_median, 'g', 'LineWidth', 1.5);
plot(xp, pred_ardM52_log.lower_count, 'g--', 'LineWidth', 1);
plot(xp, pred_ardM52_log.upper_count, 'g--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with ARD Matern 5/2 Kernel (Log Transform)');
legend('Data', 'GPR predictions (ARD Matern 5/2 kernel, median)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;
% Bottom subplot: Sqrt transform
subplot(2,1,2);
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_ardM52_sqrt.mu_count_mean, 'r', 'LineWidth', 1.5);
plot(xp, pred_ardM52_sqrt.lower_count, 'r--', 'LineWidth', 1);
plot(xp, pred_ardM52_sqrt.upper_count, 'r--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with ARD Matern 5/2 Kernel (Sqrt Transform)');
legend('Data', 'GPR predictions (ARD Matern 5/2 kernel, mean)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;

% ARD Rational Quadratic
figure;
% Top subplot: Log transform
subplot(2,1,1);
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_ardRQ_log.mu_count_median, 'g', 'LineWidth', 1.5);
plot(xp, pred_ardRQ_log.lower_count, 'g--', 'LineWidth', 1);
plot(xp, pred_ardRQ_log.upper_count, 'g--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with ARD Rational Quadratic Kernel (Log Transform)');
legend('Data', 'GPR predictions (ARD Rational Quadratic kernel, median)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;
% Bottom subplot: Sqrt transform
subplot(2,1,2);
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_ardRQ_sqrt.mu_count_mean, 'r', 'LineWidth', 1.5);
plot(xp, pred_ardRQ_sqrt.lower_count, 'r--', 'LineWidth', 1);
plot(xp, pred_ardRQ_sqrt.upper_count, 'r--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with ARD Rational Quadratic Kernel (Sqrt Transform)');
legend('Data', 'GPR predictions (ARD Rational Quadratic kernel, mean)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on; 


%% Metrics Table for Log Transform

% Compute marginal log liklihood

%I want to create a table where I track what was used to generate the
%model, and then put the calulated log liklihood


% I want to create a data set where I say which kernel was used, and then I
% calculate the following metrics (logliklihood, RMSE, and R^2) for each.
%so the columns will be: kernel (string), LogLikelihood (double), RMSE
%(double), and R2 (double)

%additionally, i'd ideally like to set up the below code in some sort of
%loop or other efficient codeblock to be able to calculate or find the
%below metrics for each gprmdl that I made above, and then put the relevant
%information in the table. 



% Store all models and their names (log transform)
models_log = {
    gprMdl_SE_M1_log, 'Squared Exponential';
    gprMdl_exp_M1_log, 'Exponential';
    gprMdl_M32_M1_log, 'Matern 3/2';
    gprMdl_M52_M1_log, 'Matern 5/2';
    gprMdl_RQ_M1_log, 'Rational Quadratic';
    gprMdl_ardSE_M1_log, 'ARD Squared Exponential';
    gprMdl_ardExp_M1_log, 'ARD Exponential';
    gprMdl_ardM32_M1_log, 'ARD Matern 3/2';
    gprMdl_ardM52_M1_log, 'ARD Matern 5/2';
    gprMdl_ardRQ_M1_log, 'ARD Rational Quadratic'
};

% Initialize arrays to store metrics
numModels = size(models_log, 1);
kernelNames_log = cell(numModels, 1);
logLikelihoods_log = zeros(numModels, 1);
RMSEs_log = zeros(numModels, 1);
R2s_log = zeros(numModels, 1);

% Calculate metrics for each model (log transform)
for i = 1:numModels
    model = models_log{i, 1};
    kernelNames_log{i} = models_log{i, 2};
    
    % Get log likelihood (property of the model)
    logLikelihoods_log(i) = model.LogLikelihood;
    
    % Get predictions on training data
    yPredict = resubPredict(model);
    yActual = datapointsM1_standardized_log;  % Actual training data (log-transformed and standardized)
    
    % Calculate residuals
    residuals = yActual - yPredict;
    
    % Calculate RMSE
    RMSEs_log(i) = sqrt(mean(residuals.^2));
    
    % Calculate R²
    SSR = sum((yPredict - mean(yActual)).^2);  % Regression sum of squares
    SST = sum((yActual - mean(yActual)).^2);   % Total sum of squares
    R2s_log(i) = SSR / SST;
end

% Create table with results (log transform)
resultsTable_log = table(kernelNames_log, logLikelihoods_log, RMSEs_log, R2s_log, ...
    'VariableNames', {'Kernel', 'LogLikelihood', 'RMSE', 'R2'});

% Display the table
disp('GPR Model Comparison Results for M1 Data (Log Transform):');
disp(resultsTable_log);

% Optionally sort by a metric (e.g., by LogLikelihood descending)
resultsTableSorted_log = sortrows(resultsTable_log, 'LogLikelihood', 'descend');
disp('Sorted by LogLikelihood (best to worst):');
disp(resultsTableSorted_log);

%% Metrics Table for Sqrt Transform
% Store all models and their names (sqrt transform)
models_sqrt = {
    gprMdl_SE_M1_sqrt, 'Squared Exponential';
    gprMdl_exp_M1_sqrt, 'Exponential';
    gprMdl_M32_M1_sqrt, 'Matern 3/2';
    gprMdl_M52_M1_sqrt, 'Matern 5/2';
    gprMdl_RQ_M1_sqrt, 'Rational Quadratic';
    gprMdl_ardSE_M1_sqrt, 'ARD Squared Exponential';
    gprMdl_ardExp_M1_sqrt, 'ARD Exponential';
    gprMdl_ardM32_M1_sqrt, 'ARD Matern 3/2';
    gprMdl_ardM52_M1_sqrt, 'ARD Matern 5/2';
    gprMdl_ardRQ_M1_sqrt, 'ARD Rational Quadratic'
};

% Initialize arrays to store metrics
kernelNames_sqrt = cell(numModels, 1);
logLikelihoods_sqrt = zeros(numModels, 1);
RMSEs_sqrt = zeros(numModels, 1);
R2s_sqrt = zeros(numModels, 1);

% Calculate metrics for each model (sqrt transform)
for i = 1:numModels
    model = models_sqrt{i, 1};
    kernelNames_sqrt{i} = models_sqrt{i, 2};
    
    % Get log likelihood (property of the model)
    logLikelihoods_sqrt(i) = model.LogLikelihood;
    
    % Get predictions on training data
    yPredict = resubPredict(model);
    yActual = datapointsM1_standardized_sqrt;  % Actual training data (sqrt-transformed and standardized)
    
    % Calculate residuals
    residuals = yActual - yPredict;
    
    % Calculate RMSE
    RMSEs_sqrt(i) = sqrt(mean(residuals.^2));
    
    % Calculate R²
    SSR = sum((yPredict - mean(yActual)).^2);  % Regression sum of squares
    SST = sum((yActual - mean(yActual)).^2);   % Total sum of squares
    R2s_sqrt(i) = SSR / SST;
end

% Create table with results (sqrt transform)
resultsTable_sqrt = table(kernelNames_sqrt, logLikelihoods_sqrt, RMSEs_sqrt, R2s_sqrt, ...
    'VariableNames', {'Kernel', 'LogLikelihood', 'RMSE', 'R2'});

% Display the table
disp('GPR Model Comparison Results for M1 Data (Sqrt Transform):');
disp(resultsTable_sqrt);

% Optionally sort by a metric (e.g., by LogLikelihood descending)
resultsTableSorted_sqrt = sortrows(resultsTable_sqrt, 'LogLikelihood', 'descend');
disp('Sorted by LogLikelihood (best to worst):');
disp(resultsTableSorted_sqrt);
%% 
% 
% *With Hyperparameter optimization*
% Parameters to optimize, specified as one of the following:
%% 
% * |'none'| — Do not optimize.
% * |'auto'| — Use |{'Sigma','Standardize'}|.
% * |'all'| — Optimize all eligible parameters, equivalent to |{'BasisFunction','KernelFunction','KernelScale','Sigma','Standardize'}|.
% * String array or cell array of eligible parameter names.
% * Vector of |optimizableVariable| objects, typically the output of <https://www.mathworks.com/help/stats/hyperparameters.html 
% |hyperparameters|>.
%% 
% The optimization attempts to minimize the cross-validation loss (error) for 
% |fitrgp| by varying the parameters. To control the cross-validation type and 
% other aspects of the optimization, use the <https://www.mathworks.com/help/stats/fitrgp.html#butnn96_sep_shared-HyperparameterOptimizationOptions 
% |HyperparameterOptimizationOptions|> name-value argument. When you use |HyperparameterOptimizationOptions|, 
% you can use the (compact) model size instead of the cross-validation loss as 
% the optimization objective by setting the |ConstraintType| and |ConstraintBounds| 
% options.
% 
% 
% 
% Note: The values of |OptimizeHyperparameters| override any values you specify 
% using other name-value arguments. For example, setting |OptimizeHyperparameters| 
% to |"auto"| causes |fitrgp| to optimize hyperparameters corresponding to the 
% |"auto"| option and to ignore any specified values for the hyperparameters.
% 
% 
% 
% 
% 
% COME BACK TO : Aggregate optimization results for multiple optimization problems, 
% returned as an <https://www.mathworks.com/help/stats/classreg.learning.paramoptim.aggregatebayesianoptimization.html 
% |AggregateBayesianOptimization|> object. To return |AggregateOptimizationResults|, 
% you must specify <https://www.mathworks.com/help/stats/fitrgp.html#butnn96-OptimizeHyperparameters 
% |OptimizeHyperparameters|> and <https://www.mathworks.com/help/stats/fitrgp.html#butnn96_sep_shared-HyperparameterOptimizationOptions 
% |HyperparameterOptimizationOptions|>. You must also specify the |ConstraintType| 
% and |ConstraintBounds| options of |HyperparameterOptimizationOptions|. For an 
% example that shows how to produce this output, see <https://www.mathworks.com/help/stats/classreg.learning.paramoptim.aggregatebayesianoptimization.html#mw_18fb71bf-a209-43e6-8055-f198157c9fe3 
% Hyperparameter Optimization with Multiple Constraint Bounds>.


% Compute marginal log liklihood
% 
% For later: cross validation
% For later: optimizing with minfunc
% 
% 
% 
%