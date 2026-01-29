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

%TO DO: 
% try defining / constraining initial hyperparameters



%% Data

clear; clc; close all; 
%samples
tpoints_M1 = [0,0,1,1,2,2,3,3,3,3,3,5,7,7,7,7,14,14,14]'; %time points where experimental data is observed
datapointsM1 = [0,10,5,50,120,125,375,62,102,100,60,325,600,55,900,225,750,400,1300]'; %M1 experimental data

tpoints_M2= [0,0, 1, 1,1, 2, 2, 3,3,3,3,3, 5, 7, 7, 7, 7, 14,14,14]'; %time points where experimental data is observed
datapointsM2 = [0,10,170,15,50, 90,269,300,15,57,100,160, 800, 600,6,400,270,200,110,100]'; %M2 experimental data


%for testing
xp = linspace(0, 14, 300)'; 

%% Optimizer on/off
% Set to true to enable hyperparameter optimization for all GP fits (log and sqrt).
% Set to false to fit with default hyperparameters only.
optimizeHyperparams = false;

if optimizeHyperparams
    optimArgs = {'OptimizeHyperparameters', 'auto', 'HyperparameterOptimizationOptions', ...
        struct('ShowPlots', false, 'Verbose', 0)};
else
    optimArgs = {};
end

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

%% Pre-processing M2 (same steps as M1)
datapointsM2_log = log(datapointsM2 + 1);
datapointsM2_sqrt = sqrt(datapointsM2);
mu_train_log_M2 = mean(datapointsM2_log);
sigma_train_log_M2 = std(datapointsM2_log);
datapointsM2_standardized_log = (datapointsM2_log - mu_train_log_M2) / sigma_train_log_M2;
mu_train_sqrt_M2 = mean(datapointsM2_sqrt);
sigma_train_sqrt_M2 = std(datapointsM2_sqrt);
datapointsM2_standardized_sqrt = (datapointsM2_sqrt - mu_train_sqrt_M2) / sigma_train_sqrt_M2;
fprintf('M2: Original min = %.2f, max = %.2f | Log std = %.4f | Sqrt std = %.4f\n\n', ...
    min(datapointsM2), max(datapointsM2), std(datapointsM2_standardized_log), std(datapointsM2_standardized_sqrt));




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

% Fit GP on standardized log-transformed data (z)


%Squared exponental (hyperparams: lengthsale (l))
gprMdl_SE_M1_log = fitrgp(tpoints_M1, datapointsM1_standardized_log, ...
    'KernelFunction', 'squaredexponential', 'Standardize', false, optimArgs{:});
[ypred_SE_M1_log, std_SE_M1_log] = predict(gprMdl_SE_M1_log,xp);


%Exponential (aka Matern 1/2)
gprMdl_exp_M1_log = fitrgp(tpoints_M1, datapointsM1_standardized_log, ...
    'KernelFunction', 'exponential', 'Standardize', false, optimArgs{:});
[ypred_exp_M1_log, std_exp_M1_log] = predict(gprMdl_exp_M1_log,xp);


%Matern 3/2 [v=3/2] (hyperparams: v, l)
gprMdl_M32_M1_log = fitrgp(tpoints_M1, datapointsM1_standardized_log, ...
    'KernelFunction', 'matern32', 'Standardize', false, optimArgs{:});
[ypred_M32_M1_log, std_M32_M1_log] = predict(gprMdl_M32_M1_log,xp);

%Matern 5/2 [v=5/2] (hyperparams: v, l)
gprMdl_M52_M1_log = fitrgp(tpoints_M1, datapointsM1_standardized_log, ...
    'KernelFunction', 'matern52', 'Standardize', false, optimArgs{:});
[ypred_M52_M1_log, std_M52_M1_log] = predict(gprMdl_M52_M1_log,xp);

%Rational Quadratic (hyperparams: alpha, l)
gprMdl_RQ_M1_log = fitrgp(tpoints_M1, datapointsM1_standardized_log, ...
    'KernelFunction', 'rationalquadratic', 'Standardize', false, optimArgs{:});
[ypred_RQ_M1_log, std_RQ_M1_log] = predict(gprMdl_RQ_M1_log,xp);

%% GP Models with Square Root Transform
% Fit GP models on standardized sqrt-transformed data
% Fit GP on standardized sqrt-transformed data (z)

%Squared exponental (hyperparams: lengthsale (l))
gprMdl_SE_M1_sqrt = fitrgp(tpoints_M1, datapointsM1_standardized_sqrt, ...
    'KernelFunction', 'squaredexponential', 'Standardize', false, optimArgs{:});
[ypred_SE_M1_sqrt, std_SE_M1_sqrt] = predict(gprMdl_SE_M1_sqrt,xp);

%Exponential (aka Matern 1/2)
gprMdl_exp_M1_sqrt = fitrgp(tpoints_M1, datapointsM1_standardized_sqrt, ...
    'KernelFunction', 'exponential', 'Standardize', false, optimArgs{:});
[ypred_exp_M1_sqrt, std_exp_M1_sqrt] = predict(gprMdl_exp_M1_sqrt,xp);

%Matern 3/2 [v=3/2] (hyperparams: v, l)
gprMdl_M32_M1_sqrt = fitrgp(tpoints_M1, datapointsM1_standardized_sqrt, ...
    'KernelFunction', 'matern32', 'Standardize', false, optimArgs{:});
[ypred_M32_M1_sqrt, std_M32_M1_sqrt] = predict(gprMdl_M32_M1_sqrt,xp);

%Matern 5/2 [v=5/2] (hyperparams: v, l)
gprMdl_M52_M1_sqrt = fitrgp(tpoints_M1, datapointsM1_standardized_sqrt, ...
    'KernelFunction', 'matern52', 'Standardize', false, optimArgs{:});
[ypred_M52_M1_sqrt, std_M52_M1_sqrt] = predict(gprMdl_M52_M1_sqrt,xp);

%Rational Quadratic (hyperparams: alpha, l)
gprMdl_RQ_M1_sqrt = fitrgp(tpoints_M1, datapointsM1_standardized_sqrt, ...
    'KernelFunction', 'rationalquadratic', 'Standardize', false, optimArgs{:});
[ypred_RQ_M1_sqrt, std_RQ_M1_sqrt] = predict(gprMdl_RQ_M1_sqrt,xp);

%% GP Models M2 (Log and Sqrt transforms)
% Log transform models
gprMdl_SE_M2_log = fitrgp(tpoints_M2, datapointsM2_standardized_log, ...
    'KernelFunction', 'squaredexponential', 'Standardize', false, optimArgs{:});
[ypred_SE_M2_log, std_SE_M2_log] = predict(gprMdl_SE_M2_log,xp);
gprMdl_exp_M2_log = fitrgp(tpoints_M2, datapointsM2_standardized_log, ...
    'KernelFunction', 'exponential', 'Standardize', false, optimArgs{:});
[ypred_exp_M2_log, std_exp_M2_log] = predict(gprMdl_exp_M2_log,xp);
gprMdl_M32_M2_log = fitrgp(tpoints_M2, datapointsM2_standardized_log, ...
    'KernelFunction', 'matern32', 'Standardize', false, optimArgs{:});
[ypred_M32_M2_log, std_M32_M2_log] = predict(gprMdl_M32_M2_log,xp);
gprMdl_M52_M2_log = fitrgp(tpoints_M2, datapointsM2_standardized_log, ...
    'KernelFunction', 'matern52', 'Standardize', false, optimArgs{:});
[ypred_M52_M2_log, std_M52_M2_log] = predict(gprMdl_M52_M2_log,xp);
gprMdl_RQ_M2_log = fitrgp(tpoints_M2, datapointsM2_standardized_log, ...
    'KernelFunction', 'rationalquadratic', 'Standardize', false, optimArgs{:});
[ypred_RQ_M2_log, std_RQ_M2_log] = predict(gprMdl_RQ_M2_log,xp);
% Sqrt transform models
gprMdl_SE_M2_sqrt = fitrgp(tpoints_M2, datapointsM2_standardized_sqrt, ...
    'KernelFunction', 'squaredexponential', 'Standardize', false, optimArgs{:});
[ypred_SE_M2_sqrt, std_SE_M2_sqrt] = predict(gprMdl_SE_M2_sqrt,xp);
gprMdl_exp_M2_sqrt = fitrgp(tpoints_M2, datapointsM2_standardized_sqrt, ...
    'KernelFunction', 'exponential', 'Standardize', false, optimArgs{:});
[ypred_exp_M2_sqrt, std_exp_M2_sqrt] = predict(gprMdl_exp_M2_sqrt,xp);
gprMdl_M32_M2_sqrt = fitrgp(tpoints_M2, datapointsM2_standardized_sqrt, ...
    'KernelFunction', 'matern32', 'Standardize', false, optimArgs{:});
[ypred_M32_M2_sqrt, std_M32_M2_sqrt] = predict(gprMdl_M32_M2_sqrt,xp);
gprMdl_M52_M2_sqrt = fitrgp(tpoints_M2, datapointsM2_standardized_sqrt, ...
    'KernelFunction', 'matern52', 'Standardize', false, optimArgs{:});
[ypred_M52_M2_sqrt, std_M52_M2_sqrt] = predict(gprMdl_M52_M2_sqrt,xp);
gprMdl_RQ_M2_sqrt = fitrgp(tpoints_M2, datapointsM2_standardized_sqrt, ...
    'KernelFunction', 'rationalquadratic', 'Standardize', false, optimArgs{:});
[ypred_RQ_M2_sqrt, std_RQ_M2_sqrt] = predict(gprMdl_RQ_M2_sqrt,xp);

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

% M2 reverse transform helpers (use M2 train stats)
reverseTransformLog_M2 = @(ypred_z, std_z) struct(...
    'mu_log', ypred_z * sigma_train_log_M2 + mu_train_log_M2, ...
    'sigma_log', std_z * sigma_train_log_M2, ...
    'mu_count_mean', exp(ypred_z * sigma_train_log_M2 + mu_train_log_M2 + 0.5 * (std_z * sigma_train_log_M2).^2) - 1, ...
    'mu_count_median', exp(ypred_z * sigma_train_log_M2 + mu_train_log_M2) - 1, ...
    'lower_log', (ypred_z - 1.96 * std_z) * sigma_train_log_M2 + mu_train_log_M2, ...
    'upper_log', (ypred_z + 1.96 * std_z) * sigma_train_log_M2 + mu_train_log_M2, ...
    'lower_count', exp((ypred_z - 1.96 * std_z) * sigma_train_log_M2 + mu_train_log_M2) - 1, ...
    'upper_count', exp((ypred_z + 1.96 * std_z) * sigma_train_log_M2 + mu_train_log_M2) - 1);
reverseTransformSqrt_M2 = @(ypred_z, std_z) struct(...
    'mu_sqrt', ypred_z * sigma_train_sqrt_M2 + mu_train_sqrt_M2, ...
    'sigma_sqrt', std_z * sigma_train_sqrt_M2, ...
    'mu_count_mean', (ypred_z * sigma_train_sqrt_M2 + mu_train_sqrt_M2).^2 + (std_z * sigma_train_sqrt_M2).^2, ...
    'mu_count_median', (ypred_z * sigma_train_sqrt_M2 + mu_train_sqrt_M2).^2, ...
    'lower_sqrt', (ypred_z - 1.96 * std_z) * sigma_train_sqrt_M2 + mu_train_sqrt_M2, ...
    'upper_sqrt', (ypred_z + 1.96 * std_z) * sigma_train_sqrt_M2 + mu_train_sqrt_M2, ...
    'lower_count', max(0, (ypred_z - 1.96 * std_z) * sigma_train_sqrt_M2 + mu_train_sqrt_M2).^2, ...
    'upper_count', ((ypred_z + 1.96 * std_z) * sigma_train_sqrt_M2 + mu_train_sqrt_M2).^2);

% Transform all predictions back to count space (log transform)
pred_SE_log = reverseTransformLog(ypred_SE_M1_log, std_SE_M1_log);
pred_exp_log = reverseTransformLog(ypred_exp_M1_log, std_exp_M1_log);
pred_M32_log = reverseTransformLog(ypred_M32_M1_log, std_M32_M1_log);
pred_M52_log = reverseTransformLog(ypred_M52_M1_log, std_M52_M1_log);
pred_RQ_log = reverseTransformLog(ypred_RQ_M1_log, std_RQ_M1_log);

% Transform all predictions back to count space (sqrt transform)
pred_SE_sqrt = reverseTransformSqrt(ypred_SE_M1_sqrt, std_SE_M1_sqrt);
pred_exp_sqrt = reverseTransformSqrt(ypred_exp_M1_sqrt, std_exp_M1_sqrt);
pred_M32_sqrt = reverseTransformSqrt(ypred_M32_M1_sqrt, std_M32_M1_sqrt);
pred_M52_sqrt = reverseTransformSqrt(ypred_M52_M1_sqrt, std_M52_M1_sqrt);
pred_RQ_sqrt = reverseTransformSqrt(ypred_RQ_M1_sqrt, std_RQ_M1_sqrt);

% Transform M2 predictions back to count space
pred_SE_log_M2 = reverseTransformLog_M2(ypred_SE_M2_log, std_SE_M2_log);
pred_exp_log_M2 = reverseTransformLog_M2(ypred_exp_M2_log, std_exp_M2_log);
pred_M32_log_M2 = reverseTransformLog_M2(ypred_M32_M2_log, std_M32_M2_log);
pred_M52_log_M2 = reverseTransformLog_M2(ypred_M52_M2_log, std_M52_M2_log);
pred_RQ_log_M2 = reverseTransformLog_M2(ypred_RQ_M2_log, std_RQ_M2_log);
pred_SE_sqrt_M2 = reverseTransformSqrt_M2(ypred_SE_M2_sqrt, std_SE_M2_sqrt);
pred_exp_sqrt_M2 = reverseTransformSqrt_M2(ypred_exp_M2_sqrt, std_exp_M2_sqrt);
pred_M32_sqrt_M2 = reverseTransformSqrt_M2(ypred_M32_M2_sqrt, std_M32_M2_sqrt);
pred_M52_sqrt_M2 = reverseTransformSqrt_M2(ypred_M52_M2_sqrt, std_M52_M2_sqrt);
pred_RQ_sqrt_M2 = reverseTransformSqrt_M2(ypred_RQ_M2_sqrt, std_RQ_M2_sqrt);

%%  Plot all GP models (M1: one figure with 5 tabs; M2: separate figure with 5 tabs)

% M1: one figure with 5 tabs (each tab = one kernel, log + sqrt subplots)
pred_log_M1 = {pred_SE_log, pred_exp_log, pred_M32_log, pred_M52_log, pred_RQ_log};
pred_sqrt_M1 = {pred_SE_sqrt, pred_exp_sqrt, pred_M32_sqrt, pred_M52_sqrt, pred_RQ_sqrt};
kernelTabNames = {'Squared Exp', 'Exponential', 'Matern 3/2', 'Matern 5/2', 'Rational Quad'};

figM1 = figure('Name', 'M1 GP Kernels', 'NumberTitle', 'off');
tgM1 = uitabgroup(figM1);
for k = 1:5
    t = uitab(tgM1, 'Title', kernelTabNames{k});
    ax1 = axes('Parent', t, 'Position', [0.1 0.55 0.85 0.38]);
    plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
    hold(ax1, 'on');
    plot(ax1, xp, pred_log_M1{k}.mu_count_median, 'g', 'LineWidth', 1.5);
    plot(ax1, xp, pred_log_M1{k}.lower_count, 'g--', 'LineWidth', 1);
    plot(ax1, xp, pred_log_M1{k}.upper_count, 'g--', 'LineWidth', 1);
    ylabel(ax1, 'M1 cell count');
    title(ax1, 'Log transform (median)');
    legend(ax1, 'Data', 'GPR (median)', '95% CI', 'Location', 'best');
    hold(ax1, 'off');
    grid(ax1, 'on');
    ax2 = axes('Parent', t, 'Position', [0.1 0.08 0.85 0.38]);
    plot(ax2, tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
    hold(ax2, 'on');
    plot(ax2, xp, pred_sqrt_M1{k}.mu_count_mean, 'r', 'LineWidth', 1.5);
    plot(ax2, xp, pred_sqrt_M1{k}.lower_count, 'r--', 'LineWidth', 1);
    plot(ax2, xp, pred_sqrt_M1{k}.upper_count, 'r--', 'LineWidth', 1);
    xlabel(ax2, 'time');
    ylabel(ax2, 'M1 cell count');
    title(ax2, 'Sqrt transform (mean)');
    legend(ax2, 'Data', 'GPR (mean)', '95% CI', 'Location', 'best');
    hold(ax2, 'off');
    grid(ax2, 'on');
end

% M2: separate figure with 5 tabs (same structure)
pred_log_M2 = {pred_SE_log_M2, pred_exp_log_M2, pred_M32_log_M2, pred_M52_log_M2, pred_RQ_log_M2};
pred_sqrt_M2 = {pred_SE_sqrt_M2, pred_exp_sqrt_M2, pred_M32_sqrt_M2, pred_M52_sqrt_M2, pred_RQ_sqrt_M2};

figM2 = figure('Name', 'M2 GP Kernels', 'NumberTitle', 'off');
tgM2 = uitabgroup(figM2);
for k = 1:5
    t = uitab(tgM2, 'Title', kernelTabNames{k});
    ax1 = axes('Parent', t, 'Position', [0.1 0.55 0.85 0.38]);
    plot(ax1, tpoints_M2, datapointsM2, 'b.', 'MarkerSize', 10);
    hold(ax1, 'on');
    plot(ax1, xp, pred_log_M2{k}.mu_count_median, 'g', 'LineWidth', 1.5);
    plot(ax1, xp, pred_log_M2{k}.lower_count, 'g--', 'LineWidth', 1);
    plot(ax1, xp, pred_log_M2{k}.upper_count, 'g--', 'LineWidth', 1);
    ylabel(ax1, 'M2 cell count');
    title(ax1, 'Log transform (median)');
    legend(ax1, 'Data', 'GPR (median)', '95% CI', 'Location', 'best');
    hold(ax1, 'off');
    grid(ax1, 'on');
    ax2 = axes('Parent', t, 'Position', [0.1 0.08 0.85 0.38]);
    plot(ax2, tpoints_M2, datapointsM2, 'b.', 'MarkerSize', 10);
    hold(ax2, 'on');
    plot(ax2, xp, pred_sqrt_M2{k}.mu_count_mean, 'r', 'LineWidth', 1.5);
    plot(ax2, xp, pred_sqrt_M2{k}.lower_count, 'r--', 'LineWidth', 1);
    plot(ax2, xp, pred_sqrt_M2{k}.upper_count, 'r--', 'LineWidth', 1);
    xlabel(ax2, 'time');
    ylabel(ax2, 'M2 cell count');
    title(ax2, 'Sqrt transform (mean)');
    legend(ax2, 'Data', 'GPR (mean)', '95% CI', 'Location', 'best');
    hold(ax2, 'off');
    grid(ax2, 'on');
end

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
    gprMdl_RQ_M1_log, 'Rational Quadratic'
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
    gprMdl_RQ_M1_sqrt, 'Rational Quadratic'
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

%% Metrics Table for M2 (Log and Sqrt)
% M2 Log transform
models_log_M2 = {
    gprMdl_SE_M2_log, 'Squared Exponential';
    gprMdl_exp_M2_log, 'Exponential';
    gprMdl_M32_M2_log, 'Matern 3/2';
    gprMdl_M52_M2_log, 'Matern 5/2';
    gprMdl_RQ_M2_log, 'Rational Quadratic'
};
kernelNames_log_M2 = cell(numModels, 1);
logLikelihoods_log_M2 = zeros(numModels, 1);
RMSEs_log_M2 = zeros(numModels, 1);
R2s_log_M2 = zeros(numModels, 1);
for i = 1:numModels
    model = models_log_M2{i, 1};
    kernelNames_log_M2{i} = models_log_M2{i, 2};
    logLikelihoods_log_M2(i) = model.LogLikelihood;
    yPredict = resubPredict(model);
    yActual = datapointsM2_standardized_log;
    residuals = yActual - yPredict;
    RMSEs_log_M2(i) = sqrt(mean(residuals.^2));
    SSR = sum((yPredict - mean(yActual)).^2);
    SST = sum((yActual - mean(yActual)).^2);
    R2s_log_M2(i) = SSR / SST;
end
resultsTable_log_M2 = table(kernelNames_log_M2, logLikelihoods_log_M2, RMSEs_log_M2, R2s_log_M2, ...
    'VariableNames', {'Kernel', 'LogLikelihood', 'RMSE', 'R2'});
disp('GPR Model Comparison Results for M2 Data (Log Transform):');
disp(resultsTable_log_M2);
disp('Sorted by LogLikelihood (best to worst):');
disp(sortrows(resultsTable_log_M2, 'LogLikelihood', 'descend'));

% M2 Sqrt transform
models_sqrt_M2 = {
    gprMdl_SE_M2_sqrt, 'Squared Exponential';
    gprMdl_exp_M2_sqrt, 'Exponential';
    gprMdl_M32_M2_sqrt, 'Matern 3/2';
    gprMdl_M52_M2_sqrt, 'Matern 5/2';
    gprMdl_RQ_M2_sqrt, 'Rational Quadratic'
};
kernelNames_sqrt_M2 = cell(numModels, 1);
logLikelihoods_sqrt_M2 = zeros(numModels, 1);
RMSEs_sqrt_M2 = zeros(numModels, 1);
R2s_sqrt_M2 = zeros(numModels, 1);
for i = 1:numModels
    model = models_sqrt_M2{i, 1};
    kernelNames_sqrt_M2{i} = models_sqrt_M2{i, 2};
    logLikelihoods_sqrt_M2(i) = model.LogLikelihood;
    yPredict = resubPredict(model);
    yActual = datapointsM2_standardized_sqrt;
    residuals = yActual - yPredict;
    RMSEs_sqrt_M2(i) = sqrt(mean(residuals.^2));
    SSR = sum((yPredict - mean(yActual)).^2);
    SST = sum((yActual - mean(yActual)).^2);
    R2s_sqrt_M2(i) = SSR / SST;
end
resultsTable_sqrt_M2 = table(kernelNames_sqrt_M2, logLikelihoods_sqrt_M2, RMSEs_sqrt_M2, R2s_sqrt_M2, ...
    'VariableNames', {'Kernel', 'LogLikelihood', 'RMSE', 'R2'});
disp('GPR Model Comparison Results for M2 Data (Sqrt Transform):');
disp(resultsTable_sqrt_M2);
disp('Sorted by LogLikelihood (best to worst):');
disp(sortrows(resultsTable_sqrt_M2, 'LogLikelihood', 'descend'));

%% Hyperparameter tables (KernelParameters + Sigma) for each model
% Build hyperparameter string from model.KernelInformation.KernelParameterNames,
% model.KernelInformation.KernelParameters, and model.Sigma
kernelNames_hp = cell(numModels, 1);
hyperparamStrs = cell(numModels, 1);

% M1 Log
for i = 1:numModels
    model = models_log{i, 1};
    kernelNames_hp{i} = models_log{i, 2};
    pNames = cellstr(model.KernelInformation.KernelParameterNames);
    pVals = model.KernelInformation.KernelParameters;
    parts = arrayfun(@(j) sprintf('%s=%g', pNames{j}, pVals(j)), 1:numel(pNames), 'UniformOutput', false);
    hyperparamStrs{i} = [strjoin(parts, ', '), sprintf(', Sigma=%g', model.Sigma)];
end
hyperparamTable_M1_log = table(kernelNames_hp, hyperparamStrs, 'VariableNames', {'Kernel', 'Hyperparameters'});
disp('Hyperparameters used (M1, Log transform):');
disp(hyperparamTable_M1_log);

% M1 Sqrt
for i = 1:numModels
    model = models_sqrt{i, 1};
    kernelNames_hp{i} = models_sqrt{i, 2};
    pNames = cellstr(model.KernelInformation.KernelParameterNames);
    pVals = model.KernelInformation.KernelParameters;
    parts = arrayfun(@(j) sprintf('%s=%g', pNames{j}, pVals(j)), 1:numel(pNames), 'UniformOutput', false);
    hyperparamStrs{i} = [strjoin(parts, ', '), sprintf(', Sigma=%g', model.Sigma)];
end
hyperparamTable_M1_sqrt = table(kernelNames_hp, hyperparamStrs, 'VariableNames', {'Kernel', 'Hyperparameters'});
disp('Hyperparameters used (M1, Sqrt transform):');
disp(hyperparamTable_M1_sqrt);

% M2 Log
for i = 1:numModels
    model = models_log_M2{i, 1};
    kernelNames_hp{i} = models_log_M2{i, 2};
    pNames = cellstr(model.KernelInformation.KernelParameterNames);
    pVals = model.KernelInformation.KernelParameters;
    parts = arrayfun(@(j) sprintf('%s=%g', pNames{j}, pVals(j)), 1:numel(pNames), 'UniformOutput', false);
    hyperparamStrs{i} = [strjoin(parts, ', '), sprintf(', Sigma=%g', model.Sigma)];
end
hyperparamTable_M2_log = table(kernelNames_hp, hyperparamStrs, 'VariableNames', {'Kernel', 'Hyperparameters'});
disp('Hyperparameters used (M2, Log transform):');
disp(hyperparamTable_M2_log);

% M2 Sqrt
for i = 1:numModels
    model = models_sqrt_M2{i, 1};
    kernelNames_hp{i} = models_sqrt_M2{i, 2};
    pNames = cellstr(model.KernelInformation.KernelParameterNames);
    pVals = model.KernelInformation.KernelParameters;
    parts = arrayfun(@(j) sprintf('%s=%g', pNames{j}, pVals(j)), 1:numel(pNames), 'UniformOutput', false);
    hyperparamStrs{i} = [strjoin(parts, ', '), sprintf(', Sigma=%g', model.Sigma)];
end
hyperparamTable_M2_sqrt = table(kernelNames_hp, hyperparamStrs, 'VariableNames', {'Kernel', 'Hyperparameters'});
disp('Hyperparameters used (M2, Sqrt transform):');
disp(hyperparamTable_M2_sqrt);

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