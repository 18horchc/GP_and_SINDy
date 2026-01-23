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

% Standardize so I can use mean 0 (try with and try without)
% Calculate mean and standard deviation of log-transformed data
mu_train = mean(datapointsM1_log);
sigma_train = std(datapointsM1_log);

% Z-score standardization: z = (y_log - mu_train) / sigma_train
datapointsM1_standardized = (datapointsM1_log - mu_train) / sigma_train;

fprintf('Original M1 data: min = %.2f, max = %.2f, mean = %.2f, std = %.2f\n', ...
    min(datapointsM1), max(datapointsM1), mean(datapointsM1), std(datapointsM1));
fprintf('Log-transformed M1 data: min = %.4f, max = %.4f, mean = %.4f, std = %.4f\n', ...
    min(datapointsM1_log), max(datapointsM1_log), mean(datapointsM1_log), std(datapointsM1_log));
fprintf('Standardized M1 data: min = %.4f, max = %.4f, mean = %.6f, std = %.6f\n', ...
    min(datapointsM1_standardized), max(datapointsM1_standardized), ...
    mean(datapointsM1_standardized), std(datapointsM1_standardized));
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
gprMdl_SE_M1 = fitrgp(tpoints_M1, datapointsM1_standardized, ...
    'KernelFunction', 'squaredexponential', 'Standardize', false);
[ypred_SE_M1, std_SE_M1] = predict(gprMdl_SE_M1,xp);


%Exponential (aka Matern 1/2)
gprMdl_exp_M1 = fitrgp(tpoints_M1, datapointsM1_standardized, ...
    'KernelFunction', 'exponential', 'Standardize', false);
[ypred_exp_M1, std_exp_M1] = predict(gprMdl_exp_M1,xp);


%Matern 3/2 [v=3/2] (hyperparams: v, l)
gprMdl_M32_M1 = fitrgp(tpoints_M1, datapointsM1_standardized, ...
    'KernelFunction', 'matern32', 'Standardize', false);
[ypred_M32_M1, std_M32_M1] = predict(gprMdl_M32_M1,xp);

%Matern 5/2 [v=5/2] (hyperparams: v, l)
gprMdl_M52_M1 = fitrgp(tpoints_M1, datapointsM1_standardized, ...
    'KernelFunction', 'matern52', 'Standardize', false);
[ypred_M52_M1, std_M52_M1] = predict(gprMdl_M52_M1,xp);

%Rational Quadratic (hyperparams: alpha, l)
gprMdl_RQ_M1 = fitrgp(tpoints_M1, datapointsM1_standardized, ...
    'KernelFunction', 'rationalquadratic', 'Standardize', false);
[ypred_RQ_M1, std_RQ_M1] = predict(gprMdl_RQ_M1,xp);

%ARD = automatic relevance determination
%ARD squared exponential
gprMdl_ardSE_M1 = fitrgp(tpoints_M1, datapointsM1_standardized, ...
    'KernelFunction', 'ardsquaredexponential', 'Standardize', false);
[ypred_ardSE_M1, std_ardSE_M1] = predict(gprMdl_ardSE_M1,xp);

%ARD exponential
gprMdl_ardExp_M1 = fitrgp(tpoints_M1, datapointsM1_standardized, ...
    'KernelFunction', 'ardexponential', 'Standardize', false);
[ypred_ardExp_M1, std_ardExp_M1] = predict(gprMdl_ardExp_M1,xp);

%ARD Matern 3/2
gprMdl_ardM32_M1 = fitrgp(tpoints_M1, datapointsM1_standardized, ...
    'KernelFunction', 'ardmatern32', 'Standardize', false);
[ypred_ardM32_M1, std_ardM32_M1] = predict(gprMdl_ardM32_M1,xp);


%ARD Matern 5/2
gprMdl_ardM52_M1 = fitrgp(tpoints_M1, datapointsM1_standardized, ...
    'KernelFunction', 'ardmatern52', 'Standardize', false);
[ypred_ardM52_M1, std_ardM52_M1] = predict(gprMdl_ardM52_M1,xp);

%ARD rational quadratic
gprMdl_ardRQ_M1 = fitrgp(tpoints_M1, datapointsM1_standardized, ...
    'KernelFunction', 'ardrationalquadratic', 'Standardize', false);
[ypred_ardRQ_M1, std_ardRQ_M1] = predict(gprMdl_ardRQ_M1,xp);

%% Reverse Transformations: Convert predictions back to original count space
% Step 1: Reverse Z-Score (Un-standardize) to get back to log-space
% μ_log = (μ_z × σ_train) + μ_train
% σ²_log = σ²_z × σ²_train (note: std_SE_M1 is σ_z, so σ²_z = std_SE_M1²)

% Helper function to reverse transform predictions
% Input: ypred_z (mean in standardized space), std_z (std in standardized space)
% Output: ypred_count (mean count), ypred_count_lower (lower bound), ypred_count_upper (upper bound)
reverseTransform = @(ypred_z, std_z) struct(...
    'mu_log', ypred_z * sigma_train + mu_train, ...
    'sigma_log', std_z * sigma_train, ...
    'mu_count_mean', exp(ypred_z * sigma_train + mu_train + 0.5 * (std_z * sigma_train).^2) - 1, ...
    'mu_count_median', exp(ypred_z * sigma_train + mu_train) - 1, ...
    'lower_log', (ypred_z - 1.96 * std_z) * sigma_train + mu_train, ...
    'upper_log', (ypred_z + 1.96 * std_z) * sigma_train + mu_train, ...
    'lower_count', exp((ypred_z - 1.96 * std_z) * sigma_train + mu_train) - 1, ...
    'upper_count', exp((ypred_z + 1.96 * std_z) * sigma_train + mu_train) - 1);

% Transform all predictions back to count space
pred_SE = reverseTransform(ypred_SE_M1, std_SE_M1);
pred_exp = reverseTransform(ypred_exp_M1, std_exp_M1);
pred_M32 = reverseTransform(ypred_M32_M1, std_M32_M1);
pred_M52 = reverseTransform(ypred_M52_M1, std_M52_M1);
pred_RQ = reverseTransform(ypred_RQ_M1, std_RQ_M1);
pred_ardSE = reverseTransform(ypred_ardSE_M1, std_ardSE_M1);
pred_ardExp = reverseTransform(ypred_ardExp_M1, std_ardExp_M1);
pred_ardM32 = reverseTransform(ypred_ardM32_M1, std_ardM32_M1);
pred_ardM52 = reverseTransform(ypred_ardM52_M1, std_ardM52_M1);
pred_ardRQ = reverseTransform(ypred_ardRQ_M1, std_ardRQ_M1);

%%To explore later %%
%non-stationary kernels
%combinging kernels through sums or products
% 
%%  Plot all GP models

% Squared Exponential
% Plot original count data and reverse-transformed predictions
figure;
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_SE.mu_count_median, 'g', 'LineWidth', 1.5);  % Median (most likely path)
plot(xp, pred_SE.lower_count, 'g--', 'LineWidth', 1);     % Lower 95% CI
plot(xp, pred_SE.upper_count, 'g--', 'LineWidth', 1);      % Upper 95% CI
xlabel('time');
ylabel('M1 cell count');
title('GPR with Squared Exponential Kernel');
legend('Data', 'GPR predictions (SE kernel, median)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;

% Exponential
figure;
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_exp.mu_count_median, 'g', 'LineWidth', 1.5);
plot(xp, pred_exp.lower_count, 'g--', 'LineWidth', 1);
plot(xp, pred_exp.upper_count, 'g--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with Exponential Kernel');
legend('Data', 'GPR predictions (Exponential kernel, median)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;

% Matern 3/2
figure;
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_M32.mu_count_median, 'g', 'LineWidth', 1.5);
plot(xp, pred_M32.lower_count, 'g--', 'LineWidth', 1);
plot(xp, pred_M32.upper_count, 'g--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with Matern 3/2 Kernel');
legend('Data', 'GPR predictions (Matern 3/2 kernel, median)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;

% Matern 5/2
figure;
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_M52.mu_count_median, 'g', 'LineWidth', 1.5);
plot(xp, pred_M52.lower_count, 'g--', 'LineWidth', 1);
plot(xp, pred_M52.upper_count, 'g--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with Matern 5/2 Kernel');
legend('Data', 'GPR predictions (Matern 5/2 kernel, median)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;

% Rational Quadratic
figure;
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_RQ.mu_count_median, 'g', 'LineWidth', 1.5);
plot(xp, pred_RQ.lower_count, 'g--', 'LineWidth', 1);
plot(xp, pred_RQ.upper_count, 'g--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with Rational Quadratic Kernel');
legend('Data', 'GPR predictions (Rational Quadratic kernel, median)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;

% ARD Squared Exponential
figure;
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_ardSE.mu_count_median, 'g', 'LineWidth', 1.5);
plot(xp, pred_ardSE.lower_count, 'g--', 'LineWidth', 1);
plot(xp, pred_ardSE.upper_count, 'g--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with ARD Squared Exponential Kernel');
legend('Data', 'GPR predictions (ARD SE kernel, median)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;

% ARD Exponential
figure;
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_ardExp.mu_count_median, 'g', 'LineWidth', 1.5);
plot(xp, pred_ardExp.lower_count, 'g--', 'LineWidth', 1);
plot(xp, pred_ardExp.upper_count, 'g--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with ARD Exponential Kernel');
legend('Data', 'GPR predictions (ARD Exponential kernel, median)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;

% ARD Matern 3/2
figure;
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_ardM32.mu_count_median, 'g', 'LineWidth', 1.5);
plot(xp, pred_ardM32.lower_count, 'g--', 'LineWidth', 1);
plot(xp, pred_ardM32.upper_count, 'g--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with ARD Matern 3/2 Kernel');
legend('Data', 'GPR predictions (ARD Matern 3/2 kernel, median)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;

% ARD Matern 5/2
figure;
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_ardM52.mu_count_median, 'g', 'LineWidth', 1.5);
plot(xp, pred_ardM52.lower_count, 'g--', 'LineWidth', 1);
plot(xp, pred_ardM52.upper_count, 'g--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with ARD Matern 5/2 Kernel');
legend('Data', 'GPR predictions (ARD Matern 5/2 kernel, median)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;

% ARD Rational Quadratic
figure;
plot(tpoints_M1, datapointsM1, 'b.', 'MarkerSize', 10);
hold on;
plot(xp, pred_ardRQ.mu_count_median, 'g', 'LineWidth', 1.5);
plot(xp, pred_ardRQ.lower_count, 'g--', 'LineWidth', 1);
plot(xp, pred_ardRQ.upper_count, 'g--', 'LineWidth', 1);
xlabel('time');
ylabel('M1 cell count');
title('GPR with ARD Rational Quadratic Kernel');
legend('Data', 'GPR predictions (ARD Rational Quadratic kernel, median)', '95% Confidence interval', 'Location', 'best');
hold off;
grid on;
% 
%% Compute marginal log liklihood

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

% Store all models and their names
models = {
    gprMdl_SE_M1, 'Squared Exponential';
    gprMdl_exp_M1, 'Exponential';
    gprMdl_M32_M1, 'Matern 3/2';
    gprMdl_M52_M1, 'Matern 5/2';
    gprMdl_RQ_M1, 'Rational Quadratic';
    gprMdl_ardSE_M1, 'ARD Squared Exponential';
    gprMdl_ardExp_M1, 'ARD Exponential';
    gprMdl_ardM32_M1, 'ARD Matern 3/2';
    gprMdl_ardM52_M1, 'ARD Matern 5/2';
    gprMdl_ardRQ_M1, 'ARD Rational Quadratic'
};

% Initialize arrays to store metrics
numModels = size(models, 1);
kernelNames = cell(numModels, 1);
logLikelihoods = zeros(numModels, 1);
RMSEs = zeros(numModels, 1);
R2s = zeros(numModels, 1);

% Calculate metrics for each model
for i = 1:numModels
    model = models{i, 1};
    kernelNames{i} = models{i, 2};
    
    % Get log likelihood (property of the model)
    logLikelihoods(i) = model.LogLikelihood;
    
    % Get predictions on training data
    yPredict = resubPredict(model);
    yActual = datapointsM1_standardized;  % Actual training data (log-transformed and standardized)
    
    % Calculate residuals
    residuals = yActual - yPredict;
    
    % Calculate RMSE
    RMSEs(i) = sqrt(mean(residuals.^2));
    
    % Calculate R²
    SSR = sum((yPredict - mean(yActual)).^2);  % Regression sum of squares
    SST = sum((yActual - mean(yActual)).^2);   % Total sum of squares
    R2s(i) = SSR / SST;
end

% Create table with results
resultsTable = table(kernelNames, logLikelihoods, RMSEs, R2s, ...
    'VariableNames', {'Kernel', 'LogLikelihood', 'RMSE', 'R2'});

% Display the table
disp('GPR Model Comparison Results for M1 Data:');
disp(resultsTable);

% Optionally sort by a metric (e.g., by LogLikelihood descending)
resultsTableSorted = sortrows(resultsTable, 'LogLikelihood', 'descend');
disp('Sorted by LogLikelihood (best to worst):');
disp(resultsTableSorted);
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