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
% normalize and standardize so I can use mean 1 (try with and try without)

% sqrt transform to ensure positive


% log transform to ensure positive

%% Understanding Basis Functions
% The basis function in GP regression represents the mean function (prior mean)
% before the GP learns from data. It's the "expected value" we assume before
% seeing any observations.
%
% Options in fitrgp:
% 1. 'constant' (default): μ(x) = β₀ (a single constant)
%    - Best when: You expect the data to fluctuate around a constant mean
%    - Use when: Data shows no clear trend, or trend is captured by the kernel
%    - Your microglia data: Good if you expect growth but kernel handles the trend
%
% 2. 'linear': μ(x) = β₀ + β₁x (linear trend)
%    - Best when: You expect a clear linear trend in your data
%    - Use when: Data shows steady increase/decrease over time
%    - Your microglia data: Good if you expect roughly linear growth over time
%
% 3. 'pureQuadratic': μ(x) = β₀ + β₁x + β₂x² (quadratic trend)
%    - Best when: You expect acceleration/deceleration (curved growth)
%    - Use when: Data shows exponential-like or saturating growth
%    - Your microglia data: Good if you expect log growth (which is curved)
%
% 4. 'none': No mean function (μ(x) = 0)
%    - Best when: Kernel alone should capture all structure
%    - Use when: You want the GP to be completely data-driven
%    - Your microglia data: Use if kernel is flexible enough
%
% How to choose for microglia data:
% - If you expect log growth: Try 'pureQuadratic' or 'linear'
% - If data is noisy around a baseline: Use 'constant' (default)
% - If unsure: Test multiple and compare LogLikelihood (higher is better)
% - The basis function affects the prior, but the GP can still learn complex
%   patterns through the kernel covariance function

%% Test Different Basis Functions for Each Kernel
% IMPORTANT: Different kernels can have different optimal basis functions!
% The kernel's smoothness and ability to capture trends interacts with the
% basis function, so we should test basis functions separately for each kernel.
%
% This section tests all basis function options for each kernel and finds
% the best basis function for each kernel type.

fprintf('\n=== Testing Basis Functions for Each Kernel ===\n');

basisOptions = {'constant', 'linear', 'pureQuadratic', 'none'};
kernelNames = {
    'squaredexponential', 'Squared Exponential';
    'exponential', 'Exponential';
    'matern32', 'Matern 3/2';
    'matern52', 'Matern 5/2';
    'rationalquadratic', 'Rational Quadratic';
    'ardsquaredexponential', 'ARD Squared Exponential';
    'ardexponential', 'ARD Exponential';
    'ardmatern32', 'ARD Matern 3/2';
    'ardmatern52', 'ARD Matern 5/2';
    'ardrationalquadratic', 'ARD Rational Quadratic'
};

% Store best basis for each kernel
bestBasisPerKernel = cell(size(kernelNames, 1), 1);
basisTestResults = cell(size(kernelNames, 1), 1);

for k = 1:size(kernelNames, 1)
    kernelFunc = kernelNames{k, 1};
    kernelDisplayName = kernelNames{k, 2};
    
    fprintf('\n--- Testing %s ---\n', kernelDisplayName);
    
    basisResults = cell(length(basisOptions), 4);
    basisResults(:,1) = basisOptions';
    
    for i = 1:length(basisOptions)
        try
            gprTest = fitrgp(tpoints_M1, datapointsM1, ...
                'KernelFunction', kernelFunc, ...
                'BasisFunction', basisOptions{i});
            
            % Get metrics
            yPredict = resubPredict(gprTest);
            residuals = datapointsM1 - yPredict;
            RMSE = sqrt(mean(residuals.^2));
            SSR = sum((yPredict - mean(datapointsM1)).^2);
            SST = sum((datapointsM1 - mean(datapointsM1)).^2);
            R2 = SSR / SST;
            
            basisResults{i,2} = gprTest.LogLikelihood;
            basisResults{i,3} = RMSE;
            basisResults{i,4} = R2;
            
            fprintf('  Basis: %-15s | LogLikelihood: %8.2f | RMSE: %8.2f | R²: %6.4f\n', ...
                basisOptions{i}, gprTest.LogLikelihood, RMSE, R2);
        catch ME
            fprintf('  Basis: %-15s | Error: %s\n', basisOptions{i}, ME.message);
            basisResults{i,2} = NaN;
            basisResults{i,3} = NaN;
            basisResults{i,4} = NaN;
        end
    end
    
    % Find best basis for this kernel (highest LogLikelihood)
    validIdx = ~isnan(cell2mat(basisResults(:,2)));
    if any(validIdx)
        validResults = basisResults(validIdx, :);
        [~, bestIdx] = max(cell2mat(validResults(:,2)));
        bestBasisPerKernel{k} = validResults{bestIdx, 1};
        
        fprintf('  >>> Best basis for %s: %s\n', kernelDisplayName, bestBasisPerKernel{k});
    else
        bestBasisPerKernel{k} = 'constant';  % Default fallback
        fprintf('  >>> Using default basis: constant\n');
    end
    
    % Store results for this kernel
    basisTestResults{k} = basisResults;
end

% Create summary table
kernelBasisSummary = table(kernelNames(:,2), bestBasisPerKernel, ...
    'VariableNames', {'Kernel', 'BestBasis'});
fprintf('\n=== Summary: Best Basis Function for Each Kernel ===\n');
disp(kernelBasisSummary);

%% *Comparing kernels* 
% *Without Bayesian optimizer*
% More info about kernels: <https://www.mathworks.com/help/stats/kernel-covariance-function-options.html 
% https://www.mathworks.com/help/stats/kernel-covariance-function-options.html> 
% 
% Not setting any inital values for kernel parameters
% 
% Basis function: Each kernel uses its optimal basis function determined above.
% This ensures each kernel is evaluated with its best-performing basis function,
% as different kernels can have different optimal basis functions.
%
% Other settings:
% * Using FitMethod = exact (default) _(Exact Gaussian process regression. This 
% value is the default if n ≤ 2000, where n is the number of observations.)_
% * Always using predict = exact _(Exact Gaussian process regression method. 
% This value is the default if n ≤ 10,000.)_
% * Standardize = false (default) - set to true if your data spans very different scales

%Squared exponental (hyperparams: lengthsale (l))
% Using kernel-specific best basis: bestBasisPerKernel{1}
gprMdl_SE_M1 = fitrgp(tpoints_M1, datapointsM1, ...
    'BasisFunction', bestBasisPerKernel{1});
[ypred_SE_M1,~] = predict(gprMdl_SE_M1,xp);


%Exponential (aka Matern 1/2)
% Using kernel-specific best basis: bestBasisPerKernel{2}
gprMdl_exp_M1 = fitrgp(tpoints_M1, datapointsM1, ...
    'KernelFunction', 'exponential', 'BasisFunction', bestBasisPerKernel{2});
[ypred_exp_M1,~] = predict(gprMdl_exp_M1,xp);


%Matern 3/2 [v=3/2] (hyperparams: v, l)
% Using kernel-specific best basis: bestBasisPerKernel{3}
gprMdl_M32_M1 = fitrgp(tpoints_M1, datapointsM1, ...
    'KernelFunction', 'matern32', 'BasisFunction', bestBasisPerKernel{3});
[ypred_M32_M1,~] = predict(gprMdl_M32_M1,xp);

%Matern 5/2 [v=5/2] (hyperparams: v, l)
% Using kernel-specific best basis: bestBasisPerKernel{4}
gprMdl_M52_M1 = fitrgp(tpoints_M1, datapointsM1, ...
    'KernelFunction', 'matern52', 'BasisFunction', bestBasisPerKernel{4});
[ypred_M52_M1,~] = predict(gprMdl_M52_M1,xp);

%Rational Quadratic (hyperparams: alpha, l)
% Using kernel-specific best basis: bestBasisPerKernel{5}
gprMdl_RQ_M1 = fitrgp(tpoints_M1, datapointsM1, ...
    'KernelFunction', 'rationalquadratic', 'BasisFunction', bestBasisPerKernel{5});
[ypred_RQ_M1,~] = predict(gprMdl_RQ_M1,xp);

    %% ARD = automatic relevance determination
%ARD squared exponential
% Using kernel-specific best basis: bestBasisPerKernel{6}
gprMdl_ardSE_M1 = fitrgp(tpoints_M1, datapointsM1, ...
    'KernelFunction', 'ardsquaredexponential', 'BasisFunction', bestBasisPerKernel{6});
[ypred_ardSE_M1,~] = predict(gprMdl_ardSE_M1,xp);

%ARD exponential
% Using kernel-specific best basis: bestBasisPerKernel{7}
gprMdl_ardExp_M1 = fitrgp(tpoints_M1, datapointsM1, ...
    'KernelFunction', 'ardexponential', 'BasisFunction', bestBasisPerKernel{7});
[ypred_ardExp_M1,~] = predict(gprMdl_ardExp_M1,xp);

%ARD Matern 3/2
% Using kernel-specific best basis: bestBasisPerKernel{8}
gprMdl_ardM32_M1 = fitrgp(tpoints_M1, datapointsM1, ...
    'KernelFunction', 'ardmatern32', 'BasisFunction', bestBasisPerKernel{8});
[ypred_ardM32_M1,~] = predict(gprMdl_ardM32_M1,xp);


%ARD Matern 5/2
% Using kernel-specific best basis: bestBasisPerKernel{9}
gprMdl_ardM52_M1 = fitrgp(tpoints_M1, datapointsM1, ...
    'KernelFunction', 'ardmatern52', 'BasisFunction', bestBasisPerKernel{9});
[ypred_ardM52_M1,~] = predict(gprMdl_ardM52_M1,xp);

%ARD rational quadratic
% Using kernel-specific best basis: bestBasisPerKernel{10}
gprMdl_ardRQ_M1 = fitrgp(tpoints_M1, datapointsM1, ...
    'KernelFunction', 'ardrationalquadratic', 'BasisFunction', bestBasisPerKernel{10});
[ypred_ardRQ_M1,std] = predict(gprMdl_ardRQ_M1,xp);

%%To explore later %%
%non-stationary kernels
%combinging kernels through sums or products
% 
% Plot


%Squared exponental 
plot(tpoints_M1,datapointsM1,'b.');
hold on;
plot(xp,ypred_SE_M1,'g','LineWidth',1.5);
plot(xp, ypred_SE_M1+std, 'g--');
plot(xp, ypred_SE_M1-std, 'g--');
xlabel('time');
ylabel('M1 data');
legend('Data','GPR predictions (SE kernel)');
hold off
% 
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

% Store all models, their names, and their optimal basis functions
models = {
    gprMdl_SE_M1, 'Squared Exponential', bestBasisPerKernel{1};
    gprMdl_exp_M1, 'Exponential', bestBasisPerKernel{2};
    gprMdl_M32_M1, 'Matern 3/2', bestBasisPerKernel{3};
    gprMdl_M52_M1, 'Matern 5/2', bestBasisPerKernel{4};
    gprMdl_RQ_M1, 'Rational Quadratic', bestBasisPerKernel{5};
    gprMdl_ardSE_M1, 'ARD Squared Exponential', bestBasisPerKernel{6};
    gprMdl_ardExp_M1, 'ARD Exponential', bestBasisPerKernel{7};
    gprMdl_ardM32_M1, 'ARD Matern 3/2', bestBasisPerKernel{8};
    gprMdl_ardM52_M1, 'ARD Matern 5/2', bestBasisPerKernel{9};
    gprMdl_ardRQ_M1, 'ARD Rational Quadratic', bestBasisPerKernel{10}
};

% Initialize arrays to store metrics
numModels = size(models, 1);
kernelNames = cell(numModels, 1);
basisFunctions = cell(numModels, 1);
logLikelihoods = zeros(numModels, 1);
RMSEs = zeros(numModels, 1);
R2s = zeros(numModels, 1);

% Calculate metrics for each model
for i = 1:numModels
    model = models{i, 1};
    kernelNames{i} = models{i, 2};
    basisFunctions{i} = models{i, 3};
    
    % Get log likelihood (property of the model)
    logLikelihoods(i) = model.LogLikelihood;
    
    % Get predictions on training data
    yPredict = resubPredict(model);
    yActual = datapointsM1;  % Actual training data
    
    % Calculate residuals
    residuals = yActual - yPredict;
    
    % Calculate RMSE
    RMSEs(i) = sqrt(mean(residuals.^2));
    
    % Calculate R²
    SSR = sum((yPredict - mean(yActual)).^2);  % Regression sum of squares
    SST = sum((yActual - mean(yActual)).^2);   % Total sum of squares
    R2s(i) = SSR / SST;
end

% Create table with results (including basis function used)
resultsTable = table(kernelNames, basisFunctions, logLikelihoods, RMSEs, R2s, ...
    'VariableNames', {'Kernel', 'BasisFunction', 'LogLikelihood', 'RMSE', 'R2'});

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