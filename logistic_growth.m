function [t, P] = logistic_growth(params, plot_results)
%LOGISTIC_GROWTH Simulates logistic population growth
%
%   [t, P] = logistic_growth() 
%       Uses default parameters and plots results
%
%   [t, P] = logistic_growth(params)
%       Uses custom parameters specified in struct:
%           params.r     - growth rate (default: 0.5)
%           params.K     - carrying capacity (default: 1000)
%           params.P0    - initial population (default: 10)
%           params.tspan - time span [t0, tf] (default: [0, 50])
%
%   [t, P] = logistic_growth(params, plot_results)
%       Set plot_results = false to suppress plotting
%
%   Model: dP/dt = r * P * (1 - P/K)
%
%   Example:
%       [t, P] = logistic_growth();
%       params.r = 1.0; params.K = 500;
%       [t, P] = logistic_growth(params, true);

    % Default parameters
    if nargin < 1 || isempty(params)
        params = struct();
    end
    if nargin < 2
        plot_results = true;
    end
    
    % Set defaults for missing fields
    if ~isfield(params, 'r'),     params.r = 0.5;        end
    if ~isfield(params, 'K'),     params.K = 1000;       end
    if ~isfield(params, 'P0'),    params.P0 = 10;        end
    if ~isfield(params, 'tspan'), params.tspan = [0 50]; end
    
    % Extract parameters
    r = params.r;
    K = params.K;
    P0 = params.P0;
    tspan = params.tspan;
    
    % Define Logistic ODE: dP/dt = r*P*(1 - P/K)
    logisticODE = @(t, P) r * P * (1 - P / K);
    
    % Solve using ODE45
    [t, P] = ode45(logisticODE, tspan, P0);
    
    % Plot results
    if plot_results
        figure;
        plot(t, P, '-', 'LineWidth', 2);
        title('Logistic Growth (ODE45 Simulation)');
        xlabel('Time');
        ylabel('Population');
        xlim([0 50]);
        grid on;
        
        % Add carrying capacity line
        hold on;
        yline(K, '--r', 'Carrying Capacity K', 'LineWidth', 1.5);
        hold off;
        
        % Display parameters
        fprintf('Logistic Growth Parameters:\n');
        fprintf('  Growth rate r = %.4f\n', r);
        fprintf('  Carrying capacity K = %.4f\n', K);
        fprintf('  Initial population P0 = %.4f\n', P0);
    end
end
