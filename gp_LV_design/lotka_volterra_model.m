function [t, prey, predators] = lotka_volterra_model(params, plot_results)
%LOTKA_VOLTERRA_MODEL Simulates the Lotka-Volterra predator-prey dynamics
%
%       dx/dt = alpha*x - beta*x*y   (prey: growth - predation)
%       dy/dt = delta*beta*x*y - gamma*y  (predator: conversion of prey eaten - mortality)
%
%   Here delta is the predator food conversion efficiency (prey -> predator growth).
%   Prey consumed per unit time = beta*x*y; predator growth = delta * (prey consumed).
%
%   [t, prey, predators] = lotka_volterra_model()
%   [t, prey, predators] = lotka_volterra_model(params)
%       Uses custom parameters specified in struct:
%           params.alpha - prey growth rate (default: 1)
%           params.beta  - predation rate / search efficiency (default: 0.2)
%           params.delta - predator growth rate per prey eaten (conversion efficiency) (default: 0.5)
%           params.gamma - predator death rate (default: 0.2)
%           params.x0    - initial prey population (default: 1)
%           params.y0    - initial predator population (default: 2)
%           params.tspan - time span [t0, tf] (default: [0, 50])
%
%   [t, prey, predators] = lotka_volterra_model(params, plot_results)
%       Set plot_results = false to suppress plotting
%   Model:
%       dx/dt = alpha*x - beta*x*y   (prey)
%       dy/dt = delta*x*y - gamma*y  (predator)%   
%Example:
%       [t, prey, pred] = lotka_volterra_model();

    % Default parameters 
    if nargin < 1 || isempty(params)
        params = struct();
    end
    if nargin < 2
        plot_results = true;
    end
    
    % Set defaults for missing fields
    if ~isfield(params, 'alpha'), params.alpha = 1;       end
    if ~isfield(params, 'beta'),  params.beta = 0.2;      end
    if ~isfield(params, 'delta'), params.delta = 0.5;     end
    if ~isfield(params, 'gamma'), params.gamma = 0.2;     end
    if ~isfield(params, 'x0'),    params.x0 = 1;          end
    if ~isfield(params, 'y0'),    params.y0 = 2;          end
    if ~isfield(params, 'tspan'), params.tspan = [0 50];  end
    
    % Extract parameters
    alpha = params.alpha;
    beta = params.beta;
    delta = params.delta;
    gamma = params.gamma;
    P0 = [params.x0; params.y0];
    tspan = params.tspan;
    
    % Define Lotka-Volterra ODE system 
    % dy/dt = delta*beta*x*y - gamma*y  (predator growth = conversion efficiency * prey consumed)
    lv_ode = @(t, P) [alpha*P(1) - beta*P(1)*P(2);           % dx/dt (prey)
                      delta*beta*P(1)*P(2) - gamma*P(2)];   % dy/dt (predator)
    
    % Solve using ODE45
    [t, P] = ode45(lv_ode, tspan, P0);
    
    % Extract populations
    prey = P(:, 1);
    predators = P(:, 2);
    
    % Plot results
    if plot_results
        figure;
        plot(t, prey, 'b-', 'LineWidth', 2);
        hold on;
        plot(t, predators, 'r-', 'LineWidth', 2);
        xlabel('Time');
        ylabel('Population');
        legend('Prey (x)', 'Predators (y)');
        title('Lotka-Volterra Population Dynamics Over Time');
        xlim([0 50]);
        grid on;
        
        % Display parameters
        fprintf('Lotka-Volterra Parameters:\n');
        fprintf('  Prey growth rate (alpha) = %.4f\n', alpha);
        fprintf('  Predation rate (beta) = %.4f\n', beta);
        fprintf('  Predator growth rate (delta) = %.4f\n', delta);
        fprintf('  Predator death rate (gamma) = %.4f\n', gamma);
        fprintf('  Initial prey (x0) = %.4f, predators (y0) = %.4f\n', params.x0, params.y0);
    end
end
