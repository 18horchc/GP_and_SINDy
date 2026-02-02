function [t, prey, predators] = lotka_volterra_model(params, plot_results)
%LOTKA_VOLTERRA_MODEL Simulates the Lotka-Volterra predator-prey dynamics
%
%   [t, prey, predators] = lotka_volterra_model()
%       Uses default parameters and plots results
%
%   [t, prey, predators] = lotka_volterra_model(params)
%       Uses custom parameters specified in struct:
%           params.alpha - prey growth rate (default: 1.1)
%           params.beta  - predation rate (default: 0.4)
%           params.delta - predator growth rate per prey eaten (default: 0.4)
%           params.gamma - predator death rate (default: 0.4)
%           params.x0    - initial prey population (default: 20)
%           params.y0    - initial predator population (default: 5)
%           params.tspan - time span [t0, tf] (default: [0, 50])
%
%   [t, prey, predators] = lotka_volterra_model(params, plot_results)
%       Set plot_results = false to suppress plotting
%
%   Model:
%       dx/dt = alpha*x - beta*x*y   (prey)
%       dy/dt = delta*x*y - gamma*y  (predator)
%
%   Example:
%       [t, prey, pred] = lotka_volterra_model();
%       params.alpha = 1.5; params.beta = 0.5;
%       [t, prey, pred] = lotka_volterra_model(params, true);

    % Default parameters
    if nargin < 1 || isempty(params)
        params = struct();
    end
    if nargin < 2
        plot_results = true;
    end
    
    % Set defaults for missing fields
    if ~isfield(params, 'alpha'), params.alpha = 1.1;      end
    if ~isfield(params, 'beta'),  params.beta = 0.4;       end
    if ~isfield(params, 'delta'), params.delta = 0.4;      end
    if ~isfield(params, 'gamma'), params.gamma = 0.4;      end
    if ~isfield(params, 'x0'),    params.x0 = 20;          end
    if ~isfield(params, 'y0'),    params.y0 = 5;           end
    if ~isfield(params, 'tspan'), params.tspan = [0 50];   end
    
    % Extract parameters
    alpha = params.alpha;
    beta = params.beta;
    delta = params.delta;
    gamma = params.gamma;
    P0 = [params.x0; params.y0];
    tspan = params.tspan;
    
    % Define Lotka-Volterra ODE system
    lv_ode = @(t, P) [alpha*P(1) - beta*P(1)*P(2);      % dx/dt (prey)
                      delta*P(1)*P(2) - gamma*P(2)];    % dy/dt (predator)
    
    % Solve using ODE45
    [t, P] = ode45(lv_ode, tspan, P0);
    
    % Extract populations
    prey = P(:, 1);
    predators = P(:, 2);
    
    % Plot results
    if plot_results
        figure;
        
        % Plot population dynamics over time
        subplot(2,1,1);
        plot(t, prey, 'b-', 'LineWidth', 2);
        hold on;
        plot(t, predators, 'r-', 'LineWidth', 2);
        xlabel('Time');
        ylabel('Population');
        legend('Prey (x)', 'Predators (y)');
        title('Lotka-Volterra Population Dynamics Over Time');
        grid on;
        
        % Plot the phase space (predators vs. prey)
        subplot(2,1,2);
        plot(prey, predators, 'k-', 'LineWidth', 2);
        xlabel('Prey Population');
        ylabel('Predator Population');
        title('Phase Space Plot (Predators vs. Prey)');
        grid on;
        
        % Display parameters
        fprintf('Lotka-Volterra Parameters:\n');
        fprintf('  Prey growth rate (alpha) = %.4f\n', alpha);
        fprintf('  Predation rate (beta) = %.4f\n', beta);
        fprintf('  Predator growth rate (delta) = %.4f\n', delta);
        fprintf('  Predator death rate (gamma) = %.4f\n', gamma);
        fprintf('  Initial prey (x0) = %.4f\n', params.x0);
        fprintf('  Initial predators (y0) = %.4f\n', params.y0);
    end
end
