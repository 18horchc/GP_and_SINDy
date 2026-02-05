function [t, M1, M2, M_eq] = sindy_mcmc_microglia(params, plot_results)
%SINDY_MCMC_MICROGLIA SINDy+MCMC model for microglial cell dynamics
%
%   Reference: Amato & Arnold (2025) "Data-driven modeling and prediction 
%              of microglial cell dynamics in the ischemic penumbra"
%              https://arxiv.org/abs/2404.10915
%              https://github.com/anarnold/Data-Driven-Modeling-Microglia
%
%   [t, M1, M2, M_eq] = sindy_mcmc_microglia()
%       Uses default parameters from paper and plots results
%
%   [t, M1, M2, M_eq] = sindy_mcmc_microglia(params)
%       Uses custom parameters specified in struct:
%           params.theta - coefficient vector [θ1, θ2, θ3, θ4, θ5, θ6]
%                          (default: MCMC posterior means from paper)
%           params.M0    - initial conditions [M1(0), M2(0)] (default: [5, 5])
%           params.tspan - time span [t0, tf] (default: [0, 50])
%
%   [t, M1, M2, M_eq] = sindy_mcmc_microglia(params, plot_results)
%       Set plot_results = false to suppress plotting
%
%   Model (continuous-time ODE, equations 3.1-3.3):
%       dM1/dt = θ1 + θ2*M1 + θ3*M2
%       dM2/dt = θ4 + θ5*M1 + θ6*M2
%
%   MCMC posterior mean estimates:
%       θ1 = 32.9650, θ2 = -0.1377, θ3 = 0.3150
%       θ4 = 78.5130, θ5 = -0.1550, θ6 = 0.0217
%
%   Example:
%       [t, M1, M2] = sindy_mcmc_microglia();

    % Default parameters (MCMC posterior means from paper)
    if nargin < 1 || isempty(params)
        params = struct();
    end
    if nargin < 2
        plot_results = true;
    end
    
    % Default coefficients from eq 3.3 (MCMC posterior means)
    default_theta = [32.9650, -0.1377, 0.3150, ...  % θ1, θ2, θ3 (M1 equation)
                     78.5130, -0.1550, 0.0217];     % θ4, θ5, θ6 (M2 equation)
    
    % Set defaults for missing fields
    if ~isfield(params, 'theta'), params.theta = default_theta; end
    if ~isfield(params, 'M0'),    params.M0 = [5; 5];           end
    if ~isfield(params, 'tspan'), params.tspan = [0 50];        end
    
    % Extract parameters
    theta = params.theta;
    M0 = params.M0(:);  % Ensure column vector
    tspan = params.tspan;
    
    % Build coefficient matrix A and constant vector b
    % dM/dt = A*M + b
    A = [theta(2), theta(3);    % [θ2, θ3]
         theta(5), theta(6)];   % [θ5, θ6]
    b = [theta(1); theta(4)];   % [θ1; θ4]
    
    % Define ODE
    sindy_ode = @(t, M) A * M + b;
    
    % Solve ODE
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    [t, M] = ode45(sindy_ode, tspan, M0, opts);
    
    M1 = M(:, 1);
    M2 = M(:, 2);
    
    % Compute equilibrium point: A*M + b = 0  =>  M = -A\b
    M_eq = -A \ b;
    
    % Plot results
    if plot_results
        figure;
        plot(t, M1, 'k-', 'LineWidth', 2); hold on;
        plot(t, M2, 'r-', 'LineWidth', 2);
        xlabel('Time (Days)', 'FontSize', 14);
        ylabel('cells/mm^2', 'FontSize', 14);
        title('SINDy+MCMC Model: Microglial Cell Dynamics', 'FontSize', 14);
        legend('M1', 'M2', 'Location', 'best');
        xlim([0 50]); 
        ylim([0 1500]);
        grid on;
        set(gca, 'FontSize', 12);
        
        % Display model info
        fprintf('\n=== SINDy+MCMC Model (Amato & Arnold, 2025) ===\n');
        fprintf('Continuous-time ODE:\n');
        fprintf('  dM1/dt = %.4f + %.4f*M1 + %.4f*M2\n', theta(1), theta(2), theta(3));
        fprintf('  dM2/dt = %.4f + %.4f*M1 + %.4f*M2\n', theta(4), theta(5), theta(6));
        fprintf('Equilibrium: M1* = %.2f, M2* = %.2f cells/mm^2\n', M_eq(1), M_eq(2));
    end
end
