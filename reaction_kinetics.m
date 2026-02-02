function [t, C] = reaction_kinetics(params, plot_results)
%REACTION_KINETICS Simulates the Robertson stiff chemical reaction system
%
%   The Robertson system is a classic stiff ODE test problem:
%       dx/dt = -0.04*x + 10^4*y*z
%       dy/dt =  0.04*x - 10^4*y*z - 3*10^7*y^2
%       dz/dt =  3*10^7*y^2
%
%   [t, C] = reaction_kinetics()
%       Uses default parameters and plots results. C = [x, y, z].
%
%   [t, C] = reaction_kinetics(params)
%       Uses custom parameters specified in struct:
%           params.x0    - initial x (default: 1)
%           params.y0    - initial y (default: 0)
%           params.z0    - initial z (default: 0)
%           params.tspan - time span [t0, tf] (default: [0, 50])
%
%   [t, C] = reaction_kinetics(params, plot_results)
%       Set plot_results = false to suppress plotting
%
%   Note: Uses ode15s (stiff solver) due to large disparity in time scales.
%
%   Example:
%       [t, C] = reaction_kinetics();

    % Default parameters
    if nargin < 1 || isempty(params)
        params = struct();
    end
    if nargin < 2
        plot_results = true;
    end
    
    % Set defaults for missing fields
    if ~isfield(params, 'x0'),    params.x0 = 1;           end
    if ~isfield(params, 'y0'),    params.y0 = 0;           end
    if ~isfield(params, 'z0'),    params.z0 = 0;           end
    if ~isfield(params, 'tspan'), params.tspan = [0 50];   end
    
    % Extract parameters
    x0 = params.x0;
    y0 = params.y0;
    z0 = params.z0;
    tspan = params.tspan;
    
    % Initial conditions: [x, y, z]
    C0 = [x0; y0; z0];
    
    % Robertson system ODE:
    % dx/dt = -0.04*x + 10^4*y*z
    % dy/dt =  0.04*x - 10^4*y*z - 3*10^7*y^2
    % dz/dt =  3*10^7*y^2
    robertson_ode = @(t, C) [-0.04*C(1) + 1e4*C(2)*C(3);           % dx/dt
                              0.04*C(1) - 1e4*C(2)*C(3) - 3e7*C(2)^2;  % dy/dt
                              3e7*C(2)^2];                         % dz/dt
    
    % Solve using ode15s (stiff solver)
    [t, C] = ode15s(robertson_ode, tspan, C0);
    
    % Plot results
    if plot_results
        figure;
        plot(t, C(:,1), '-r', 'LineWidth', 2); hold on;
        plot(t, C(:,2), '-g', 'LineWidth', 2);
        plot(t, C(:,3), '-b', 'LineWidth', 2);
        xlabel('Time');
        ylabel('Concentration');
        legend('x', 'y', 'z');
        title('Robertson Stiff Reaction System');
        xlim([0 50]);
        grid on;
        
        % Display parameters
        fprintf('Robertson System Parameters:\n');
        fprintf('  Initial conditions: x(0)=%.2f, y(0)=%.2f, z(0)=%.2f\n', x0, y0, z0);
        fprintf('  Equations: dx/dt = -0.04*x + 10^4*y*z\n');
        fprintf('             dy/dt = 0.04*x - 10^4*y*z - 3*10^7*y^2\n');
        fprintf('             dz/dt = 3*10^7*y^2\n');
    end
end
