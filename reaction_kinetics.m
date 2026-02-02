function [t, C] = reaction_kinetics(params, plot_results)
%REACTION_KINETICS Simulates stiff chemical reaction kinetics A -> B -> C
%
%   [t, C] = reaction_kinetics()
%       Uses default parameters and plots results
%
%   [t, C] = reaction_kinetics(params)
%       Uses custom parameters specified in struct:
%           params.k1    - fast reaction rate A->B (default: 1e3)
%           params.k2    - slow reaction rate B->C (default: 1e-2)
%           params.C0    - initial concentrations [A0, B0, C0] (default: [1, 0, 0])
%           params.tspan - time span [t0, tf] (default: [0, 100])
%
%   [t, C] = reaction_kinetics(params, plot_results)
%       Set plot_results = false to suppress plotting
%
%   Model (stiff system):
%       dA/dt = -k1 * A          (fast decay)
%       dB/dt = k1 * A - k2 * B  (intermediate)
%       dC/dt = k2 * B           (slow formation)
%
%   Note: Uses ode15s for stiff systems due to large difference in k1/k2
%
%   Example:
%       [t, C] = reaction_kinetics();
%       params.k1 = 1e4; params.k2 = 1e-1;
%       [t, C] = reaction_kinetics(params, true);

    % Default parameters
    if nargin < 1 || isempty(params)
        params = struct();
    end
    if nargin < 2
        plot_results = true;
    end
    
    % Set defaults for missing fields
    if ~isfield(params, 'k1'),    params.k1 = 1e3;         end
    if ~isfield(params, 'k2'),    params.k2 = 1e-2;        end
    if ~isfield(params, 'C0'),    params.C0 = [1.0, 0, 0]; end
    if ~isfield(params, 'tspan'), params.tspan = [0 100];  end
    
    % Extract parameters
    k1 = params.k1;
    k2 = params.k2;
    C0 = params.C0;
    tspan = params.tspan;
    
    % Define stiff reaction kinetics ODE system
    stiff_ode = @(t, C) [-k1 * C(1);              % dA/dt
                          k1 * C(1) - k2 * C(2);  % dB/dt
                          k2 * C(2)];             % dC/dt
    
    % Solve using ode15s (stiff solver)
    [t, C] = ode15s(stiff_ode, tspan, C0);
    
    % Plot results
    if plot_results
        figure;
        plot(t, C(:,1), '-r', 'LineWidth', 2); hold on;
        plot(t, C(:,2), '-g', 'LineWidth', 2);
        plot(t, C(:,3), '-b', 'LineWidth', 2);
        xlabel('Time');
        ylabel('Concentration');
        legend('A', 'B', 'C');
        title('Stiff Reaction Kinetics: A \rightarrow B \rightarrow C');
        grid on;
        
        % Display parameters
        fprintf('Reaction Kinetics Parameters:\n');
        fprintf('  Fast rate k1 (A->B) = %.2e\n', k1);
        fprintf('  Slow rate k2 (B->C) = %.2e\n', k2);
        fprintf('  Stiffness ratio k1/k2 = %.2e\n', k1/k2);
        fprintf('  Initial concentrations: A0=%.2f, B0=%.2f, C0=%.2f\n', C0(1), C0(2), C0(3));
    end
end
