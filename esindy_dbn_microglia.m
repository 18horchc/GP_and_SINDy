function [t, M1, M2] = esindy_dbn_microglia(params, plot_results)
%ESINDY_DBN_MICROGLIA ESINDy+DBN model for microglial cell dynamics
%
%   Reference: Amato & Arnold (2025) "Data-driven modeling and prediction 
%              of microglial cell dynamics in the ischemic penumbra"
%              https://arxiv.org/abs/2404.10915
%              https://github.com/anarnold/Data-Driven-Modeling-Microglia
%
%   [t, M1, M2] = esindy_dbn_microglia()
%       Uses default parameters from paper and plots results
%
%   [t, M1, M2] = esindy_dbn_microglia(params)
%       Uses custom parameters specified in struct:
%           params.theta - coefficient matrix [θ1, θ2, θ3; θ4, θ5, θ6]
%                          (default: averaged coefficients from paper)
%           params.M0    - initial conditions [M1(0), M2(0)] (default: [5, 5])
%           params.t_end - end time in days (default: 50)
%
%   [t, M1, M2] = esindy_dbn_microglia(params, plot_results)
%       Set plot_results = false to suppress plotting
%
%   Model (discrete-time, equations 3.4-3.5):
%       M1(t+1) = θ1 + θ2*M1(t) + θ3*M2(t)
%       M2(t+1) = θ4 + θ5*M1(t) + θ6*M2(t)
%
%   Averaged coefficients from 10 proxy datasets:
%       θ1 = 30.5432, θ2 = 0.8960, θ3 = 0.3242
%       θ4 = 83.4946, θ5 = -0.1096, θ6 = 0.8692
%
%   Example:
%       [t, M1, M2] = esindy_dbn_microglia();

    % Default parameters (averaged coefficients from paper)
    if nargin < 1 || isempty(params)
        params = struct();
    end
    if nargin < 2
        plot_results = true;
    end
    
    % Default coefficients from eq 3.4-3.5 (averaged across 10 proxy datasets)
    default_theta = [30.5432, 0.8960, 0.3242;    % [θ1, θ2, θ3] for M1
                     83.4946, -0.1096, 0.8692];  % [θ4, θ5, θ6] for M2
    
    % Set defaults for missing fields
    if ~isfield(params, 'theta'), params.theta = default_theta; end
    if ~isfield(params, 'M0'),    params.M0 = [5; 5];           end
    if ~isfield(params, 't_end'), params.t_end = 50;            end
    
    % Extract parameters
    theta = params.theta;
    M0 = params.M0(:);  % Ensure column vector
    t_end = params.t_end;
    
    % Time points (discrete, daily)
    t = (0:1:t_end)';
    n_steps = length(t);
    
    % Initialize arrays
    M1 = zeros(n_steps, 1);
    M2 = zeros(n_steps, 1);
    
    % Initial conditions
    M1(1) = M0(1);
    M2(1) = M0(2);
    
    % Iterate the discrete-time model
    for i = 2:n_steps
        M1(i) = theta(1,1) + theta(1,2)*M1(i-1) + theta(1,3)*M2(i-1);
        M2(i) = theta(2,1) + theta(2,2)*M1(i-1) + theta(2,3)*M2(i-1);
    end
    
    % Plot results
    if plot_results
        figure;
        plot(t, M1, 'k-', 'LineWidth', 2); hold on;
        plot(t, M2, 'r-', 'LineWidth', 2);
        xlabel('Time (Days)', 'FontSize', 14);
        ylabel('cells/mm^2', 'FontSize', 14);
        title('ESINDy+DBN Model: Microglial Cell Dynamics', 'FontSize', 14);
        legend('M1', 'M2', 'Location', 'best');
        xlim([0 t_end]); 
        ylim([0 1500]);
        grid on;
        set(gca, 'FontSize', 12);
        
        % Display model info
        fprintf('\n=== ESINDy+DBN Model (Amato & Arnold, 2025) ===\n');
        fprintf('Discrete-time model:\n');
        fprintf('  M1(t+1) = %.4f + %.4f*M1(t) + %.4f*M2(t)\n', theta(1,:));
        fprintf('  M2(t+1) = %.4f + %.4f*M1(t) + %.4f*M2(t)\n', theta(2,:));
    end
end
