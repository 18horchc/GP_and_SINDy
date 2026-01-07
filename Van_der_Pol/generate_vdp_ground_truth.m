function [t_true, X_true, dXdt_true, params] = generate_vdp_ground_truth(t_span, x0, mu)
%GENERATE_VDP_GROUND_TRUTH Generate ground truth solution for Van der Pol oscillator
%
% Inputs:
%   t_span: Time vector (e.g., 0:0.01:10) or [t_start, t_end] for ode45
%   x0: Initial conditions [x(0); y(0)] (default: [2; 0])
%   mu: Van der Pol parameter (default: 1.0)
%
% Outputs:
%   t_true: Time vector
%   X_true: State matrix [x(t), y(t)] of size [length(t), 2]
%   dXdt_true: Derivative matrix [dx/dt, dy/dt] of size [length(t), 2]
%   params: Structure with parameters used
%
% Van der Pol equations:
%   dx/dt = y
%   dy/dt = mu*(1 - x^2)*y - x

    if nargin < 3
        mu = 1.0;  % Default: moderate nonlinearity
    end
    if nargin < 2
        x0 = [2; 0];  % Default initial conditions
    end
    if nargin < 1
        t_span = 0:0.01:10;  % Default time span
    end
    
    % Store parameters
    params.mu = mu;
    params.x0 = x0;
    
    % If t_span is a range, use ode45; if it's a vector, interpolate
    if length(t_span) == 2
        % Use ode45 with adaptive time stepping
        options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
        [t_true, X_true] = ode45(@(t, x) vdp_ode(t, x, mu), t_span, x0, options);
        % ode45 returns [N x 2], which is what we want
    else
        % Use ode45 on full span, then interpolate to requested times
        t_full = [min(t_span), max(t_span)];
        options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
        [t_ode, X_ode] = ode45(@(t, x) vdp_ode(t, x, mu), t_full, x0, options);
        % Interpolate to requested time points
        X_true = interp1(t_ode, X_ode, t_span, 'pchip');
        t_true = t_span(:);
    end
    
    % Compute true derivatives
    dXdt_true = zeros(size(X_true));
    for i = 1:length(t_true)
        dXdt_true(i, :) = vdp_ode(t_true(i), X_true(i, :)', mu)';
    end
    
    % Ensure column vectors
    t_true = t_true(:);
end

function dxdt = vdp_ode(t, x, mu)
%VDP_ODE Van der Pol oscillator ODE
%   dx/dt = y
%   dy/dt = mu*(1 - x^2)*y - x
    
    dxdt = zeros(2, 1);
    dxdt(1) = x(2);  % dx/dt = y
    dxdt(2) = mu * (1 - x(1)^2) * x(2) - x(1);  % dy/dt = mu*(1-x^2)*y - x
end

