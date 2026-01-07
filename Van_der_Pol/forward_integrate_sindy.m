function [t_sim, X_sim] = forward_integrate_sindy(Xi, library_names, t_span, x0, options)
%FORWARD_INTEGRATE_SINDY Forward integrate discovered SINDy model
%
% Inputs:
%   Xi: Coefficient matrix [num_terms x num_states]
%   library_names: Cell array of library term names
%   t_span: Time span [t_start, t_end] or time vector
%   x0: Initial conditions [x0; y0]
%   options: (optional) ODE solver options
%
% Outputs:
%   t_sim: Time vector
%   X_sim: Simulated state matrix [length(t) x num_states]

    if nargin < 5
        options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
    end
    
    % Create ODE function handle
    sindy_ode = @(t, x) sindy_rhs(t, x, Xi, library_names);
    
    % Integrate
    if length(t_span) == 2
        % Time span provided
        [t_sim, X_sim] = ode45(sindy_ode, t_span, x0, options);
    else
        % Time vector provided
        [t_sim, X_sim] = ode45(sindy_ode, [min(t_span), max(t_span)], x0, options);
        % Interpolate to requested time points
        X_sim = interp1(t_sim, X_sim, t_span, 'pchip');
        t_sim = t_span(:);
    end
end

function dxdt = sindy_rhs(t, x, Xi, library_names)
%SINDY_RHS Right-hand side function for SINDy-discovered model
%   Evaluates: dx/dt = Theta(x) * Xi
    
    % Determine polynomial order from library names
    % Count terms to infer order
    num_terms = length(library_names);
    if num_terms <= 3
        polyOrder = 1;
    elseif num_terms <= 6
        polyOrder = 2;
    else
        polyOrder = 3;
    end
    
    % Build library at current state (x is column vector [x; y])
    Theta_x = build_library_vdp(x', polyOrder);
    
    % Compute derivative
    dxdt = Theta_x * Xi;
    dxdt = dxdt(:);  % Ensure column vector
end

