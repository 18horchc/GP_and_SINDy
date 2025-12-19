function [t, X, dXdt, true_Xi] = generate_ground_truth(t_span, x0, params)
    % GENERATE_GROUND_TRUTH Produces clean ODE data and the true SINDy coefficients.
    %
    % Inputs:
    %   t_span - [start, end, dt] or [start, end]
    %   x0     - Initial conditions [prey_0; predator_0]
    %   params - Struct with alpha, beta, delta, gamma
    
    % 1. Define the ODE system (Lotka-Volterra)
    % dx/dt = alpha*x - beta*x*y
    % dy/dt = delta*x*y - gamma*y
    lv_ode = @(t, x) [params.alpha * x(1) - params.beta * x(1) * x(2);
                      params.delta * x(1) * x(2) - params.gamma * x(2)];

    % 2. Solve the ODE
    % Using a tight tolerance to ensure this is "Ground Truth"
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    [t, X] = ode45(lv_ode, t_span, x0, options);

    % 3. Generate the true derivatives at those points
    dXdt = zeros(size(X));
    for i = 1:length(t)
        dXdt(i, :) = lv_ode(t(i), X(i, :)');
    end

    % 4. Construct the "True Xi" Matrix
    % Assuming a library order: [1, x, y, x^2, xy, y^2]
    % This is what SINDy will try to recover.
    true_Xi = zeros(6, 2); 
    % Column 1: dx/dt
    true_Xi(2, 1) = params.alpha;   % term: x
    true_Xi(5, 1) = -params.beta;   % term: xy
    
    % Column 2: dy/dt
    true_Xi(3, 2) = -params.gamma;  % term: y
    true_Xi(5, 2) = params.delta;   % term: xy
end