function dXdt = compute_finite_differences(t, X, method)
%COMPUTE_FINITE_DIFFERENCES Compute derivatives using finite differences
%
% Inputs:
%   t: Time vector [N x 1]
%   X: State matrix [N x D] where columns are state variables
%   method: (optional) 'central' (default), 'forward', or 'backward'
%
% Outputs:
%   dXdt: Derivative matrix [N x D] where columns are derivatives

    if nargin < 3
        method = 'central';
    end
    
    [N, D] = size(X);
    dXdt = zeros(N, D);
    
    % Ensure time is sorted
    [t_sorted, sort_idx] = sort(t);
    X_sorted = X(sort_idx, :);
    
    % Compute derivatives for each state variable
    for d = 1:D
        x = X_sorted(:, d);
        
        switch lower(method)
            case 'central'
                % Central difference for interior points
                % Forward difference for first point
                % Backward difference for last point
                dXdt(1, d) = (x(2) - x(1)) / (t_sorted(2) - t_sorted(1));
                for i = 2:N-1
                    dt = t_sorted(i+1) - t_sorted(i-1);
                    if dt > 0
                        dXdt(i, d) = (x(i+1) - x(i-1)) / dt;
                    else
                        dXdt(i, d) = 0;  % Handle duplicate time points
                    end
                end
                dXdt(N, d) = (x(N) - x(N-1)) / (t_sorted(N) - t_sorted(N-1));
                
            case 'forward'
                % Forward difference
                for i = 1:N-1
                    dt = t_sorted(i+1) - t_sorted(i);
                    if dt > 0
                        dXdt(i, d) = (x(i+1) - x(i)) / dt;
                    else
                        dXdt(i, d) = 0;
                    end
                end
                % Use last forward difference for last point
                if N > 1
                    dXdt(N, d) = dXdt(N-1, d);
                else
                    dXdt(N, d) = 0;
                end
                
            case 'backward'
                % Backward difference
                % Use first backward difference for first point
                if N > 1
                    dXdt(1, d) = (x(2) - x(1)) / (t_sorted(2) - t_sorted(1));
                else
                    dXdt(1, d) = 0;
                end
                for i = 2:N
                    dt = t_sorted(i) - t_sorted(i-1);
                    if dt > 0
                        dXdt(i, d) = (x(i) - x(i-1)) / dt;
                    else
                        dXdt(i, d) = 0;
                    end
                end
                
            otherwise
                error('Unknown method: %s. Use ''central'', ''forward'', or ''backward''', method);
        end
    end
    
    % Restore original ordering if needed
    if ~issorted(t)
        [~, unsort_idx] = sort(sort_idx);
        dXdt = dXdt(unsort_idx, :);
    end
end

