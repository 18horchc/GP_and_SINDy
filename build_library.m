%% 4a. Library Function
function Theta = build_library(X, polyOrder)
    % Builds a library of functions: [1, x, y, x^2, xy, y^2, ...]
    [n, d] = size(X);
    Theta = ones(n, 1); % Start with constant term
    
    % First order terms
    Theta = [Theta, X];
    
    % Second order terms (for 2D system: x^2, xy, y^2)
    if polyOrder >= 2
        for i = 1:d
            for j = i:d
                Theta = [Theta, X(:,i).*X(:,j)];
            end
        end
    end
    % Note: You can add sin(X) or 1./X here if your model needs them
end