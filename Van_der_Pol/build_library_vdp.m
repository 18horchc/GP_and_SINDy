function [Theta, library_names] = build_library_vdp(X, polyOrder)
%BUILD_LIBRARY_VDP Build polynomial library matrix for Van der Pol oscillator
%
% Inputs:
%   X: State matrix [N x 2] where columns are [x, y]
%   polyOrder: Maximum polynomial order (default: 3)
%              Order 2: [1, x, y, x², xy, y²]
%              Order 3: [1, x, y, x², xy, y², x³, x²y, xy², y³]
%
% Outputs:
%   Theta: Library matrix [N x num_terms]
%   library_names: Cell array of strings describing each term

    if nargin < 2
        polyOrder = 3;  % Default: up to cubic terms
    end
    
    [N, D] = size(X);
    if D ~= 2
        error('Van der Pol requires 2D state (x, y)');
    end
    
    x = X(:, 1);
    y = X(:, 2);
    
    Theta = [];
    library_names = {};
    
    % Order 0: Constant
    Theta = [Theta, ones(N, 1)];
    library_names{end+1} = '1';
    
    % Order 1: Linear terms
    Theta = [Theta, x];
    library_names{end+1} = 'x';
    Theta = [Theta, y];
    library_names{end+1} = 'y';
    
    if polyOrder >= 2
        % Order 2: Quadratic terms
        Theta = [Theta, x.^2];
        library_names{end+1} = 'x^2';
        Theta = [Theta, x.*y];
        library_names{end+1} = 'xy';
        Theta = [Theta, y.^2];
        library_names{end+1} = 'y^2';
    end
    
    if polyOrder >= 3
        % Order 3: Cubic terms
        Theta = [Theta, x.^3];
        library_names{end+1} = 'x^3';
        Theta = [Theta, x.^2.*y];
        library_names{end+1} = 'x^2y';
        Theta = [Theta, x.*y.^2];
        library_names{end+1} = 'xy^2';
        Theta = [Theta, y.^3];
        library_names{end+1} = 'y^3';
    end
    
    library_names = library_names(:);  % Ensure column vector
end

