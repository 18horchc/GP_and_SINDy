function [Xi_true, library_names] = compute_true_vdp_coefficients(mu, polyOrder)
%COMPUTE_TRUE_VDP_COEFFICIENTS Compute true Van der Pol coefficients for library
%
% Inputs:
%   mu: Van der Pol parameter (default: 1.0)
%   polyOrder: Polynomial order of library (default: 3)
%
% Outputs:
%   Xi_true: True coefficient matrix [num_terms x 2]
%   library_names: Library term names (for reference)

    if nargin < 2
        polyOrder = 3;
    end
    if nargin < 1
        mu = 1.0;
    end
    
    % Build library names (same as build_library_vdp)
    library_names = {};
    library_names{end+1} = '1';
    library_names{end+1} = 'x';
    library_names{end+1} = 'y';
    
    if polyOrder >= 2
        library_names{end+1} = 'x^2';
        library_names{end+1} = 'xy';
        library_names{end+1} = 'y^2';
    end
    
    if polyOrder >= 3
        library_names{end+1} = 'x^3';
        library_names{end+1} = 'x^2y';
        library_names{end+1} = 'xy^2';
        library_names{end+1} = 'y^3';
    end
    
    library_names = library_names(:);
    num_terms = length(library_names);
    
    % Initialize coefficient matrix
    Xi_true = zeros(num_terms, 2);
    
    % Van der Pol equations:
    % dx/dt = y
    % dy/dt = mu*(1 - x^2)*y - x = mu*y - mu*x^2*y - x
    
    % Find indices for each term
    idx_1 = find(strcmp(library_names, '1'));
    idx_x = find(strcmp(library_names, 'x'));
    idx_y = find(strcmp(library_names, 'y'));
    idx_x2y = find(strcmp(library_names, 'x^2y'));
    
    % dx/dt = y
    if ~isempty(idx_y)
        Xi_true(idx_y, 1) = 1.0;
    end
    
    % dy/dt = mu*y - mu*x^2*y - x
    if ~isempty(idx_x)
        Xi_true(idx_x, 2) = -1.0;  % -x term
    end
    if ~isempty(idx_y)
        Xi_true(idx_y, 2) = mu;  % mu*y term
    end
    if ~isempty(idx_x2y)
        Xi_true(idx_x2y, 2) = -mu;  % -mu*x^2*y term
    end
end

