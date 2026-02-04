function K = periodicKernel(XN, XM, theta)
%PERIODICKERNEL Periodic kernel for GP regression.
%   K(x_a, x_b) = sigmaF^2 * exp(-(2/l^2) * sin^2(pi * |x_a - x_b| / p))
%
%   theta(1) = sigmaF (signal std)
%   theta(2) = period (p)
%   theta(3) = lengthScale (l)
%
%   K = periodicKernel(XN, XM, theta)
%   XN is n x d, XM is m x d; returns K of size n x m.

sigmaF = theta(1);
p = theta(2);
l = theta(3);

% Pairwise distance |x_a - x_b|: n x m so K(i,j) = k(XN(i,:), XM(j,:))
dist = abs(bsxfun(@minus, XN, XM.'));
K = sigmaF^2 * exp(-2 * sin(pi * dist / p).^2 / l^2);
end
