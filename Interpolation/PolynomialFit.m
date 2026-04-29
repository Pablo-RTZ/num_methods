function [coeffs, E] = PolynomialFit(xi, fi)

%PolynomialFit Least squares polynomial regression.
%   Polynomial coefficients and error calculation for polynomial fitting
%
%   coeffs = PolynomialFit(xi, fi)
%   Returns polynomial coefficients [a0, a1, ..., an] where 
%   S = a0 + a1*x + a2*x^2 + ... + an*x^n
%
%   Inputs:
%       xi  - independent variable data points
%       fi  - dependent variable data points
%
%   Outputs:
%       coeffs  - Coefficients for the polynomial
%       E       - Quadratic error at the fitting nodes

arguments
    xi (:,1) double
    fi (:,1) double
end

if length(xi) ~= length(fi)
    error('xi and fi must have the same length.');
end

if length(xi) < 2
    error('At least 2 data points are required for polynomial regression.');
end

% Main program
xi = xi(:);
fi = fi(:);
n = length(xi);

exp = 0:n-1;
M = xi.^exp;


a = M \ fi;
coeffs = a(:)';

E = sum((fi - M * a).^2);

% Optional outputs handling

if nargout < 2, clear E; end

end