function [coeffs, E] = LinearFit(xi, fi)

%LinarFit Least squares linear regression.
%   Polynomial coefficients and error calculation for linear fitting
%
%   coeffs = LinearFit(xi, fi)
%   Returns polynomial coefficients [a1, a0] where S = a0 + a1*x
%
%   Inputs:
%       xi  - independent variable data points
%       fi  - dependent variable data points
%
%   Outputs:
%       coeffs  - Coefficients for the line
%       E       - Quadratic error at the fitting nodes

arguments
    xi (:,1) double
    fi (:,1) double
end

if length(xi) ~= length(fi)
    error('xi and fi must have the same length.');
end

if length(xi) < 2
    error('At least 2 data points are required for linear regression.');
end

% Main program

n = length(xi);
a1 = (n*sum(xi.*fi) - sum(xi)*sum(fi)) / (n*sum(xi.^2) - sum(xi)^2);
a0 = (sum(fi) - a1*sum(xi)) / n;

coeffs = [a1, a0];

E = sum((fi - (a0 + a1*xi)).^2);

% Optional outputs handling

if nargout < 2, clear E; end


end