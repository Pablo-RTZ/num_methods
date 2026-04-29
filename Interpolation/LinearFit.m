function [coeffs, E] = LinearFit(xi, fi)

%LinarFit Least squares linear regression.
%   Polynomial coefficients and error calculation for linear fitting
%
%   coeffs = LinearFit(xi, fi)
%   Returns polynomial coefficients [a0, a1] where S = a0 + a1*x
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

xi2 = xi.^2;
xifi = xi.*fi;
n = length(xi);

a0 = (sum(xi2)*sum(fi) - sum(xifi)*sum(xi)) / (n*sum(xi2) - sum(xi)^2);
a1 = (n*sum(xifi) - sum(xi)*sum(fi)) / (n*sum(xi2) - sum(xi)^2);

coeffs = [a0, a1];

E = 0;
for i = 1:n
    E = E + (fi(i) - a0 - a1*xi(i))^2;
end

% Optional outputs handling

if nargout < 2, clear E; end


end