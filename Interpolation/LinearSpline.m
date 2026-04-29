function [S, L] = LinearSpline(xi, fi)

%LinearSpline Linear spline interpolation.
%   Returns a piecewise symbolic function for linear spline interpolation
%
%   S = LinearSpline(xi, fi)
%   Returns a symbolic piecewise function S(x)
%
%   Inputs:
%       xi  - independent variable data points (must be strictly increasing)
%       fi  - dependent variable data points
%
%   Outputs:
%       S   - Symbolic piecewise linear spline function
%       L   - Array of symbolic linear functions for each interval

arguments
    xi (:,1) double
    fi (:,1) double
end

if length(xi) ~= length(fi)
    error('xi and fi must have the same length.');
end

if length(xi) < 2
    error('At least 2 data points are required for linear spline interpolation.');
end

if any(diff(xi) <= 0)
    error('xi must be strictly increasing for spline interpolation.');
end

% Initalization

xi = xi(:);
fi = fi(:);
n = length(xi) - 1;

syms x
L = sym(ones(1, n));
S = sym(0);

% Main program

for i = 1:n
    m = (fi(i+1) - fi(i)) / (xi(i+1) - xi(i));
    
    L(i) = fi(i) + m * (x - xi(i));
    
    S = piecewise(x >= xi(i) & x <= xi(i+1), L(i), S);
end

% Optional outputs handling

if nargout < 2, clear L; end

end