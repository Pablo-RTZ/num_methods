function [S, L] = QuadraticSpline(xi, fi)

%QuadraticSpline Quadratic spline interpolation.
%   Returns a piecewise symbolic function for quadratic spline interpolation
%   with C1 continuity (continuous first derivatives)
%
%   S = QuadraticSpline(xi, fi)
%   Returns a symbolic piecewise function S(x)
%
%   Inputs:
%       xi  - independent variable data points (must be strictly increasing)
%       fi  - dependent variable data points
%
%   Outputs:
%       S   - Symbolic piecewise quadratic spline function
%       L   - Coefficient matrix where each column [a; b; c] represents
%             the quadratic: a + b*(x - xi(i)) + c*(x - xi(i))^2

arguments
    xi (:,1) double
    fi (:,1) double
end

if length(xi) ~= length(fi)
    error('xi and fi must have the same length.');
end

if length(xi) < 3
    error('At least 3 data points are required for quadratic spline interpolation.');
end

if any(diff(xi) <= 0)
    error('xi must be strictly increasing for spline interpolation.');
end

% Initialization

xi = xi(:);
fi = fi(:);
n = length(xi) - 1;

syms x
L = sym(zeros(3, n));
S = sym(0);

M = sym(zeros(2*n, 2*n));
b = zeros(2*n, 1);

% Main program

% First condition: S1'(x1) = 0 (natural spline condition at start)
M(1, 1) = 1;

% Last condition: Sn'(xn+1) = 0 (natural spline condition at end)
delta_xn = xi(n+1) - xi(n);
M(2*n, 2*n-1) = 1;
M(2*n, 2*n) = 2 * delta_xn;

% Build the system of equations
for i = 1:n
    delta_xi = xi(i+1) - xi(i);
    
    % Continuity of function values at right endpoint
    M(2*i, 2*i-1) = delta_xi;
    M(2*i, 2*i) = delta_xi^2;
    b(2*i) = fi(i+1) - fi(i);
    
    % Continuity of first derivatives at interior nodes
    if i < n
        M(2*i+1, 2*i-1) = 1;
        M(2*i+1, 2*i) = 2 * delta_xi;
        M(2*i+1, 2*i+1) = -1;
    end
end

s = M \ b;

% Construct the quadratic spline segments
for i = 1:n
    L(1, i) = fi(i);
    L(2, i) = s(2*i - 1);
    L(3, i) = s(2*i);
    
    quadratic_piece = L(1, i) + L(2, i) * (x - xi(i)) + L(3, i) * (x - xi(i))^2;
    S = piecewise(x >= xi(i) & x <= xi(i+1), quadratic_piece, S);
end

% Optional outputs handling

if nargout < 2, clear L; end

end