function [S, L] = CubicSpline(xi, fi, cond)

%CubicSplineBoundary Cubic spline interpolation with specified boundary conditions.
%   Returns a piecewise symbolic function for cubic spline interpolation
%   with custom first derivative boundary conditions. If no boundary
%   conditions are given, it defaults to natural splines
%   (derivative 0 at endpoints).
%
%   S = CubicSpline(xi, fi, cond)
%   Returns a symbolic piecewise function S(x)
%
%   [S,L] = CubicSpline(xi, fi)
%   Returns a symbolic piecewise function S(x) for the natural splines, as
%   well as the different pieces separately.
%
%   Inputs:
%       xi   - independent variable data points (must be strictly increasing)
%       fi   - dependent variable data points
%       cond - vector with two elements [d0, dn] specifying the first
%              derivative values at the boundaries:
%
%   Outputs:
%       S   - Symbolic piecewise cubic spline function
%       L   - Coefficient matrix where each column [a; b; c; d] represents
%             the cubic: a + b*(x - xi(i)) + c*(x - xi(i))^2 + d*(x - xi(i))^3

arguments
    xi (:,1) double
    fi (:,1) double
    cond (1,2) double = [0,0]
end

if length(xi) ~= length(fi)
    error('xi and fi must have the same length.');
end

if length(xi) < 4
    error('At least 4 data points are required for cubic spline interpolation.');
end

if any(diff(xi) <= 0)
    error('xi must be strictly increasing for spline interpolation.');
end

% Initialization

xi = xi(:);
fi = fi(:);
n = length(xi) - 1;
h = diff(xi);

syms x
L = sym(zeros(4, n));
S = sym(0);

% Boundary conditions
d0 = cond(1);
dn = cond(2);
m = n + 1;

% Main program

% Build the tridiagonal system for second derivatives (c coefficients)
dP = zeros(1, m);
dS = zeros(1, m-1);
dI = zeros(1, m-1);
b = zeros(1, m);

% First equation
dP(1) = 2 * h(1);
dS(1) = h(1);
b(1) = 3 * ((fi(2) - fi(1)) / h(1) - d0);

% Interior equations
for i = 2:m-1
    dI(i-1) = h(i-1);
    dP(i) = 2 * (h(i-1) + h(i));
    
    if i < m-1
        dS(i) = h(i);
    end
    
    b(i) = 3 * ((fi(i+1) - fi(i)) / h(i) - (fi(i) - fi(i-1)) / h(i-1));
end

% Last equation
dI(m-1) = h(n);
dP(m) = 2 * h(n);
b(m) = 3 * (dn - (fi(n+1) - fi(n)) / h(n));

A = diag(dP) + diag(dI, -1) + diag(dS, 1);

% Solve for second derivatives
c = A \ b(:);
c = c(:)';

b_coeff = zeros(1, n);
d_coeff = zeros(1, n);

for i = 1:n
    h_i = h(i);
    c_i = c(i);
    c_next = c(i+1);
    
    b_coeff(i) = (fi(i+1) - fi(i)) / h_i - h_i / 3 * (2*c_i + c_next);
    
    d_coeff(i) = (c_next - c_i) / (3 * h_i);
end


% Each column: [a; b; c; d] for interval i
L(1, :) = fi(1:n);
L(2, :) = b_coeff;
L(3, :) = c(1:n);
L(4, :) = d_coeff;

for i = 1:n
    cubic_piece = L(1, i) + L(2, i) * (x - xi(i)) + ...
                  L(3, i) * (x - xi(i))^2 + L(4, i) * (x - xi(i))^3;
    S = piecewise(x >= xi(i) & x <= xi(i+1), cubic_piece, S);
end

% Optional outputs handling

if nargout < 2, clear L; end

end