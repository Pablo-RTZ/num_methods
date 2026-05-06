function [X,Y] = FiniteDifferences(p, q, r, a, b, bc, n)

%FiniteDifferences Finite differences method for solving second order linear BVP.
% Solves problems of the form y'' + p(x)y' + q(x)y = r(x)
%
%   [X, Y] = FiniteDifferences(p, q, r, a, b, bc, n)
%   Solves the BVP and returns function at evaluation nodes
%
%   Inputs:
%       p   - Coefficient function for y' term
%       q   - Coefficient function for y term
%       r   - Forcing function
%       a   - Interval start
%       b   - Interval end
%       bc  - Boundary conditions [y(a), y(b)]
%       n   - Number of subintervals (optional, default=10)

arguments
    p (:,1) function_handle
    q (:,1) function_handle
    r (:,1) function_handle
    a (1,1) double
    b (1,1) double {mustBeGreaterThan(b,a)}
    bc (2,1) double
    n (1,1) double {mustBeInteger, mustBePositive} = 10
end

% Initialization

alpha = bc(1);
beta = bc(2);

h = (b-a)/n;
X = a:h:b;
X = X(:);
x = X(2:end-1);

px = p(x);
qx = q(x);
rx = r(x);

% Main program

dp = 2+h^2*qx;
ds = -1+h/2*px(1:end-1);
di = -1-h/2*px(2:end);
d = -h^2*rx;
d(1) = d(1) + (1+h/2*px(1))*alpha;
d(end) = d(end) + (1-h/2*px(end))*beta;

A = diag(dp) + diag(ds,1) + diag(di,-1);
y = A\d;
Y = [alpha; y; beta];

end