function [x, y] = LinearShooting(p, q, r, a, b, bc, n)

%LinearShooting Linear shooting method for solving second order linear BVP.
% Solves problems of the form y'' + p(x)y' + q(x)y = r(x)
%
%   [X, Y] = LinearShooting(p, q, r, a, b, bc, n)
%   Solves the BVP and returns function at evaluation nodes
%
%   Inputs:
%       p   - Coefficient function for y' term
%       q   - Coefficient function for y term  
%       r   - Forcing function
%       a   - Interval start
%       b   - Interval end
%       bc  - Boundary conditions [y(a), y(b)]
%       n   - Number of subintervals

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
h = (b-a)/n;
x = a:h:b;
alpha = bc(1);
beta = bc(2);

ode1 = @(x, y) [y(2); r(x) - p(x)*y(2) - q(x)*y(1)];
ode2 = @(x, y) [y(2); -p(x)*y(2) - q(x)*y(1)];

y0_1 = [alpha; 0];
y0_2 = [0; 1];

[x, y1] = ode45(ode1, x, y0_1);
[x, y2] = ode45(ode2, x, y0_2);

C = (beta - y1(end,1)) / y2(end,1);
y = y1(:,1) + C * y2(:,1);

end