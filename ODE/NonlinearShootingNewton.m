function [x, y, t, iter, dif] = NonlinearShootingNewton(f, a, b, bc, opts)

%NonlinearShootingNewton Nonlinear shooting method with Newton's method for solving second order nonlinear BVP.
% Solves problems of the form y'' = f(x, y, y') with Dirichlet boundary conditions
%
%   [x, y, t, iter, incr] = NonlinearShootingNewton(f, a, b, bc)
%   Solves the BVP using Newton's method and returns function at evaluation nodes
%
%   Inputs:
%       f       - Function handle for the nonlinear ODE system
%       a       - Interval start
%       b       - Interval end
%       bc      - Boundary conditions [y(a), y(b)]
%       n       - Number of subintervals
%       tol     - Tolerance for convergence
%       maxiter - Maximum number of iterations
%
%   Outputs:
%       x       - Evaluation nodes
%       y       - Solution values at nodes
%       t       - Final initial slope estimate
%       iter    - Number of iterations performed
%       dif    - Final increment value

arguments
    f (:,1) function_handle
    a (1,1) double
    b (1,1) double {mustBeGreaterThan(b,a)}
    bc (2,1) double
    opts.n (1,1) double {mustBeInteger, mustBePositive} = 10
    opts.tol (1,1) double = 1e-6
    opts.maxiter (1,1) double {mustBeInteger, mustBePositive} = 100
end

% Initialization

h = (b-a)/opts.n;
x = a:h:b;
alpha = bc(1);
beta = bc(2);
iter = 1;

% Initial guess for slope
t = (beta - alpha)/(b - a);

[x, y] = ode45(f, x, [alpha, t, 0, 1]);
dif = abs(y(end,1) - beta);

while dif > opts.tol && iter < opts.maxiter

    t = t - (y(end,1) - beta) / y(end,3);
    [x, y] = ode45(f, x, [alpha, t, 0, 1]);
    dif = abs(y(end,1) - beta);
    iter = iter + 1;
end

if iter == opts.maxiter
    error("Maximum number of iterations reached without convergence");
end

end
