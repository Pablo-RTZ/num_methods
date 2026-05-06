function [x, y, t, iter, dif] = NonlinearShootingSecant(f, a, b, bc, opts)
%NonlinearShootingSecant Nonlinear shooting method with Secant method for solving second order nonlinear BVP.
% Solves problems of the form y'' = f(x, y, y') with Dirichlet boundary conditions
%
%   [x, y, t, iter, incr] = NonlinearShootingSecant(f, a, b, bc, n, tol, maxiter)
%   Solves the BVP using the Secant method and returns function at evaluation nodes
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

% Initial guesses for slope (t1 and t2)
t1 = (beta - alpha)/(b - a);
t2 = t1 + 0.1;

[x, y1] = ode45(f, x, [alpha, t1, 0, 1]);
[x, y2] = ode45(f, x, [alpha, t2, 0, 1]);

dif = abs(y2(end,1) - beta);

% Main program

while dif > opts.tol && iter < opts.maxiter
    t = t2 - (y2(end,1) - beta) * (t2 - t1) / (y2(end,1) - y1(end,1));
    
    t1 = t2;
    y1 = y2;
    t2 = t;
    
    [x, y2] = ode45(f, x, [alpha, t2, 0, 1]);
    
    dif = abs(y2(end,1) - beta);
    iter = iter + 1;
end

y = y2(:,1);

t = t2;

if iter == opts.maxiter && dif > opts.tol
    error("Maximum number of iterations reached without convergence");
end

end