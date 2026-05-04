function [sol, dif, iter, ACOC] = Secant(f, opts)

%Secant Secant method for solving nonlinear equations.
% Uses function handles, and uses step size as stopping criterion.
% This is a derivative-free method that approximates the derivative using previous points.
%
%   sol = Secant(f)
%   Uses 0 and 1 as initial guesses, a 1e-10 tolerance, and 50 iterations by
%   default
%
%   sol = Secant(f, "x0", 0, "x1", 2, "tol", 1e-8, "maxiter", 100)
%   allows name-value pair inputs in any order.
%
%   Outputs:
%       sol   - root approximation
%       dif   - last step size
%       iter  - number of iterations
%       ACOC  - Approximate Computational Order of Convergence

arguments
    f (1,1) function_handle

    opts.x0 (1,1) {mustBeA(opts.x0, {'double','sym'})} = 0
    opts.x1 (1,1) double = 1
    opts.tol (1,1) double = 1e-10
    opts.maxiter (1,1) double = 50
end

% Initialization

x0 = opts.x0;
x1 = opts.x1;
tol = opts.tol;
maxiter = opts.maxiter;

iter = 0;
dif = tol + 1;
I = [];

% Main loop

while (dif > tol) && (iter < maxiter)
    x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
    dif = abs(x2 - x1);
    I(end+1) = dif;
    x0 = x1;
    x1 = x2;
    iter = iter + 1;
end

% Stopping criterion

if iter >= maxiter
    sol = [];
    ACOC = [];
    disp("The method has not converged to a root")
else
    sol = x2;
    
    if numel(I) >= 3
        ACOC = log(I(3:end)./I(2:end-1)) ./ log(I(2:end-1)./I(1:end-2));
    else
        ACOC = [];
    end
end

% Optional outputs handling

if nargout < 4, clear ACOC; end
if nargout < 3, clear iter; end
if nargout < 2, clear dif; end

end