function [sol, dif, iter, ACOC] = Traub(f, df, opts)

% Traub Traub's method for solving nonlinear equations.
% Uses function handles, and uses step size as stopping criterion.
% Traub's method is a two-step Newton-like method with third-order convergence.
%
%   sol = Traub(f, df)
%   Uses 0 as initial guess, a 1e-10 tolerance, and 50 iterations by
%   default
%
%   sol = Traub(f, df, "x0", 1, "tol", 1e-8, "maxiter", 100)
%   allows name-value pair inputs in any order.
%
%   Outputs:
%       sol   - root approximation
%       dif   - last step size
%       iter  - number of iterations
%       ACOC  - Approximate Computational Order of Convergence

arguments
    f (1,1) function_handle
    df (1,1) function_handle

    opts.x0 (1,1) double = 0
    opts.tol (1,1) double = 1e-10
    opts.maxiter (1,1) double = 50
end

% Initialization

x0 = opts.x0;
tol = opts.tol;
maxiter = opts.maxiter;

iter = 0;
dif = tol + 1;
I = [];

% Main loop

while (dif > tol) && (iter < maxiter)
    y0 = x0 - f(x0)/df(x0);
    x1 = y0 - f(y0)/df(x0);
    dif = abs(x1 - x0);
    I(end+1) = dif;
    x0 = x1;
    iter = iter + 1;
end

% Stopping criterion

if iter >= maxiter
    sol = [];
    ACOC = [];
    disp("The method has not converged to a root")
else
    sol = x1;
    
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