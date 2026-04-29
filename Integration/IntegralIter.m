function [I,dif,iter] = IntegralIter(f,method,opts)

%IntegralIter Applies an iterative scheme for any integral method.
% Method must have arguments (f,a,b,n) or (f,n)
%
%   I = IntegralIter(f,method,"a",-1,"b",1)
%   Integrates using the method given in the finite interval
%
%   I = IntegralIter(f,method,"tol",1e-8,"maxiter",200)
%   Integrates in infinite intervals (and sets tolerance and maxiter)
%
%   Inputs:
%       f       - function handle
%       method  - method name
%       a       - interval start
%       b       - interval end
%       tol     - difference between iterations
%       maxiter - maximum number of iterations allowed
%
%   Outputs:
%       I - Integral aproximation
%       dif - difference between iterations
%       iter - Number of iterations needed

arguments
    f (1,1) function_handle
    method (1,1) str
    opts.a (1,1) double
    opts.b (1,1) double
    opts.tol (1,1) double = 1e-10
    opts.maxiter (1,1) int = 100
end

filename = method + ".m";
if ~isfile(filename)
    error("Method not found on current directory");
end

hasInterval = isfield(opts,"a") && isfield(opts,"b");

if hasInterval && ~(opts.a < opts.b)
    error('a must be less than b.');
end

% Initialization

dif=tol+1;
iter=1;
n=1;

% Main loop

while dif > opts.tol && iter < opts.maxiter
    n = 2*n;
    iter = iter + 1;

    if hasInterval
        I_new = method(f, opts.a, opts.b, n);
    else
        I_new = method(f, n);
    end

    if iter > 2
        dif = abs(I_new - I);
    end

    I = I_new;
end

% Stopping criterion

if dif>tol
    disp("The method has not converged within the required tolerance")
end

% Optional outputs handling

if nargout < 3, clear iter; end
if nargout < 2, clear dif; end

end