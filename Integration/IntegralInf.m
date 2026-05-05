function [I, dif, iter] = IntegralInf(f,method, a,sign,opts)

%IntegralInf Applies an iterative scheme to aproximate infinite integrals.
% Method must have arguments (f,a,b,n)
%
%   I = IntegralInf(f,@method, 0,"+")
%   Integrates using the method given in the interval [0,+Inf]
%
%   I = IntegralInf(f,@method,0,"-","tol",1e-8,"maxiter",200)
%   Integrates in [-Inf, 0] with the given tolerance and maximum iterations
%
%   Inputs:
%       f       - function handle
%       method  - method name
%       a       - interval start
%       sign    - Sign of the infinite
%       tol     - difference between iterations
%       maxiter - maximum number of iterations allowed
%
%   Outputs:
%       I - Integral aproximation
%       dif - difference between last two iterations
%       iter - Number of iterations needed

arguments
    f (1,1) function_handle
    method (1,1) function_handle
    a (1,1) double
    sign (1,1) string {mustBeMember(sign, ["+","-"])} = "+"
    opts.tol (1,1) double = 1e-10
    opts.maxiter (1,1) double {mustBeInteger, mustBeNonnegative} = 100
end

% Initialization

iter = 1;
dif = opts.tol +1;
I = 0;
x = 1;
left = a;

while abs(dif) > opts.tol && iter < opts.maxiter
    if sign == "+"
        right = a + x;
        dif = IntegralIter(f, method, "a", left, "b", right, "tol", 0.1*opts.tol);
        left = right;
    else
        right = a - x;
        dif = IntegralIter(f, method, "a", right, "b", left, "tol", 0.1*opts.tol);
        left = right;
    end

    I = I + dif;
    x = 2*x;
    iter = iter + 1;
end

% Optional outputs handling

if nargout < 3, clear iter; end
if nargout < 2, clear dif; end

end