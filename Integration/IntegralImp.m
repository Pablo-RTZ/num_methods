function [I, dif, iter] = IntegralImp(f,method, a, b, sign,opts)

%IntegralImp Applies an iterative scheme to aproximate improper integrals.
% Method must have arguments (f,a,b,n)
%
%   I = IntegralInf(f,@method, 0,1,"+")
%   Integrates using the method given in the interval [0,1] with
%   discontinuity at 1
%
%   I = IntegralInf(f,@method,0,1,"-","tol",1e-8,"maxiter",200)
%   Integrates in [0,1] with discontinuity at 0 with the given tolerance
%   and maximum iterations
%
%   Inputs:
%       f       - function handle
%       method  - method name
%       a       - interval start
%       b       - interval end
%       sign    - Side with the discontinuity
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
    b (1,1) double
    sign (1,1) string {mustBeMember(sign, ["+","-"])} = "+"
    opts.tol (1,1) double = 1e-2
    opts.maxiter (1,1) double {mustBeInteger, mustBeNonnegative} = 20
end

if ~(a < b)
    error('a must be less than b.');
end

% Initialization

iter = 1;
dif = opts.tol +1;
Iprev = 0;
x = 10;

while abs(dif) > opts.tol && iter < opts.maxiter

    if sign == "+"
        right = b - 1/x;
        I = IntegralIter(f, method, "a", a, "b", right, "tol", 0.1*opts.tol);

    else
        left = a + 1/x;
        I = IntegralIter(f, method, "a", left, "b", b, "tol", 0.1*opts.tol);
    end

    dif = abs(I - Iprev);
    Iprev = I;

    x = 2*x;
    iter = iter + 1;
end

% Optional outputs handling

if nargout < 3, clear iter; end
if nargout < 2, clear dif; end

end