function [I, dif, iter] = IntegralInf(f,method, a,sign,opts)

%IntegralInf Applies an iterative scheme to aproximate infinite integrals.
% Method must have arguments (f,a,b,n)
%
%   I = IntegralInf(f,method, 0,+,)
%   Integrates using the method given in the interval [0,+Inf]
%
%   I = IntegralIter(f,method,0,-,"tol",1e-8,"maxiter",200)
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
    method (1,1) str
    a (1,1) double
    sign (1,1) str {mustBeMember(sign,{"+","-"})} = "+"
    opts.tol (1,1) double = 1e-10
    opts.maxiter (1,1) int = 100
end

filename = method + ".m";
if ~isfile(filename)
    error("Method not found on current directory");
end

% Initialization

iter = 1;
dif = tol +1;
x = a;
I = 0;

while abs(dif) < tol && iter < maxiter
    if sign == "+"
        dif = IntegralIter(f,method,"a",x,"b",2*x,"tol",0.1*tol);
    else
        dif = IntegralIter(f,method,"a",2*x,"b",x,"tol",0.1*tol);
    end
    x = 2*x;
    iter = iter+1;
    I = I + dif;
end

% Optional outputs handling

if nargout < 3, clear iter; end
if nargout < 2, clear dif; end

end