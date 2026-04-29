function [I, dif, iter] = IntegralImp(f,method, a, b, sign,opts)

%IntegralImp Applies an iterative scheme to aproximate improper integrals.
% Method must have arguments (f,a,b,n)
%
%   I = IntegralInf(f,method, 0,1,+)
%   Integrates using the method given in the interval [0,1] with
%   discontinuity at 1
%
%   I = IntegralInf(f,method,0,1,-,"tol",1e-8,"maxiter",200)
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
    method (1,1) str
    a (1,1) double
    b (1,1) double
    sign (1,1) str {mustBeMember(sign,{"+","-"})} = "+"
    opts.tol (1,1) double = 1e-10
    opts.maxiter (1,1) int = 100
end

filename = method + ".m";
if ~isfile(filename)
    error("Method not found on current directory");
end

if ~(a < b)
    error('a must be less than b.');
end

% Initialization

iter = 1;
dif = tol +1;
Iprev = 0;

while abs(dif) < tol && iter < maxiter
    if sign == "+"
        I = IntegralIter(f,method,"a",a,"b",(1-1./iter)*b,"tol",0.1*tol);
    else
        I = IntegralIter(f,method,"a",(1-1./iter)*a,"b",b,"tol",0.1*tol);
    end
    iter = iter+1;
    dif = abs(I-Iprev);
    Iprev=I;
end

% Optional outputs handling

if nargout < 3, clear iter; end
if nargout < 2, clear dif; end

end