function [I,R,dif,iter]=Romberg(f,a,b,opts)

%Romberg Romberg's method for iteratively aproximating definite integrals.
% Uses function handles
%
%   I = Romberg(f,a,b)
%   Integrates in interval [a,b] using to a 1e-10 tolerance, max 100
%   iterations
%
%   I = Romberg(f,a,b, "tol", 1e-8,"maxiter",100)
%   allows name-value pair inputs in any order.
%
%   Inputs:
%       f       - function handle
%       a       - interval start
%       b       - interval end
%       tol     - difference between iterations
%       maxiter - maximum number of iterations allowed
%
%   Outputs:
%       I - Integral aproximation
%       R - Romberg matrix
%       dif - difference between iterations
%       iter - Number of iterations needed


arguments
    f (1,1) function_handle
    a (1,1) double
    b (1,1) double
    opts.tol (1,1) double = 1e-10
    opts.maxiter (1,1) int = 100
end

if ~(a < b)
    error('a must be less than b.');
end

% Initialization

dif=tol+1;
iter=1;
n=1;
k=1;
R(1,1)=Trapezoidal(f,a,b,n);

% Main loop

while and(dif>opts.tol, iter<opts.maxiter)
    k=k+1; n=2*n; iter=iter+1;
    R(k,1)=Trapezoidal(f,a,b,n);
    for j=2:k
        R(k,j)=(4^(j-1)*R(k,j-1)-R(k-1,j-1))/(4^(j-1)-1);
    end
    dif=abs(R(k,k)-R(k,k-1));
end

% Stopping criterion

if dif>tol
    disp("The method has not converged within the required tolerance")
else
    I=R(k,k);

end

% Optional outputs handling

if nargout < 4, clear iter; end
if nargout < 3, clear dif; end
if nargout < 2, clear R; end

end
