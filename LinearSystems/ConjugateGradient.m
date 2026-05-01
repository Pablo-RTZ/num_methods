function [sol,dif,res,iter,ACOC] = ConjugateGradient(A,b,opts)

%ConjugateGradient Iteratively solves Ax=b using conjugate gradient descent
% Works for real, positively defined matrices
%
%   sol = ConjugateGradient(A,b)
%   Solves Ax=b using gradient descent
%
%   Inputs:
%       A       - coefficient matrix
%       b       - independent term
%       x0      - initial guess
%       tol     - tolerance
%       maxiter - maximum number of iterations
%
%   Outputs:
%       sol     - aproximate solution of Ax=b
%       dif     - difference between iterations (vector)
%       res     - error at solution |F(x)| (vector)
%       iter    - number of iterations
%       ACOC    - aproximate computational order of convergence

arguments
    A (:,:) double
    b (:,1) double
    opts.x0 (:,1) double = zeros(len(b),1)
    opts.tol (1,1) double = 1e-8
    opts.maxiter (1,1) int = 50
end

if size(A,1) ~= size(A,2)
    error('The coefficient matrix must be square')
end

[~,p] = chol(A);
if p ~= 0
    error('The matrix is not positively defined')
end

% Initialization

x0 = x0(:);
b = b(:);
r0 = b-A*x0;
d0 = r0;
t0 = r0'*d0/(d0'*A*d0);

x1 = x0+t0*d0;
r1 = b-A*x1;
d1 = r1-(d0'*A*r1)/(d0'*A*d0)*d0;

res = norm(r1);
dif = norm(x1-x0);
iter = 1;

% Main program

while res(end)+dif(end)>tol && iter<=maxiter
    t = r1'*d1/(d1'*A*d1);
    x = x1+t*d1;
    r = b-A*x;
    d = r-(d1'*A*r)/(d1'*A*d1)*d1;

    iter = iter+1;
    res(iter) = norm(r);
    dif(iter) = norm(x-x1);
    x1=x;
    d1=d;
    r1=r;

end

ACOC = log(dif(3:end)./dif(2:end-1))./log(dif(2:end-1)./dif(1:end-2));

if res(end)+dif(end)>tol
    disp('The method has not converged')
else
    sol = x;
end

% Optional outputs handling

if nargout < 5, clear ACOC; end
if nargout < 4, clear iter; end
if nargout < 3, clear res; end
if nargout < 2, clear dif; end

end