function [sol,dif,res,iter,ACOC] = JOR(A,b,opts)

%JOR Iteratively solves Ax=b using Jacobi's Over Relaxation method
% Works for complex, strictly diagonal dominant matrices
%
%   sol = Jacobi(A,b)
%   Solves Ax=b using Jacobi's method
%
%   Inputs:
%       A       - coefficient matrix
%       b       - independent term
%       x0      - initial guess
%       tol     - tolerance
%       maxiter - maximum number of iterations
%       w       - parameter (convergence is ensured for 0<w<1)
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
    opts.w (1,1) double = 0.5
end

if size(A,1) ~= size(A,2)
    error('The coefficient matrix must be square')
end

D = abs(diag(A));
S = sum(abs(A), 2) - D;
if ~all(D > S)
    error('The matrix is not strictly diagonal dominant')
end

% Initialization

iter = 0;

res=tol;
dif=1;

L = tril(A,-1);
U = triu(A,1);
Dm1 = diag(1./diag(A));

% Main program

while res(end)+dif(end)>tol && iter<=maxiter
    x = ((1-w)*eye(n)-w*Dm1*(L+U))*x0+w*Dm1*b;
    iter = iter+1;
    res(iter) = norm(A*x-b);
    dif(iter) = norm(x-x0);
    x0=x;
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