function [sol,dif,res,iter,ACOC] = SOR2(A,b,opts)

%SOR2 Iteratively solves Ax=b using another formulation for the SOR method
% Works for complex matrices
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
%       w       - parameter
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

% Initialization

iter = 0;

res=tol;
dif=1;

L = tril(A,-1);
U = triu(A,1);
D = diag(diag(A));

M=U*w+D;

% Main program

while res(end)+dif(end)>tol && iter<=maxiter
    x = SI(M,((1-w)*D-w*L)*x0+b*w);

    iter=iter+1;
    res(iter)=norm(A*x-b,2);
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