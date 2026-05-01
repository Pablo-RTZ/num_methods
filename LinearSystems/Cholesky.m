function [sol, L] = Cholesky(A,b)

%Cholesky Solves Ax=b using cholesky factorization, and returns LL^T factorization
% Works for real, positively defined matrices
%
%   sol = Cholesky(A,b)
%   Solves Ax=b using Cholesky factorization
%
%   [sol, L] = Gauss(A,b)
%   Solves Ax=b and returns LL^T factorization
%
%   Inputs:
%       A   - coefficient matrix
%       b   - independent term

arguments
    A (:,:) double
    b (:,1) double
end

if size(A,1) ~= size(A,2)
    error('The coefficient matrix must be square')
end

[~,p] = chol(A);
if p ~= 0
    error('The matrix is not positively defined')
end

n = length(b);
L = zeros(n,n);
L(1,1) = sqrt(A(1,1));

% Main program

for i=2:n
    L(i,1) = A(i,1)/L(1,1);
end

for i = 2:n-1
    L(i,i) = sqrt(A(i,i)-sum(L(i,1:i-1).^2));
    for j = i+1:n
        L(j,i) = (A(i,j)-sum(L(i,1:i-1).*L(j,1:i-1)))./L(i,i);
    end
end

L(n,n) = sqrt(A(n,n)-sum(L(n,1:n-1).^2));
z = SD(L,b);
sol = SI(L',z);

if nargout <2, clear L; end

end