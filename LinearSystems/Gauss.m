function [sol, L,U] = Gauss(A,b)

%Gauss Solves Ax=b using gaussian elimination, and returns LU factorization
% Works for complex matrices
%
%   sol = Gauss(A,b)
%   Solves Ax=b using gaussian elimination
%
%   [sol, L,U] = Gauss(A,b)
%   Solves Ax=b and returns LU factorization
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

[~,n] = size(A);
for i=1:n
    if det(A(1:i,1:i)) == 0
        error('The matrix has zero pivots')
    end
end

L = eye(n);
U = A;

% Main program

for k = 1:n-1
    for i = k+1:n
        L(i,k) = U(i,k) / U(k,k);
        for j = k:n
            U(i,j) = U(i,j) - L(i,k) * U(k,j);
        end
    end
end

z = SD(L,b);
sol = SI(U,z);

% Optional output handling

if nargout == 2, clear sol; end
if nargout < 3, clear L U; end

end