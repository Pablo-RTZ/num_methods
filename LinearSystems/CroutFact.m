function [sol,L,D] = CroutFact(A,b)

%CroutFact Solves Ax=b using Crout's factorization, and returns LDL^T factorization
% Works for complex matrices
%
%   sol = CroutFact(A,b)
%   Solves Ax=b using gaussian elimination
%
%   [sol, L,D] = CroutFact(A,b)
%   Solves Ax=b and returns LDL^T factorization
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

if not(rank(A) == n)
    error('Matrix is singular to working precission')
end

% Initialization

n = length(b);
L = zeros(n,n);
D = zeros(n,n);

L(1,1)=1;
D(1,1)=A(1,1);

% Main program

for i=2:n
    L(i,1)=A(i,1)/D(1,1);
    for j=2:i-1
        sum = 0;
        for k=1:j-1
            sum = sum + L(i,k)*D(k,k)*L(j,k);
        end
        L(i,j) = (A(i,j)-sum)/D(j,j);
    end

    L(i,i) = 1;
    
    sum = 0;
    for k = 1:i-1
        sum = sum + L(i,k)^2*D(k,k);
    end
    D(i,i) = A(i,i) - sum;
end

z = SD(L,b);
y = z./diag(D);
sol = SI(L.',y);

% Optional output handling

if nargout == 2, clear sol; end
if nargout < 3, clear L D; end

end