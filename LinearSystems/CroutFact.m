function [sol,L,D] = CroutFact(A,b)

%CroutFact Solves Ax=b using Crout's factorization, and returns LDL^T factorization
% Works for complex, symmetric matrices.
%
%   sol = CroutFact(A,b)
%   Solves Ax=b using LDL^T factorization
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

n = length(b);

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

for j = 1:n
    % Compute D(j,j)
    sum_k = 0;
    for k = 1:j-1
        sum_k = sum_k + L(j,k)^2 * D(k,k);
    end
    D(j,j) = A(j,j) - sum_k;

    % Set diagonal of L to 1
    L(j,j) = 1;

    % Compute L(i,j) for i > j
    for i = j+1:n
        sum_k = 0;
        for k = 1:j-1
            sum_k = sum_k + L(i,k) * D(k,k) * L(j,k);
        end
        L(i,j) = (A(i,j) - sum_k) / D(j,j);
    end
end

z = DS(L,b);
y = z./diag(D);
sol = IS(L.',y);


% Optional output handling

if nargout == 2, clear sol; end
if nargout < 3, clear L D; end

end