function [sol,Q,R] = QR(A,b)

%QR Solves Ax=b using QR factorization, and returns QR
% Works for complex matrices, if columns are LI
%
%   sol = QR(A,b)
%   Solves Ax=b using QR factorization
%
%   [sol, Q,R] = QR(A,b)
%   Solves Ax=b and returns QR factorization
%
%   Inputs:
%       A   - coefficient matrix
%       b   - independent term

arguments
    A (:,:) double
    b (:,1) double
end

if rank(A) ~= size(A, 2)
    error('Columns arent lineally independent')
end

b = b(:);
[m,n] = size(A);
Q = zeros(m,n);
R = zeros(n,n);

% Main program

for i = 1:n
    v = A(:,i);

    for j = 1:i-1
        R(j,i) = dot(Q(:,j), A(:,i));
        v = v - R(j,i) * Q(:,j);
    end
    R(i,i) = norm(v);
    Q(:,i) = v / R(i,i);
end

z = Q'*b;
sol = SI(R,z);

end