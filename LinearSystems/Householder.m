function [sol,Q,R] = Householder(A,b)

%Householder Solves Ax=b using Householder transformations, and returns QR factorization
% Works for complex matrices, rectangular (if not overrestricted)
%
%   sol = Householder(A,b)
%   Solves Ax=b using QR factorization
%
%   [sol, Q,R] = Householder(A,b)
%   Solves Ax=b and returns QR factorization
%
%   Inputs:
%       A   - coefficient matrix
%       b   - independent term

arguments
    A (:,:) double
    b (:,1) double
end

[m, n] = size(A);
Q = eye(m);
A1 = A;

% Main program

for i = 1:min(m-1,n)
    % Column to triangularize
    v = A1(i:m, i);

    % Householder vector
    v(1) = v(1) + norm(v);

    % Householder matrix
    M = eye(m-i+1) - (2/norm(v)^2) * (v*v');
    H = [eye(i-1), zeros(i-1, m-i+1); zeros(m-i+1, i-1), M];

    % Applies transform
    A1 = H * A1;
    Q = Q * H;
end

R=A1;

z = Q'*b;
sol = SI(R,z);

% Optional output handling

if nargout == 2, clear sol; end
if nargout < 3, clear Q R; end

end