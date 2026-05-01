function sol = Crout(A,b)

%Crout Solves Ax=b using Crout's method
% Works for complex, tridiagonal matrices
%
%   sol = Crout(A,b)
%   Solves Ax=b using Crout's method
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

% Initialization

% Extracts diagonals
dS = diag(A,1);
dP = diag(A);
dI = diag(A,-1);

l(1)=dP(1);
u(1)=dS(1)/l(1);

% Main program

for i=2:n-1
    l(i)=dP(i)-dI(i-1)*u(i-1);
    u(i)=dS(i)/l(i);
end
l(n)=dP(n)-dI(n-1)*u(n-1);

% Solves Lz = d (using DS)
z(1)=b(1)/l(1);
for i=2:n
    z(i)=(1/l(i))*(b(i)-dI(i-1)*z(i-1));
end

% Solves Ux = z (using IS)
sol(n)=z(n);
for i=n-1:-1:1
    sol(i)=z(i)-u(i)*sol(i+1);
end

sol=sol(:);

end