function [sol,L,U] = CroutPenta(A,b)

%CroutPenta Solves Ax=b using Crout's method, and returns LU
% Works for complex, pentadiagonal matrices
%
%   sol = Crout(A,b)
%   Solves Ax=b using Crout's method extended to 5 diagonals
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

n = length(b);

dIi = diag(A,-2);
dI = diag(A,-1);
d = diag(A);
dS = diag(A,1);
dSs = diag(A,2);

% Initializes diagonals for L
p = zeros(n-2,1);
q = zeros(n-1,1);
r = zeros(n,1);
s = zeros(n-1,1);
t = zeros(n-2,1);

% Main program

% First rows
r(1) = d(1);
if n >= 2
    s(1) = dS(1) / r(1);
end
if n >= 3
    t(1) = dSs(1) / r(1);
end

if n >= 2
    q(2) = dI(1);
    r(2) = d(2) - q(2) * s(1);
    if n >= 3
        s(2) = (dS(2) - q(2) * t(1)) / r(2);
    end
    if n >= 4
        t(2) = dSs(2) / r(2);
    end
end

% Intermediate rows
for i = 3:n-2
    p(i) = dIi(i-2);
    q(i) = dI(i-1) - p(i) * s(i-2);
    r(i) = d(i) - p(i) * t(i-2) - q(i) * s(i-1);
    s(i) = (dS(i) - q(i) * t(i-1)) / r(i);
    t(i) = dSs(i) / r(i);
end

% Last rows
if n >= 3
    i = n-1;
    p(i) = dIi(i-2);
    q(i) = dI(i-1) - p(i) * s(i-2);
    r(i) = d(i) - p(i) * t(i-2) - q(i) * s(i-1);
    if n >= 2
        s(i) = (dS(i) - q(i) * t(i-1)) / r(i);
    end
end

if n >= 2
    i = n;
    p(i) = dIi(i-2);
    q(i) = dI(i-1) - p(i) * s(i-2);
    r(i) = d(i) - p(i) * t(i-2) - q(i) * s(i-1);
end

% Constructs LU
L = zeros(n,n);
U = eye(n);

for i = 1:n
    L(i,i) = r(i);
    if i >= 2
        L(i,i-1) = q(i);
    end
    if i >= 3
        L(i,i-2) = p(i);
    end
    
    if i <= n-1
        U(i,i+1) = s(i);
    end
    if i <= n-2
        U(i,i+2) = t(i);
    end
end

z = SD(L, b);
sol = SI(U, z);

% Optional arguments clearing

if nargout == 2, clear sol; end
if nargout < 3, clear L U; end

end