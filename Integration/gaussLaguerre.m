function I = gaussLaguerre(f,n)

%gaussLaguerre Gauss-Laguerre quadrature for approximating indefinite integrals.
% Uses function handles and Laguerre polynomials
%
%   I = gaussLaguerre(f)
%   Integrates in [0, Inf) using 5 points by default
%
%   I = gaussLaguerre(f, 10)
%   changes number of quadrature points
%
%   Inputs:
%       f   - function handle
%       n   - number of quadrature points

arguments
    f (1,1) function_handle
    n (1,1) int = 5
end

% Main program

xi = roots(Laguerre(n));

syms x
Num0 = x - xi;
L = [];

% Manually calculates weights (no closed formula)

for i = 1:n
    Num = Num0;
    Num(i) = [];

    Den = xi(i)-xi;
    Den(i)=[];

    L = prod(Num)./prod(Den);
    L = L.*exp(-x);
    L = matlabFunction(L);
    c(i) = integral(L,0,Inf);

end

I = sum(c'.*f(xi));
end