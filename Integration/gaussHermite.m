function I = gaussHermite(f,n)

%gaussHermite Gauss-Hermite quadrature for approximating indefinite integrals.
% Uses function handles and Hermite polynomials
%
%   I = gaussHermite(f)
%   Integrates in (-Inf, Inf) using 5 points by default
%
%   I = gaussHermite(f, 10)
%   changes number of quadrature points
%
%   Inputs:
%       f   - function handle
%       n   - number of quadrature points

arguments
    f (1,1) function_handle
    n (1,1) double {mustBeInteger, mustBeNonnegative} = 5
end

% Main program

xi = roots(Hermite(n));

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
    L = L.*exp(-x.^2);
    L = matlabFunction(L);
    c(i) = integral(L,-Inf,Inf);

end

I = sum(c'.*f(xi));
end