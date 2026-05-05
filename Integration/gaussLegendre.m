function I = gaussLegendre(f,a,b,n)

%gaussLegendre Gauss-Legendre quadrature for approximating definite integrals.
% Uses function handles and Legendre polynomials
%
%   I = gaussLegendre(f,a,b)
%   Integrates in interval [a,b] using 5 points by default
%
%   I = gaussLegendre(f,a,b, 10)
%   changes number of quadrature points
%
%   Inputs:
%       f   - function handle
%       a   - interval start
%       b   - interval end
%       n   - number of quadrature points

arguments
    f (1,1) function_handle
    a (1,1) double
    b (1,1) double
    n (1,1) double {mustBeInteger, mustBeNonnegative} = 5
end

if ~(a < b)
    error('a must be less than b.');
end

% Main program

xi = roots(Legendre(n));
xi = sort(xi);
syms x
polinomio = poly2sym(Legendre(n));
dp = diff(polinomio, x);
p = subs(dp, x, xi);
ci = 2./((1-xi.^2).*(p.^2));
t = ((b - a)/2) * xi + (a + b)/2;
I = (b - a)./2 .* sum(ci .* f(t));

end