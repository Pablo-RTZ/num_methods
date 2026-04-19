function I = gaussLobatto(f,a,b,n)
%gaussLobatto Gauss-Lobatto quadrature for approximating definite integrals.
% Uses function handles and Legendre polynomials, but using endpoints as
% nodes
%
%   I = gaussLobatto(f,a,b)
%   Integrates in interval [a,b] using 5 points by default
%
%   I = gaussLobatto(f,a,b, 10)
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
    n (1,1) int = 5
end

if ~(a < b)
    error('a must be less than b.');
end

% Main program

xi=roots(polyder(Legendre(n-1)));
xi=[-1; xi; 1];
syms x
polinomio=poly2sym(Legendre(n-1));
p=subs(polinomio,x,xi);
ci = 2./(n*(n-1).*p.^2);
t = ((b - a)/2) * xi + (a + b)/2;
y=f(t);
I=(b-a)./2.*sum(ci.*y);

end