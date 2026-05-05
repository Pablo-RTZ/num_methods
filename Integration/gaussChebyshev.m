function I = gaussChebyshev(f,a,b,n)

%gaussChebyshev Gauss-Chebyshev quadrature for approximating definite integrals.
% Uses function handles and Chebyshev polynomials
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

xi = roots(Cheby(n));
xi = sort(xi);
ci=pi/n;
t = ((b - a)/2) * xi + (a + b)/2;
I = ((b - a)/2) * sum(ci .* f(t));

end