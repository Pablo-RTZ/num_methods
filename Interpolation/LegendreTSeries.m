function [LT, c] = LegendreTSeries(f, N, a, b)

%LegendreTSeries Compute Legendre series expansion of a function.
%   Returns a symbolic function representing the N-term Truncated Legendre polynomial
%   series expansion of f(x) over the interval [a, b].
%
%   [LT, c] = serieLegendreT(f, N, a, b)
%   Returns a symbolic function LT(x) that approximates f(x) on [a, b],
%   as well as the Fourier-Legendre coefficients.
%
%   Inputs:
%       f   - Symbolic function to be expanded (in terms of 'x')
%       N   - Number of terms in the Legendre series (non-negative integer)
%       a   - Left endpoint of the interval
%       b   - Right endpoint of the interval (b > a)
%
%   Outputs:
%       LT  - Symbolic function representing the Legendre series approximation
%       c   - Vector of Fourier-Legendre coefficients c_k for k = 0,...,N

arguments
    f sym
    N (1,1) int
    a (1,1) double
    b (1,1) double {mustBeGreaterThan(b, a)}
end

% Initialization

syms x
ft = subs(f, x, (b-a)/2 * x + (b+a)/2);
LT = sym(0);
c = zeros(1, N+1);

% Main program

for k = 0:N
    lk = poly2sym(legendreT(k), x);
    ck = (2*k+1)/2 * int(ft * lk, x, -1, 1);
    c(k+1) = ck;
    LT = LT + ck * lk;
end

LT = subs(LT, x, (2*x - b - a)/(b - a));

% Optional outputs handling

if nargout < 2, clear c; end

end