function [CT, c] = ChebyshevTSeries(f, N, a, b)

%ChebyshevTSeries Compute Chebyshev series expansion of a function.
%   Returns a symbolic function representing the N-term Truncated Chebyshev polynomial
%   series expansion of f(x) over the interval [a, b].
%
%   [CT, c] = serieChebyshevT(f, N, a, b)
%   Returns a symbolic function CT(x) that approximates f(x) on [a, b],
%   as well as the Chebyshev coefficients.
%
%   Inputs:
%       f   - Symbolic function to be expanded (in terms of 'x')
%       N   - Number of terms in the Chebyshev series (non-negative integer)
%       a   - Left endpoint of the interval
%       b   - Right endpoint of the interval (b > a)
%
%   Outputs:
%       CT  - Symbolic function representing the Chebyshev series approximation
%       c   - Vector of Chebyshev coefficients c_k for k = 0,...,N

arguments
    f sym
    N (1,1) int
    a (1,1) double
    b (1,1) double {mustBeGreaterThan(b, a)}
end

% Initialization

syms x
ft = subs(f, x, (b-a)/2 * x + (b+a)/2);
CT = sym(0);
c = zeros(1, N+1);

% Main program

for k = 0:N
    tk = poly2sym(chebyT(k), x);
    ck = 2/pi * int((ft * tk) / sqrt(1 - x^2), x, -1, 1);
    c(k+1) = ck;
    
    if k == 0
        CT = CT + ck/2;
    else
        CT = CT + ck * tk;
    end
end

CT = subs(CT, x, (2*x - b - a)/(b - a));

% Optional outputs handling

if nargout < 2, clear c; end

end