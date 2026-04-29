function F = TrigSeries(f, N, a, b)

%TrigSeries Compute trigonometric Fourier series of a function.
%   Returns a symbolic function representing the N-term Fourier series
%   expansion of f(x) over the interval [a, b].
%
%   F = serieTrigonometrica(f, N, a, b)
%   Returns a symbolic function F(x) that approximates f(x) on [a, b].
%
%   Inputs:
%       f   - Symbolic function to be expanded (in terms of 'x')
%       N   - Number of terms in the Fourier series (non-negative integer)
%       a   - Left endpoint of the interval
%       b   - Right endpoint of the interval (b > a)
%
%   Outputs:
%       F   - Symbolic function representing the Fourier series approximation

arguments
    f sym
    N (1,1) int
    a (1,1) double
    b (1,1) double {mustBeGreaterThan(b, a)}
end

syms x
ft = subs(f, x, (b-a)/(2*pi)*x + (a+b)/2);

for k = 0:N
    ak = 1/pi * int(ft * cos(k*x), x, -pi, pi);
    bk = 1/pi * int(ft * sin(k*x), x, -pi, pi);
    
    switch k
        case 0
            F = ak/2;
        case N
            F = F + ak * cos(N*x);
        otherwise
            F = F + ak * cos(k*x) + bk * sin(k*x);
    end
end

F = subs(F, x, (2*pi*x - pi*(a+b))/(b-a));

end