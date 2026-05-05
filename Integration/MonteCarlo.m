function I = MonteCarlo(f,a,b,n)

%MonteCarlo Monte-Carlo method for aproximating definite integrals.
% Uses function handles
%
%   I = MonteCarlo(f,a,b)
%   Integrates in interval [a,b] using 500 points by default
%
%   I = MonteCarlo(f,a,b, 100)
%   changes number of points
%
%   Inputs:
%       f   - function handle
%       a   - interval start
%       b   - interval end
%       n   - number of points

arguments
    f (1,1) function_handle
    a (1,1) double
    b (1,1) double
    n (1,1) double {mustBeInteger, mustBeNonnegative} = 500
end

if ~(a < b)
    error('a must be less than b.');
end

% Main program

nodes = a+(b-a)*rand(1,n);
I = (b-a)/n*sum(f(nodes));

end