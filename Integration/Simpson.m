function I = Simpson(f,a,b,n)

%Simpson Compound simpson's method for aproximating definite integrals.
% Uses function handles
%
%   I = Simpson(f,a,b)
%   Integrates in interval [a,b] using 10 points by default
%
%   I = Simpson(f,a,b, 100)
%   changes interval number
%
%   Inputs:
%       f   - function handle
%       a   - interval start
%       b   - interval end
%       n   - number of intervals

arguments
    f (1,1) function_handle
    a (1,1) double
    b (1,1) double
    n (1,1) double {mustBeInteger, mustBeNonnegative} = 10
end

if ~(a < b)
    error('a must be less than b.');
end

% Main program

h = (b-a)/n;
nodes = a:h:b;
weights = 2*ones(1,n+1);
weights(1)=1;
weights(end)=1;
weights(2:2:end-1)=4;
I = h/3*sum(weights.*f(nodes));

end