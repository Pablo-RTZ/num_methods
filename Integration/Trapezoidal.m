function I = Trapezoidal(f,a,b,n)

%Trapezoidal Compound trapezoidal method for aproximating definite integrals.
% Uses function handles
%
%   I = Trapezoidal(f,a,b)
%   Integrates in interval [a,b] using 50 points by default
%
%   I = Trapezoidal(f,a,b, 100)
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
    n (1,1) int = 50
end

if ~(a < b)
    error('a must be less than b.');
end

% Main program

h = (b-a)/n;
nodes = a:h:b;
weights = 2*ones(1,n+1);
weights(1)=1;
weights(n+1)=1;
I = h/2*sum(weights.*f(nodes));

end