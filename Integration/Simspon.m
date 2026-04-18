function I = Simspon(f,a,b,n)

%Simpson Compound simpson method for aproximating definite integrals.
% Uses function handles
%
%   I = Simpson(f,a,b)
%   Integrates in interval [a,b] using 50 points by default
%
%   I = Simpson(f,a,b, "n", 100)
%   allows name-value pair inputs in any order.
%
%   Outputs:
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
weights(end)=1;
weights(2:2:end-1)=4;
I = h/3*sum(weights.*f(nodes));

end