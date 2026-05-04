function [X,Y] = Euler(f, a,b, iv, n)

%Euler Euler's method for solving IVP.
%
%   [X, Y] = Euler(f,a,b, iv)
%   Solves the problem and returns function at evaluation nodes
%
%   Inputs:
%       f   - ODE
%       a   - Interval start
%       b   - Interval end
%       iv  - Initial value
%       n   - Number of nodes to interpolate

arguments
    f (:,1)
    a (1,1) double
    b (1,1) double {mustBeGreaterThan(b,a)}
    iv (:,1) double
    n (1,1) int = 10
end

% Initialization

h = (b-a)/n;
x = a:h:b;
y = zeros(n+1,length(iv));
y(1,:) = iv;

% Main program

for k = 1:n
    y(k+1,:) = y(k,:)+h.*feval(f,x(k),y(k,:))';
end

end