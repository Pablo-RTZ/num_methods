function [X,Y] = RungeKutta(f, a,b, iv, n)

%RungeKutta RungeKutta's method (RK4) for solving IVP.
%
%   [X, Y] = RungeKutta(f,domain, iv)
%   Solves the problem and returns function at evaluation nodes
%
%   Inputs:
%       f   - ODE
%       a   - Interval start
%       b   - Interval end
%       iv  - Initial value
%       n   - number of nodes to interpolate

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
    k1 = feval(f,t(k),y(k,:))';
    k2 = feval(f,t(k)+h/2,y(k,:)+h/2*k1)';
    k3 = feval(f,t(k)+h/2,y(k,:)+h/2*k2)';
    k4 = feval(f,t(k+1),y(k,:)+h*k3)';
    y(k+1,:) = y(k,:)+h/6*(k1+2*k2+2*k3+k4);
end

end