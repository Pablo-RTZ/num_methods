function [X,Y] = MilneTrap(f, a,b, iv, n)

%MilneTrap Milne's method for solving IVP, using Trapezoidal as integrator.
%
%   [X, Y] = MilneTrap(f,a,b, iv)
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

% Heun as a predictor for the first iterations
for k=1:3
    k1 = feval(f,t(k),y(k,:))';
    k2 = feval(f,t(k+1),y(k,:)+h*k1)';
    y(k+1,:) = y(k,:)+h/2*(k1+k2);
end

% Main program

for k = 4:n
    % Milne as a predictor
    yp = y(k-3,:) + 4*h/3 * (2*feval(f,t(k),y(k,:)) - feval(f,t(k-1),y(k-1,:)) + 2*feval(f,t(k-2),y(k-2,:)));
    
    % Implicit trapezoidal as a corrector
    y(k+1,:) = y(k,:)+h/2*(feval(f,t(k),y(k,:)) + feval(f,t(k+1),yp))';
end

end