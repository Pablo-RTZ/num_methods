function [X,Y] = MilneSimpson(f, a,b, iv, n)

%MilneSimpson Milne's method for solving IVP, using Simpson as integrator.
%
%   [X, Y] = MilneSimpson(f,a,b, iv)
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

% RK4 as a predictor for the first iterations
for k=1:3
    k1 = feval(f,t(k),y(k,:))';
    k2 = feval(f,t(k)+h/2,y(k,:)+h/2*k1)';
    k3 = feval(f,t(k)+h/2,y(k,:)+h/2*k2)';
    k4 = feval(f,t(k+1),y(k,:)+h*k3)';
    y(k+1,:) = y(k,:)+h/6*(k1+2*k2+2*k3+k4);
end

% Main program

for k = 4:n
    % Milne as a predictor
    f_k = feval(f,t(k),y(k,:))';
    f_k_1 = feval(f,t(k-1),y(k-1,:))';
    f_k_2 = feval(f,t(k-2),y(k-2,:))';
    yp = y(k-3,:) + 4*h/3*(2*f_k - f_k_1 + 2*f_k_2);

    % Simpson as a corrector
    y(k+1,:) = y(k-1,:)+h/3*(feval(f,t(k+1),yp) + 4*feval(f,t(k),y(k,:)) + feval(f,t(k-1),y(k-1,:)))';
end

end