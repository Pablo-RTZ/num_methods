function [sol, dif, iter, ACOC] = Bisection(f, a, b, opts)

%Bisection Bisection method for solving nonlinear equations.
% Uses function handles, and uses function error as stopping criterion.
%
%   sol = Bisection(f, a, b)
%   Uses default tolerance 1e-10 and 50 iterations
%
%   sol = Bisection(f, a, b, "tol", 1e-8, "maxiter", 100)
%   allows name-value pair inputs in any order.
%
%   Outputs:
%       sol   - root approximation
%       dif   - last step size
%       iter  - number of iterations
%       ACOC  - Approximate Computational Order of Convergence

arguments
    f (1,1) function_handle
    a (1,1) double
    b (1,1) double

    opts.tol (1,1) double = 1e-10
    opts.maxiter (1,1) double = 50
end

% Initialization

tol = opts.tol;
maxiter = opts.maxiter;

iter = 0;
dif = tol + 1;
I = [];

% Check initial interval

if f(a)*f(b) > 0
    sol = [];
    ACOC = [];
    disp("Interval does not have opposite signs at the endpoints")
    return
elseif f(a)*f(b) == 0
    if f(a) == 0
        sol = a;
        dif = 0;
        iter = 1;
        ACOC = [];
        disp("Root found at left endpoint")
        return
    elseif f(b) == 0
        sol = b;
        dif = 0;
        iter = 1;
        ACOC = [];
        disp("Root found at right endpoint")
        return
    end
end

% Main loop

x_left = a;
x_right = b;
x_mid = x_left;
x_mid_old = x_left;

while (abs(f(x_mid)) > tol) && (iter < maxiter)
    x_mid = (x_left + x_right) / 2;
    dif = abs(x_mid_old - x_mid);
    I(end+1) = dif;
    
    if f(x_left)*f(x_mid) < 0
        x_right = x_mid;
    else
        x_left = x_mid;
    end
    
    x_mid_old = x_mid;

    iter = iter + 1;
end

% Stopping criterion

if iter >= maxiter
    sol = [];
    ACOC = [];
    disp("The method has not converged to a root")
else
    sol = x_mid;
    
    if numel(I) >= 3
        ACOC = log(I(3:end)./I(2:end-1)) ./ log(I(2:end-1)./I(1:end-2));
    else
        ACOC = [];
    end
end

% Optional outputs handling

if nargout < 4, clear ACOC; end
if nargout < 3, clear iter; end
if nargout < 2, clear dif; end

end