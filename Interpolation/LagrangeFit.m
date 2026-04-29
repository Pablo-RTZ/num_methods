function p = LagrangeFit(xi, fi)

%LagrangeFit Lagrange polynomial interpolation.
%   Returns a symbolic polynomial for Lagrange interpolation
%
%   p = LagrangeFit(xi, fi)
%   Returns a symbolic function p(x) that interpolates the given points
%
%   Inputs:
%       xi  - independent variable data points
%       fi  - dependent variable data points
%
%   Outputs:
%       p   - Symbolic Lagrange polynomial function

arguments
    xi (:,1) double
    fi (:,1) double
end

if length(xi) ~= length(fi)
    error('xi and fi must have the same length.');
end

if length(xi) < 2
    error('At least 2 data points are required for Lagrange interpolation.');
end

% Initialization

xi = xi(:);
fi = fi(:);
n = length(xi);
syms x

L = [];

Num0 = x - xi;

% Main program

for i = 1:n
    Num = Num0;
    Num(i) = [];
    
    Den = xi(i) - xi;
    Den(i) = [];

    L = [L prod(Num) / prod(Den)];
end

p = sum(L(:) .* fi);

end