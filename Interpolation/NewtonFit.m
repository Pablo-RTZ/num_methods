function p = NewtonFit(xi, fi)

%NewtonFit Newton polynomial interpolation.
%   Returns a symbolic polynomial using Newton's divided differences
%
%   p = NewtonFit(xi, fi)
%   Returns a symbolic function p(x) that interpolates the given points
%
%   Inputs:
%       xi  - independent variable data points
%       fi  - dependent variable data points
%
%   Outputs:
%       p   - Symbolic Newton polynomial function

arguments
    xi (:,1) double
    fi (:,1) double
end

if length(xi) ~= length(fi)
    error('xi and fi must have the same length.');
end

if length(xi) < 2
    error('At least 2 data points are required for Newton interpolation.');
end

% Initialization

syms x
xi = xi(:);
fi = fi(:);
n = length(xi);

vx = sym(ones(1, n));
D = fi;

% Main program

for col = 2:n
    for fil = col:n
        D(fil, col) = (D(fil, col-1) - D(fil-1, col-1)) / ...
                      (xi(fil) - xi(fil-col+1));
    end
    vx(col) = vx(col-1) * (x - xi(col-1));
end


dD = diag(D);
p = sum(dD .* vx(:));

end