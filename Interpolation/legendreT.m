function ln = legendreT(n)
% ln = legendreT(n)
% Returns coefficients for the n-th Truncated Legendre polynomial

l0 = 1;
ln = l0;
degree = 1;

while degree <= n
    if degree == 1
        l1 = [1, 0];
        ln = l1;
        degree = degree + 1;
    else
        if degree >= 2
            l2 = (2*degree-1)/degree * conv([1, 0], l1) - ...
                 (degree-1)/degree * [0, 0, l0];
            ln = l2;
            l0 = l1;
            l1 = l2;
            degree = degree + 1;
        end
    end
end

end