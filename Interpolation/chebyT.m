function tTn = chebyT(n)
% tTn = chebyT(n)
% Returns coefficients for the n-th Truncated Chebyshev polynomial

tT0 = 1;
tTn = tT0;
degree = 1;

while degree <= n
    if degree == 1
        tT1 = [1, 0];
        tTn = tT1;
        degree = degree + 1;
    else
        if degree >= 2
            tT2 = 2 * conv([1, 0], tT1) - [0, 0, tT0];
            tTn = tT2;
            tT0 = tT1;
            tT1 = tT2;
            degree = degree + 1;
        end
    end
end

end