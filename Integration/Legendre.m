function Pn = Legendre(n)
% Pn = Legendre(n)
% Returns coefficients for the n-th Legengre polynomial

if n == 0
    Pn = 1;
    return;
elseif n == 1
    Pn = [1 0];
    return;
end

Pnm1 = 1;
Pn = [1 0];

for k = 2:n
    Pn_x = conv([1 0], Pn);
    term1 = (2*k - 1) * Pn_x;
    len_diff = length(term1) - length(Pnm1);
    if len_diff > 0
        Pnm1 = [zeros(1, len_diff), Pnm1];
    elseif len_diff < 0
        term1 = [zeros(1, -len_diff), term1];
    end
    term2 = (k - 1) * Pnm1;

    Pk = (term1 - term2) / k;

    Pnm1 = Pn;
    Pn = Pk;
end
end
