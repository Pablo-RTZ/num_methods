function Ln = Laguerre(n)
% Ln = Laguerre(n)
% Returns coefficients for the n-th Laguerre polynomial

if n == 0
    Ln = 1;
    return;
elseif n == 1
    Ln = [-1 1];
    return;
end

L0 = 1;
L1 = [-1 1];

for k = 2:n
    L2 = conv([-1, 2*(k-2) + 3], L1) - ((k-2) + 1)^2 * [0, 0, L0];
    L0 = L1;
    L1 = L2;
end

Ln = L2;

end