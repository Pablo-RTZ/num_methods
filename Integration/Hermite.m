function Hn = Hermite(n)
% Hn = Hermite(n)
% Returns coefficients for the n-th Hermite polynomial


if n == 0
    Hn = 1;
    return
elseif n == 1
    Hn = [1 0];
    return
end

H0 = 1;
H1 = [1 0];

for k = 2:n
    H2 = [H1, 0] - (k-1)*[0, 0, H0];
    H0 = H1;
    H1 = H2;
end

Hn = H2;
end